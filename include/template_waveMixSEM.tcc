

//constructor 
template<int dim>
   waveMixSEM<dim>::waveMixSEM(const Parameters &prm):
    cg_fe(dii::FE_Q<dim>(prm.p_degree) ^ dim, dii::FE_Nothing<dim>()^dim),
    dg_fe(dii::FE_Nothing<dim>()^dim, dii::FE_DGQ<dim>(prm.p_degree) ^ dim),
    dof_handler(mesh),
    param(prm),
    stiffness_tensor(get_stiffness_tensor<dim>(param.vp,param.vs,param.rho)),
    source_function(param),
    timestep_number(1),
    mesh_op(prm)
    
   { 
    fe_collection.push_back(cg_fe);
    fe_collection.push_back(dg_fe);

    qcollection.push_back(dii::QGaussLobatto<dim>(param.p_degree + 1));
    qcollection.push_back(dii::QGaussLobatto<dim>(param.p_degree + 1));

    face_qcollection.push_back(dii::QGaussLobatto<dim - 1>(param.p_degree + 1));
    face_qcollection.push_back(dii::QGaussLobatto<dim - 1>(param.p_degree + 1));
    out.open ("waveMixSEM_log.txt");

    std::cout << " Source position: " << source_function.source_pos << std::endl;
    std::cout<< "constructor called" <<std::endl;

    out << " Source position: " << source_function.source_pos << std::endl;
   }
   

template<int dim>
  void waveMixSEM<dim>::make_mesh()
  {
    //creating a uniform mesh  
    mesh_op.mesh_create(mesh);
    const double h_bg = mesh.begin_active()->diameter();
    double h = mesh.begin_active()->extent_in_direction(1);
    std::cout <<"Lx: " << param.Lx << " Ly: " << param.Ly << std::endl;

    std::vector<dii::GridTools::PeriodicFacePair<typename dii::Triangulation<dim>::cell_iterator>>
                      tria_periodic_faces;
    
    if(param.source == 4)
    {
      dii::GridTools::collect_periodic_faces(mesh,
                                            boundary_id_left,
                                            boundary_id_right,
                                            direction,
                                            tria_periodic_faces);
    
      mesh.add_periodicity(tria_periodic_faces);
    }

    // Get the fracture coordinate in physical mesh units
  
    double frac_coord=mesh_op.get_fcoord();
    std::cout << "Suggested fracture coordinate: " << frac_coord
                << std::endl;
   //Identifiying cg and dg regions (material_id and settting fe type)
   //cg_id-> fe_index(0) and dg1_id and dg2_id -> fe_index(1)
   const int rows_above=std::max(1, param.dg_rows/2);
   double top_dg_y=0.0, bottom_dg_y=0.0;
   cg_dg_cells<dim>(dof_handler, frac_coord, rows_above,
                     cg_id, dg1_id, dg2_id,
                     top_dg_y, bottom_dg_y);

    // Refining DG rows (and 2 CG rows above/below) if requested

    mesh_op.print_mesh(mesh);


   // Obtaining receivers' coordinates
    const unsigned int noffset= 1;
    receivers_pos = mesh_op.receiver_coordinates(h, frac_coord, bottom_dg_y,
                                                 noffset, param.receivers_x);


    std::cout << "final fracture coordinate at dg interface: " << frac_coord
                << std::endl << std::endl;
   
    std::cout << "Receivers' coordinates:" << std::endl;
                for (unsigned int i = 0; i < receivers_pos.size(); ++i)
                  std::cout << "  [" << i << "] " << receivers_pos[i] << std::endl;

    std::cout << std::endl;
    std::cout << "Number of active cells: " << mesh.n_active_cells()
    << std::endl;


    // Defining recevers' coordinates

    out << "Mesh Lx: " << mesh_op.Lx <<std::endl;
    out << "Mesh Ly: " << mesh_op.Ly <<std::endl;
    out << "fracture  coordinate: " << frac_coord <<std::endl;
    out << "top_dg_y: " << top_dg_y <<std::endl;
    out << "bottom_dg_y: " << bottom_dg_y <<std::endl;
    out << "Number of active cells: " << mesh.n_active_cells()<< std::endl << std::endl;
    out << "background cell diameter: " << h_bg << std::endl;
    out << "Receivers' coordinates:" << std::endl;
                for (unsigned int i = 0; i < receivers_pos.size(); ++i)
                  out << "  [" << i << "] " << receivers_pos[i] << std::endl;

    out << std::endl;
  }
  

template<int dim>
  void waveMixSEM<dim>::setup_system()

  {    
     dof_handler.distribute_dofs(fe_collection);

  
    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
                << std::endl
                << std::endl;
    
     //periodic Boundary conditions


    dii:: GridTools::collect_periodic_faces(dof_handler,
                                            boundary_id_left,
                                            boundary_id_right,
                                            direction,   // 0 for x, 1 for y
                                            periodic_faces);
    //applying periodic boundary conditions
    constraints.clear(); 

    if(param.source == 4)
        dii::DoFTools::make_periodicity_constraints<dim,dim>(periodic_faces,
                                                             constraints); 
    constraints.close();  
  
    //Flux sparsity pattern
    dii::DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    dii::DoFTools::make_flux_sparsity_pattern(dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from(dsp);

    //Volume sparsity pattern
    
    dii::DynamicSparsityPattern dsp_volume(dof_handler.n_dofs(), dof_handler.n_dofs());
    dii::DoFTools::make_sparsity_pattern(dof_handler, dsp_volume, constraints, false);
    sparsity_pattern_volume.copy_from(dsp_volume);
    
    //initializing global matrices

    stiffness_matrix.reinit(sparsity_pattern);
    boundary_matrix.reinit(sparsity_pattern_volume);

    // initializing  global vectors
    mass_vector.reinit(dof_handler.n_dofs());
    inverse_mass_vector.reinit(dof_handler.n_dofs());
    solution_u.reinit(dof_handler.n_dofs());
    solution_v.reinit(dof_handler.n_dofs());
    old_solution_u.reinit(dof_handler.n_dofs());
    old_solution_v.reinit(dof_handler.n_dofs());
    source_term.reinit(dof_handler.n_dofs());
   // writing to file

  out << "Number of degrees of freedom: " << dof_handler.n_dofs()
   << std::endl
   << std::endl;

  }

template<int dim>
  void waveMixSEM<dim>::assemble_MKmatrices()
  {
    std::cout << " assembly begins " << std::endl;
    out << " assembly begins " << std::endl;
    
    
    //Assembly of  mass and stiffness matrices
    /**************************** */
    dii::hp::FEValues<dim> hp_fe_values(fe_collection, qcollection,
                                  dii::update_values | dii::update_gradients |
                                  dii::update_quadrature_points | dii::update_JxW_values);

    dii::FullMatrix<double> cell_stiffness_matrix;
    std::vector<dii::types::global_dof_index> local_dof_indices ;
    dii::Vector<double> cell_mass_vector;
    //assembly per cell

    for (const auto &cell : dof_handler.active_cell_iterators())
        {
        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
        cell_stiffness_matrix.reinit(dofs_per_cell,dofs_per_cell);
        cell_stiffness_matrix =0.0;
        cell_mass_vector.reinit(dofs_per_cell);
        cell_mass_vector =0.0;

        hp_fe_values.reinit(cell);
        const dii::FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();
       
        for (const unsigned int q_point : fe_values.quadrature_point_indices())
        {
          //Precompute strain
          std::vector<dii::SymmetricTensor<2, dim>> eps_phi(dofs_per_cell);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            eps_phi[i] = get_strain(fe_values, i, q_point);

          for (const unsigned int i : fe_values.dof_indices())
          { 
            double phi_i = fe_values.shape_value(i,q_point);
            cell_mass_vector(i) +=  phi_i * phi_i *fe_values.JxW(q_point);

            for (const unsigned int j : fe_values.dof_indices())                 
                cell_stiffness_matrix(i, j) += (eps_phi[i] *  stiffness_tensor *  eps_phi[j]) *                    
                                               fe_values.JxW(q_point);        
          }
        }

        //Trasfer local cell matrix to global matrix
        local_dof_indices.resize(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);

        constraints.distribute_local_to_global(cell_stiffness_matrix, cell_mass_vector,
                     local_dof_indices, stiffness_matrix, mass_vector);

        }//end of loop by cell dof   
  
    // inverse mass vector
    const double mass_eps =
    1.0e-20 * std::max(mass_vector.linfty_norm(), 1.0);
  
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
       inverse_mass_vector(i) =
       (mass_vector(i) > mass_eps) ? (1.0 / mass_vector(i)) : 0.0;
  
    constraints.distribute(inverse_mass_vector); 


  }

    template<int dim>
  void waveMixSEM<dim>::include_interface_fluxes()
  { 
    std::cout << " flux assembly begins " << std::endl;
    //Assembly of interface flux terms into stiffness matrix using MeshWorker
    /**************************** */
    const dii::SymmetricTensor<2, dim> fracture_stiffness = 
                                          get_fracture_stiffness<dim>(param.Zn, param.Zt);

    // Scratch and copy data objects

    // Define update flags for faces
    const dii::UpdateFlags face_flags = dii::update_values | dii::update_gradients |
                                        dii::update_quadrature_points |
                                        dii::update_normal_vectors | dii::update_JxW_values;  

    // Strain functions for interface and boundarys
    auto interface_strain = [](const auto& fiv, const bool& ht, const auto& si, const auto& qp, const unsigned int comp) 
    { return get_strain_comp<dim>(fiv, ht, si, qp, comp); };
    
    auto boundary_strain = [](const auto& fv, const auto& si, const auto& qp) 
    { return get_strain_comp<dim>(fv, si, qp); };

    typename Mw_helper<dim>::ScratchData scratch_data(fe_collection, face_qcollection, face_flags,
                                                      param.vp, param.rho, param.penalty, param.mix_penalty, param.p_degree,
                                                      cg_id, dg1_id, dg2_id, param.fracture,
                                                      stiffness_tensor, fracture_stiffness,
                                                      interface_strain, boundary_strain);

    typename Mw_helper<dim>::CopyData copy_data;

    // MeshWorker workers
    // Empty cell worker
    using CellIteratorType = decltype(dof_handler.begin_active());
    std::function<void(const CellIteratorType &, 
                       typename Mw_helper<dim>::ScratchData &, 
                       typename Mw_helper<dim>::CopyData &)>  empty_cell_worker; 
    
    //Bpoundary worker (empty for Neumann boundary conditions)
    const auto boundary_worker =[this](const typename dii::DoFHandler<dim>::cell_iterator & cell,
                                 const unsigned int & face_n,
                                 typename Mw_helper<dim>::ScratchData & scratch_data,
                                 typename Mw_helper<dim>::CopyData & copy_data)
    {
      if(param.source == 4)
        {
         const unsigned int bid = cell->face(face_n)->boundary_id();
         if (bid != boundary_id_left && bid != boundary_id_right)
            return; // optional: Mw_helper<dim>::boundary_worker(...)

         for (const auto &pair : periodic_faces)
           {
            if (cell == pair.cell[0] && face_n == pair.face_idx[0])
              {
                constexpr unsigned int sf  = 0;  // set from pair if needed (refinement)
                constexpr unsigned int nsf = 0;
                Mw_helper<dim>::face_flux_worker(cell,
                                                 face_n,
                                                 sf,
                                                 pair.cell[1],
                                                 pair.face_idx[1],
                                                 nsf,
                                                 scratch_data,
                                                 copy_data);
                return;
              }
           }
        }
      else
          return;
    };

    //Face worker
    const auto face_worker = [](const typename dii::DoFHandler<dim>::cell_iterator &cell,
                                 const unsigned int &f,
                                 const unsigned int &sf,
                                 const typename dii::DoFHandler<dim>::cell_iterator &ncell,
                                 const unsigned int &nf,
                                 const unsigned int &nsf,
                                 typename Mw_helper<dim>::ScratchData &scratch_data,
                                 typename Mw_helper<dim>::CopyData &copy_data)
    {
        Mw_helper<dim>::face_flux_worker(cell, f, sf, ncell, nf, nsf, scratch_data, copy_data);
    };

  

    const auto copier = [&](const typename Mw_helper<dim>::CopyData &c)
    { 
        constraints.distribute_local_to_global(c.cell_matrix,
                                              c.local_dof_indices,
                                              stiffness_matrix);

        for (const typename Mw_helper<dim>::CopyDataFace &cdf : c.face_data)
          {
            constraints.distribute_local_to_global(cdf.cell_matrix,
                                                   cdf.joint_dof_indices,
                                                   stiffness_matrix);
          }
    };

    // MeshWorker loop over interior faces
    diiM::mesh_loop(dof_handler.begin_active(),
                    dof_handler.end(),
                    empty_cell_worker,
                    copier,
                    scratch_data,    
                    copy_data, 
                    diiM::assemble_boundary_faces |
                    diiM::assemble_own_interior_faces_once,
                    boundary_worker,
                    face_worker
                    );
    
                    
    std::cout << " flux assembly completed " << std::endl;
   
  }

  template <int dim>
  void waveMixSEM<dim>::assemble_boundary_matrix()
  {
    const double vp= param.vp;
    const double vs= param.vs;
    const double rho= param.rho;
    //check if source is isotropic

    const bool is_boundary_source =(param.source == 4);
    const bool point_source_abc =(! is_boundary_source);


    dii::hp::FEFaceValues<dim> hp_fe_facevalues (fe_collection, face_qcollection, 
                                                 dii::update_values | 
                                                 dii::update_quadrature_points |
                                                 dii::update_normal_vectors | 
                                                 dii::update_JxW_values);

    dii::FullMatrix<double> cell_boundary_matrix;
    std::vector<dii::types::global_dof_index> local_dof_indices;

    for (const auto &cell : dof_handler.active_cell_iterators())
         for (unsigned int f = 0; f < dii::GeometryInfo<dim>::faces_per_cell; ++f)
             if (cell->at_boundary(f))
                {

                  unsigned int bid = cell->face(f)->boundary_id();
                  bool is_periodic = (bid == boundary_id_left || bid == boundary_id_right);
                  bool boundary_source_abc= (is_boundary_source && !is_periodic);

                  if (point_source_abc || boundary_source_abc)
                  {
                  const auto &fe = cell->get_fe();
                  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
                  cell_boundary_matrix.reinit(dofs_per_cell,dofs_per_cell);
                  cell_boundary_matrix =0.0;
          
                  hp_fe_facevalues.reinit(cell, f);
                  const dii::FEFaceValues<dim> &fe_facevalues = hp_fe_facevalues.get_present_fe_values();
                  const std::vector<dii::Tensor<1, dim>> &normals = fe_facevalues.get_normal_vectors();

                  for (const unsigned int q_point : fe_facevalues.quadrature_point_indices())
                 {
                   const dii::Tensor<1, dim> &normal = normals[q_point];

                   for (const unsigned int i: fe_facevalues.dof_indices()) 
                    {
                      const unsigned int comp_i =
                                        fe.system_to_component_index(i).first;

                      const unsigned int mcomp_i = comp_i % dim;
                      double ni = normal[mcomp_i];
                      double phi_i = fe_facevalues.shape_value(i, q_point);

                      for (const unsigned int j: fe_facevalues.dof_indices())
                          {
                            const unsigned int comp_j =
                                        fe.system_to_component_index(j).first;                           
                            const unsigned int mcomp_j = comp_j % dim;                          
                            double nj = normal[mcomp_j];
                            double phi_j = fe_facevalues.shape_value(j, q_point);
                            // creating a Kronecker delta
                            double delta_ij = (mcomp_i == mcomp_j) ? 1.0 : 0.0;
                            double geometry_term = rho* (vp - vs)*ni *nj  + rho* vs*delta_ij;

                            cell_boundary_matrix(i, j) += phi_i * 
                                                          phi_j * 
                                                          geometry_term *
                                                          fe_facevalues.JxW(q_point) ;
                          }
                    }
                  }

                  local_dof_indices.resize(dofs_per_cell);
                  cell->get_dof_indices(local_dof_indices);
                  constraints.distribute_local_to_global(cell_boundary_matrix,
                                                         local_dof_indices,
                                                         boundary_matrix);


                  }
            }
  }



  template<int dim>
  void waveMixSEM<dim>::assemble_point_source()
  {
        //Mapping source to reference cell
   dii::types::fe_index fe_id ;

   if (!source_initialized) //only once
      {
        dii::Point<dim> p = source_function.source_pos;
          
        auto result =
        dii::GridTools::find_active_cell_around_point(mapping,
                                                   dof_handler,
                                                   p,
                                                   {},
                                                   1e-12);
         
        source_initialized = true;
        source_cell = result.first;
        source_pos_ref=result.second;
    
        std::cout << "p_ref: " << source_pos_ref << std::endl;
        out << "p_ref: " << source_pos_ref << std::endl;
          // Print all vertices in real (physical) coordinates
        for (unsigned int v = 0; v < dii::GeometryInfo<dim>::vertices_per_cell; ++v)
           {
             const dii::Point<dim> &vertex = source_cell->vertex(v);
             std::cout << "vertex " << v << ": " << vertex << std::endl;
             out << "vertex " << v << ": " << vertex << std::endl;
           }
       
           fe_id = source_cell->active_fe_index();
           if (fe_id != 0)
            throw std::runtime_error("Source cell is not a CG cell");
      }   
        // 1) Build a "quadrature" with exactly the reference point
    std::vector<dii::Point<dim>> quad_points(1);
    quad_points[0] = source_pos_ref;
    dii::Quadrature<dim> point_quadrature(quad_points);
          // 2) FEValues to evaluate shape function gradients at that point

    dii::FEValues<dim> fe_point_values(mapping,
                                      source_cell->get_fe(),
                                      point_quadrature,
                                      dii::update_gradients);

    fe_point_values.reinit(source_cell);
   
    const unsigned int dofs_per_cell = source_cell->get_fe().n_dofs_per_cell();
    std::vector<dii::types::global_dof_index> local_dof_indices(dofs_per_cell);
    source_cell->get_dof_indices(local_dof_indices);


    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      {
        const unsigned int comp_i =  
              source_cell->get_fe().system_to_component_index(i).first;

        double force_i =0.0;     
        if (param.source == 1 || param.source == 3)
        {
          //get the moment tensor
          dii::SymmetricTensor<2,dim> M = source_function.get_MT();
          // 2) ∂shape_i / ∂x_j for j=0..dim-1
          dii::Tensor<1,dim> grad_phi_i = 
                          fe_point_values.shape_grad_component(i, 0, comp_i);
 
          for (unsigned int j = 0; j < dim; ++j)
              force_i += M[comp_i][j] * grad_phi_i[j];  
        }
        else if (param.source == 2 )
        {
          const dii::Tensor<1,dim> force_vector =
                            source_function.get_force_vector();
          const double phi = 
              source_cell->get_fe().shape_value_component(i, source_pos_ref, comp_i);  
          force_i = force_vector[comp_i] * phi;
        }
        else
          throw std::runtime_error("Invalid source type");

        source_term(local_dof_indices[i]) += force_i;    

      }
  }
  template<int dim>
  void waveMixSEM<dim>::Neumann_plane_wave()
  {
    //labeling the extreme boundary faces with bid_min and bid_max to taper the source
    //boundary_cellsf: vector of pairs of cell iterators and face indices at boundary_id
    std::vector<std::pair<typename dii::DoFHandler<dim>::cell_iterator,
    unsigned int>> boundary_cellsf;
    const unsigned int boundary_id = 3, axis = 0;
    double Lmin = 0.0, Lmax = 0.0;
    const double fraction_nf = param.taper;
   
    mesh_op.label_extreme_boundary_faces(dof_handler, 
                                         boundary_id, bid_min, bid_max,
                                         axis, fraction_nf, Lmin, Lmax, boundary_cellsf);

    dii::hp::FEFaceValues<dim>  hp_fe_facevalues (fe_collection, face_qcollection,
                                                  dii::update_values | 
                                                  dii::update_quadrature_points |
                                                  dii::update_JxW_values);

    dii::Vector<double> cell_source_vector;
    std::vector<dii::types::global_dof_index> local_dof_indices;
    const dii::Tensor<1,dim> force_vector = source_function.get_force_vector();

    for (const auto &it : boundary_cellsf)
 
        {
          const auto &cell = it.first;
          const unsigned int f = it.second;

          const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
          cell_source_vector.reinit(dofs_per_cell);
          cell_source_vector =0.0;
  
          hp_fe_facevalues.reinit(cell, f);
          const dii::FEFaceValues<dim> &fe_facevalues = hp_fe_facevalues.get_present_fe_values();

          for (const unsigned int q_point : fe_facevalues.quadrature_point_indices())
            {
             for (const unsigned int i: fe_facevalues.dof_indices())
              {
                const unsigned int comp_i = cell->get_fe().system_to_component_index(i).first;
                const unsigned int mcomp_i = comp_i % dim;
                
                double phi_i = fe_facevalues.shape_value(i,q_point);
                cell_source_vector(i) +=  phi_i *force_vector[mcomp_i] 
                                        * fe_facevalues.JxW(q_point);
              }
            }
            local_dof_indices.resize(dofs_per_cell);
            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_source_vector,
                                                   local_dof_indices,
                                                   source_term);

        }
  }



  template <int dim>
   void waveMixSEM<dim>::timestep_loop()
  {

    const double rho=param.rho;
    const double dtf=param.dt_factor;
     //intial values - use hanging node constraints from setup_system
    
    dii::VectorTools::project (dof_handler, 
                               constraints,
                               qcollection,
                               InitialValuesU<dim>(),
                               old_solution_u);

    dii::VectorTools::project (dof_handler,
                               constraints,
                               qcollection,
                               InitialValuesV<dim>(),                       
                               old_solution_v);
    // constraining old solution vectors
    constraints.distribute(old_solution_u);
    constraints.distribute(old_solution_v);

    dii::Vector<double> tmp (solution_u.size());
    dii::Vector<double> forcing_terms (solution_u.size());

    source_term = 0.0;

    if (param.source > 4 )//point sourc
         throw std::runtime_error("Invalid source type");
    else if (param.source == 4)//Neumann boundary source
         Neumann_plane_wave ();
    else 
        assemble_point_source ();

    // Open output files for receivers' data (one per y-level)
    const std::string base = param.outputfile + "-";
    std::ofstream out_yup((base + "r_yup.dat").c_str());
    std::ofstream out_ydg((base + "r_ydg.dat").c_str());
    std::ofstream out_ydown((base + "r_ydown.dat").c_str());

    // Write headers
    write_receiver_header(out_yup, out_ydg, out_ydown);

    // time calculations
    // const double min_cellsize= extreme_cellsize(mesh);
    const double min_cellsize = dii::GridTools::minimal_cell_diameter(mesh);
    std::cout << "hmin = " << min_cellsize <<std::endl;

    double time_step =  dtf* calc_time_interval(min_cellsize, param.vp,param.p_degree);
    const double total_time=param.total_time;

    std::cout <<"time_step:  " <<time_step << std::endl;
    std::cout << "   Total_time: " << total_time << std::endl<< std::endl;

    out << "hmin = " << min_cellsize <<std::endl;
    out <<"time_step:  " <<time_step << std::endl;
    out << "   Total_time: " << total_time << std::endl<< std::endl;

    for (double time=0.0;
         time<=total_time;
         time+=time_step, ++timestep_number)
         {
   
          std::cout <<"timestep_number:  " <<timestep_number << std::endl;
          std::cout <<"current time:  " <<time << std::endl;

          out <<"timestep_number:  " <<timestep_number << std::endl;
          out <<"current time:  " <<time << std::endl;
          // ***Solving for v**
          solution_v = old_solution_v;

          stiffness_matrix.vmult (tmp, old_solution_u);
          tmp *= -time_step/rho;
          tmp.scale(inverse_mass_vector);
          solution_v += tmp;

        //Absorbing  Boundary condition for v
         boundary_matrix.vmult (tmp, old_solution_v);
         tmp *= -time_step/rho;
         tmp.scale(inverse_mass_vector);
         solution_v += tmp;

          //Include the source
          if (time <= source_function.source_duration())
          {          
            source_function.set_time(time); 
            double source_value= source_function.value();
            forcing_terms.equ(source_value * time_step/rho, source_term);//forcing_terms = source_value * time_step/rho * source_term
            forcing_terms.scale(inverse_mass_vector);

          }
          else
           forcing_terms =0.0;

          solution_v += forcing_terms;
          constraints.distribute(solution_v);

          //***Solving for u**
          solution_u = old_solution_u;
          solution_u.add(time_step, solution_v);
          constraints.distribute(solution_u);

          //calculating energ
          dealii::Vector<double> tmp(solution_v); 
          tmp.scale(mass_vector); 
          double total_energy=0.5*(solution_v * tmp) + 
                              0.5*stiffness_matrix.matrix_norm_square(solution_u);

          std::cout << "   Total energy: "<< total_energy << std::endl<<std::endl;
          out << "   Total energy: "<< total_energy << std::endl;
          
          if (total_energy>1e6) //stop the simulation if the total energy is too high (1e6 is a arbitrary value )
          {
           std::cout << "solution does not converge: stopping the simulation"
           << std::endl;
           out << "solution does not converge: stopping the simulation"
           << std::endl;
           break;
          }
       
          if (timestep_number % param.vtk_step == 0 && param.vtk)
            output_results();

          old_solution_u = solution_u;
          old_solution_v = solution_v;

          write_receiver_data(time, out_yup, out_ydg, out_ydown);
          //if (timestep_number >= 10)
             //break;
         }
  }


  template <int dim>
    void waveMixSEM<dim>::output_results () const
    {
      dii::DataOut<dim> data_out;
      data_out.set_flags (dii::DataOutBase::VtkFlags(param.p_degree)); // Since fe.degree is 2
      data_out.attach_dof_handler (dof_handler);
      
      // Add individual components directly using component selection
      data_out.add_data_vector(solution_u, std::vector<std::string>{"u0_cg", "u1_cg", "u2_dg", "u3_dg"});

      // Compute norms as sqrt(u0^2 + u1^2) and sqrt(u2^2 + u3^2) per support point
      const unsigned int n_dofs = solution_u.size();
      dii::Vector<double> norm_cg(n_dofs), norm_dg(n_dofs);
      norm_cg = 0.0;
      norm_dg = 0.0;

      for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
        const unsigned int n_support_points = dofs_per_cell / 2;  // 2 components per point (CG or DG)
        std::vector<dii::types::global_dof_index> local_dof_indices(dofs_per_cell);
        cell->get_dof_indices(local_dof_indices);
        const unsigned int fe_index = cell->active_fe_index();

        // Per support point: store the two component values (vector indexed by support_point)
        std::vector<std::array<double, 2>> sp_vals(n_support_points, {{0.0, 0.0}});

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int global_dof = local_dof_indices[i];
          const unsigned int comp = cell->get_fe().system_to_component_index(i).first;
          const unsigned int support_point = cell->get_fe().system_to_component_index(i).second;
          const double value = solution_u[global_dof];

          if (fe_index == 0 && (comp == 0 || comp == 1))  // CG cell
            sp_vals[support_point][comp] = value;
          else if (fe_index == 1 && (comp == 2 || comp == 3))  // DG cell
            sp_vals[support_point][comp - 2] = value;
        }

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int global_dof = local_dof_indices[i];
          const unsigned int comp = cell->get_fe().system_to_component_index(i).first;
          const unsigned int support_point = cell->get_fe().system_to_component_index(i).second;

          if (fe_index == 0 && (comp == 0 || comp == 1))  // CG
          {
            const double u0 = sp_vals[support_point][0];
            const double u1 = sp_vals[support_point][1];
            norm_cg[global_dof] = std::sqrt(u0 * u0 + u1 * u1);
          }
          else if (fe_index == 1 && (comp == 2 || comp == 3))  // DG
          {
            const double u2 = sp_vals[support_point][0];
            const double u3 = sp_vals[support_point][1];
            norm_dg[global_dof] = std::sqrt(u2 * u2 + u3 * u3);
          }
        }
      }

      // Add norms
      data_out.add_data_vector(norm_cg, "norm_cg");
      data_out.add_data_vector(norm_dg, "norm_dg");


      const unsigned int n_subdivisions = param.p_degree; 

      data_out.build_patches(n_subdivisions);
      
      const std::string output_dir = "./vtk/";
      const std::string filename = output_dir +
                                  param.outputfile + "-" +
                                  dii::Utilities::int_to_string (timestep_number, 3) +
                                  ".vtk";
      std::ofstream output (filename.c_str());
      data_out.write_vtk (output);
    }
  

template<int dim>
  void waveMixSEM<dim>::write_receiver_header(std::ofstream &out_yup,
                                            std::ofstream &out_ydg,
                                            std::ofstream &out_ydown)
  {
    const unsigned int nr = param.receivers_x.size();
    if (receivers_pos.empty() || nr == 0 || receivers_pos.size() != 3 * nr)
      return;

    auto write_header = [&](std::ofstream &f, unsigned int offset, const char *label)
    {
      f << "# " << label << "  y = " << receivers_pos[offset][1] << std::endl;
      f << "# x_coords:";
      for (unsigned int i = 0; i < nr; ++i)
        f << " " << receivers_pos[offset + i][0];
      f << std::endl;
      f << "# time";
      for (unsigned int i = 0; i < nr; ++i)
        f << "  ux_r" << i << "  uy_r" << i;
      f << std::endl;
    };

    write_header(out_yup,      0,      "CG yup");
    write_header(out_ydg,      nr,     "DG ydg");
    write_header(out_ydown,  2 * nr,   "CG ydown");
  }

template<int dim>
  void waveMixSEM<dim>::write_receiver_data(const double time,
                                         std::ofstream &out_yup,
                                         std::ofstream &out_ydg,
                                         std::ofstream &out_ydown)
  {
    const unsigned int nr = param.receivers_x.size();
    if (receivers_pos.empty() || nr == 0)
      return;

    // receivers_pos layout:
    //   [0    .. nr-1]    : CG at yup    (components 0,1)
    //   [nr   .. 2*nr-1]  : DG at ydg    (components 2,3)
    //   [2*nr .. 3*nr-1]  : CG at ydown  (components 0,1)

    dii::Vector<double> value(4);

    // yup (CG): time ux_r0 uy_r0 ux_r1 uy_r1 ...
    out_yup << time;
    for (unsigned int i = 0; i < nr; ++i)
    {
      dii::VectorTools::point_value(dof_handler, solution_u, receivers_pos[i], value);
      out_yup << " " << value(0) << " " << value(1);
    }
    out_yup << std::endl;

    // ydg (DG): time ux_r0 uy_r0 ux_r1 uy_r1 ...
    out_ydg << time;
    for (unsigned int i = nr; i < 2 * nr; ++i)
    {
      dii::VectorTools::point_value(dof_handler, solution_u, receivers_pos[i], value);
      out_ydg << " " << value(2) << " " << value(3);
    }
    out_ydg << std::endl;

    // ydown (CG): time ux_r0 uy_r0 ux_r1 uy_r1 ...
    out_ydown << time;
    for (unsigned int i = 2 * nr; i < 3 * nr; ++i)
    {
      dii::VectorTools::point_value(dof_handler, solution_u, receivers_pos[i], value);
      out_ydown << " " << value(0) << " " << value(1);
    }
    out_ydown << std::endl;
  }

template<int dim>
  void waveMixSEM<dim>::run()
  {
    make_mesh();
    setup_system();
    assemble_MKmatrices();
    include_interface_fluxes();
    assemble_boundary_matrix(); 
    timestep_loop();

  }

