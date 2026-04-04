template <int dim>
  Source<dim>::Source(const Parameters  & param_ ):
  dii::Function<dim>(dim),
  param(param_),
  fo(param_.fo),
  source_rad(param_.source_rad),
  stf_type(param_.stf_type),
  M11(param_.M11), M12(param_.M12), M22(param_.M22), M21(param_.M21),
  loc_x(param_.loc_x), loc_y(param_.loc_y), loc_z(param_.loc_z),
  source_pos(get_source_position(param))
  { 

  }

  template <int dim>
  dii::Point<dim> Source<dim>::get_source_position( const Parameters & prm) 
  {
    dii::Point<dim> p;

    if constexpr (dim == 2)
      {
        p[0] = prm.loc_x;
        p[1] = prm.loc_y;
      }
    else if constexpr (dim == 3)
      {
        p[0] = prm.loc_x;
        p[1] = prm.loc_y;
        p[2] = prm.loc_z;
      }
    else
      {
        static_assert(dim == 2 || dim == 3, "Unsupported dimension");
      }

    return p;

  }

template <int dim>
  dii::SymmetricTensor<2, dim> Source<dim>::get_MT() const
  {
     dii::SymmetricTensor<2, dim> tmp;
     tmp.clear();  
         // Diagonal terms
          // Diagonals
    tmp[0][0] = M11;
    tmp[1][1] = M22;

    // Off-diagonals
    if (param.source == 1)
      {
        tmp[0][1] = M12;
        tmp[1][0] = M12;  // MUST BE EXPLICIT
      }
    else if (param.source == 3)
      {
        tmp[0][1] = M12;
        tmp[1][0] = -M12;  // MUST BE EXPLICIT
      }
    else
      {
        throw std::runtime_error("Unsupported source type");
      }
    return tmp;

  }

  template <int dim>
  dii::Tensor<1, dim> Source<dim>::get_force_vector() const
  {
    dii::Tensor<1, dim> tmp;
    const double pi = dii::numbers::PI;
    const double phi_rad = param.phi * pi / 180.0;
    tmp[0] = std::cos(phi_rad);
    tmp[1] = std::sin(phi_rad);
    return tmp;
  }

  template <int dim>
  bool Source<dim>::is_isotropic() const
  {
    bool is_isotropic = (param.source == 1 && (param.M11 == param.M22 && param.M12 == 0.0));
    bool is_Patboundary=(param.source == 4 && param.phi == -90.0);

    return (is_isotropic || is_Patboundary);
  }

template <int dim>
  bool Source<dim>::is_isotropic(const Parameters & prm)
  {
    bool is_isotropic = (prm.source == 1 && (prm.M11 == prm.M22 && prm.M12 == 0.0));
    bool is_Patboundary=(prm.source == 4 && prm.phi == -90.0);

    return (is_isotropic || is_Patboundary);
  }

template <int dim>
    double Source<dim>::source_duration() const
    {
        return 4.0/fo;
    }

template <int dim>
    double Source<dim>::source_time(const double & t )const
    {
       const double to = 2.0/fo;
       const double pi = dii::numbers::PI;
       const double f2 = pi*pi*fo*fo;
       const double time2 = (t-to)*(t-to);
       
       return (1.0-2.0*f2*time2)*std::exp(-f2*time2);
    }

template <int dim>
    double Source<dim>::source_space(const dii::Point<dim> & p)const
    { 
      const double pi = dii::numbers::PI;
      const double c  = 1.0/ (std::pow( std::sqrt(pi)*source_rad, dim));
      const double r2=source_rad*source_rad;
      const double d=source_pos.distance(p);

      return c*std::exp(- d*d/r2); 
    }
    
template <int dim>
    double Source<dim>::value(const dii::Point<dim>  & p,
                              const unsigned int component) const
    { (void) p;
      (void) component;
      Assert (component == 0, dii::ExcIndexRange(component, 0, 1));

      const double t = this->get_time();

      if (t > source_duration() )
        return 0.0;

      return source_time(t);
    }

template <int dim>
    std::vector<typename dii::DoFHandler<dim>::active_cell_iterator>
    Source<dim>::find_source_cells(const dii::DoFHandler<dim> & dof_handler) const
    {   
        const double factor = 2.0;
        std::vector<typename dii::DoFHandler<dim>::active_cell_iterator> source_cells;

        for (const auto &cell : dof_handler.active_cell_iterators())
        {
         // cell diameter
         const double h = cell->diameter();

            // compute safe cutoff distance
         const double cutoff = factor * source_rad + 0.5 * h;

         // include cell if center is within cutoff
         if (cell->center().distance(source_pos) <= cutoff)
          {
            source_cells.push_back(cell);
          }
        }

        return source_cells;
    }

template <int dim>
   Mesh_operations<dim>::Mesh_operations(const Parameters & param_):
    param(param_),
    vs (param_.vs), 
    vp(param_.vp),
    wavesamples(param_.nsamples), 
    sourcesampling(param_.sourcesampling), 
    fo(param_.fo), 
    source_rad(param_.source_rad),
    boundary_points(get_boundary_points()),
    hbg(0.0),
    outer_radius(0.0),
    Lx(param_.Lx),
    Ly(param_.Ly)
    { 
    }
  

  template <int dim>
   std::vector<dii::Point<dim>> Mesh_operations<dim>::get_boundary_points() const 
   { 
    std::vector<dii::Point<dim>> bpoints(2);
    const double Lx = param.Lx;
    const double Ly = param.Ly; 
    
     if constexpr ( dim==2)
      { 
        bpoints[0]=dii::Point<dim>(0.0, 0.0);
        bpoints[1]=dii::Point<dim>(Lx, -Ly);
      }
     
    else if constexpr ( dim==3)
      { 
        bpoints[0]=dii::Point<dim>(0.0, 0.0, 0.0);
        bpoints[1]=dii::Point<dim>(Lx, Ly, -Ly);
      
      }

    else
      {
         static_assert(dim == 2 || dim == 3, "Unsupported dimension.");
      }
    return bpoints;
   }



  template <int dim>
   double Mesh_operations<dim>::get_velocity() const
   {
    bool isotropic= Source<dim>::is_isotropic(param);
    if (isotropic)
      return vp;
    else return vs;
   }


template <int dim>
   void Mesh_operations<dim>::mesh_create(dii::Triangulation<dim> & mesh) const 
   {
    double ratio, Lmax, Lmin, eps=1e-10;
    enum direction {x,y};
    direction axis;

     if (Lx>Ly+ eps || std::abs(Lx-Ly) < eps)
       {
         Lmax = Lx;
         Lmin = Ly;
         axis = x;

       }
    else
       {
         Lmax = Ly;
         Lmin = Lx;
         axis = y;
       }

    ratio = Lmax/Lmin;
    
    // create a new mesh and refine it globally
    const double c = get_velocity();
    const double lambda_s = c/fo;
    const double h = lambda_s/wavesamples;

    if (Lmin/h > 1)
    {
      const unsigned int refine_min = static_cast<unsigned int> ( std::ceil(std::log2(Lmin/h)) );
      const unsigned int rep_min = 1u << refine_min;  //2^refine_min
      const unsigned int rep_max =
           std::max(1u, static_cast<unsigned int>(std::ceil(rep_min * (ratio))));

       std::vector<unsigned int> repetitions_dim(dim);
       repetitions_dim[0] = (axis == x) ? rep_max : rep_min;
       repetitions_dim[1] = (axis == x) ? rep_min : rep_max;


       mesh.clear();
       dii::GridGenerator::subdivided_hyper_rectangle(mesh,
                                                     repetitions_dim,
                                                     boundary_points[0], 
                                                     boundary_points[1],
                                                     true) ;
     
   }
  }

  template <int dim>
   void Mesh_operations<dim>::label_extreme_boundary_faces(dii::DoFHandler<dim> &dof_handler,
                                                           const unsigned int boundary_id,
                                                           const unsigned int bid_min,
                                                           const unsigned int bid_max,
                                                           const unsigned int axis, // 0 for x, 1 for y
                                                           const double fraction_nf,
                                                           double &Lmin, double &Lmax,
                                                           std::vector<std::pair<typename dii::DoFHandler<dim>::cell_iterator,
                                                           unsigned int>> &faces) const
   {


    if (axis >= dim)
       throw std::runtime_error("axis must be in [0, " + std::to_string(dim - 1) +
                   "] for dim=" + std::to_string(dim) + ".");

    if (fraction_nf >= 0.5)
        throw std::runtime_error(
           "fraction_nf must be < 0.5 (otherwise min/max boundary bands can overlap).");

    faces.clear();
 
  // --- Step 1: collect faces with given boundary_id ---
   for (auto &cell : dof_handler.active_cell_iterators())
     {
     for (unsigned int f = 0; f < dii::GeometryInfo<dim>::faces_per_cell; ++f)
       {
         if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == boundary_id)
         {
             faces.emplace_back(cell, f);
         }
       }
    }
 
   const unsigned int nfaces = faces.size();
   if (nfaces == 0)
        throw std::runtime_error("No faces found with boundary_id: " + 
                                 std::to_string(boundary_id) + ".");
 

       // --- Step 2: sort faces by m-th coordinate ---
   std::sort(faces.begin(), faces.end(),[axis](const auto &a, const auto &b)
          {
             const double xa = a.first->face(a.second)->center()[axis];
             const double xb = b.first->face(b.second)->center()[axis];
              return xa < xb;
           });
   
   // --- Step 3 : assign faces with bid_min→ 31 ---
   const double Lref=param.Lx*fraction_nf;
   Lmin = 0.0;

   //for (unsigned int i = 0; i < nextreme_out; ++i)
   for ( auto & it: faces)
       {
         auto &cell = it.first;
         const unsigned int f = it.second;
         auto face = cell->face(f);
   
         face->set_boundary_id(bid_min);
   
         // accumulate length along m-direction
         Lmin += cell->extent_in_direction(axis);
         if (Lmin >= Lref)
           break;
       }
   
       // --- Step 4: assign faces with bid_max→ 32 ---

   Lmax = 0.0;

   //for (unsigned int i = 0; i < nextreme_out; ++i)
   for (auto rit =faces.rbegin(); rit != faces.rend(); ++rit)  
       {
        auto &cell =rit->first;
        const unsigned int f = rit->second;
        auto face = cell->face(f);
   
        face->set_boundary_id(bid_max);
   
        Lmax += cell->extent_in_direction(axis);
        if (Lmax >= Lref)
          break;
       }
   }

template <int dim>
  void Mesh_operations<dim>::source_refine(dii::Triangulation<dim> & mesh) 
  {
    const dii::Point<dim> source_pos = Source<dim>::get_source_position(param);
    const double target_h = 2.0 * source_rad / sourcesampling;

    // Background cell size from the already-refined uniform mesh
    const double h_bg = mesh.begin_active()->diameter();

    if (h_bg <= target_h * 1.05)
    {
        hbg = h_bg;
        std::cout << "Source refinement: h_bg (" << h_bg
                  << ") <= target_h (" << target_h << "), no refinement needed." << std::endl;
        return;
    }

    // Doubling layers needed: target_h -> 2*target_h -> ... -> h_bg
    const int nlayers_auto = static_cast<int>(std::ceil(std::log2(h_bg / target_h)));

    const double safety_factor = 1.2; 
    const double inner_radius = safety_factor * source_rad;

    // Transition zone: layer k has cell size target_h*2^k, spans ~1.5x its size
    const double layer_thickness_factor = 1.5;
    const double transition = layer_thickness_factor * target_h
                              * (std::pow(2.0, nlayers_auto) - 1.0);
    outer_radius = inner_radius + transition;

    std::cout << "Source refinement: target_h=" << target_h 
              << ", h_bg=" << h_bg
              << ", nlayers=" << nlayers_auto
              << ", inner_r=" << inner_radius 
              << ", outer_r=" << outer_radius << std::endl;

    auto get_target_h_at_distance = 
    [target_h, h_bg, inner_radius, this]
    (double d) -> double 
    {
        if (d <= inner_radius)
            return target_h;

        if (d >= outer_radius)
            return h_bg;

        const double t = (d - inner_radius) / (outer_radius - inner_radius);
        return target_h * std::exp(t * std::log(h_bg / target_h));
    };

    while (true)
    {
        bool flag_refine = false;
        for (auto &cell : mesh.active_cell_iterators())
        {
            double h = cell->diameter();
            double d = cell->center().distance(source_pos);
            double target_h_local = get_target_h_at_distance(d);

            if (h > target_h_local * 1.05)
            {
                cell->set_refine_flag();
                flag_refine = true;
            }
        }

        if (flag_refine)
            mesh.execute_coarsening_and_refinement();
        else
            break;
    }
    
    for (const auto &cell : mesh.active_cell_iterators())
    {
        double d = cell->center().distance(source_pos);
        if (d >= 1.1 * outer_radius)
        {
            hbg = cell->extent_in_direction(1);
            break;
        }
    }
    
    if (hbg == 0.0)
         throw std::runtime_error(
         "No uniform mesh region found outside refinement zone. "
         "Ensure mesh domain is large enough.");

  }

template <int dim>
  void Mesh_operations<dim>::print_mesh(dii::Triangulation<dim> &mesh) const
  {
  const std::string filename = "./vtk/" + param.outputfile + "_mesh.vtk";
  std::ofstream out(filename);
  dii::GridOut grid_out;
  grid_out.write_vtk(mesh, out);
  }

  template <int dim>
  void Mesh_operations<dim>::refine_dg_rows(dii::Triangulation<dim> & mesh, 
                                            const double & top_dg_y, 
                                            const double & bottom_dg_y)
  {
    if (!param.refine_dg_rows)
      return;

    const double tol = 1e-10 * std::max(hbg, 1.0);
    const double y_inner_lo = bottom_dg_y - 1.0 * hbg - tol;
    const double y_inner_hi = top_dg_y + 1.0 * hbg + tol;
    const double y_outer_hi = top_dg_y + 2.0 * hbg + tol;
    const double y_outer_lo = bottom_dg_y - 2.0 * hbg - tol;

    auto in_inner = [&](double cy) { return cy >= y_inner_lo && cy <= y_inner_hi; };
    auto in_outer_only = [&](double cy) {
      return (cy > y_inner_hi - 2.0 * tol && cy <= y_outer_hi) ||
             (cy >= y_outer_lo && cy < y_inner_lo + 2.0 * tol);
    };

    for (auto &cell : mesh.active_cell_iterators())
      {
        const double cy = cell->center()[1];
        const bool dg = (cell->material_id() != 0);
        if (dg || in_inner(cy) || in_outer_only(cy))
          cell->set_refine_flag();
      }
    mesh.execute_coarsening_and_refinement();

    for (auto &cell : mesh.active_cell_iterators())
      {
        const double cy = cell->center()[1];
        if (cell->material_id() != 0 || in_inner(cy))
          cell->set_refine_flag();
      }
    mesh.execute_coarsening_and_refinement();

    std::cout << " refine_dg_rows: 2 levels inner (DG±1*hbg), 1 level outer (±2*hbg). Active cells: "
              << mesh.n_active_cells() << std::endl;
  }


  template <int dim>
   double Mesh_operations<dim>::get_fcoord() const
  { 
    const bool is_boundary_source =(param.source == 4 );
    double pos_y=0.0;
     
    if(!is_boundary_source)
      pos_y = Source<dim>::get_source_position(param)[1];

    const double c = get_velocity();

    const double min_dist=2.9*c/fo;

    const double fcoord = pos_y - min_dist;

    if (std::abs(fcoord) > Ly)
      throw std::runtime_error("Fracture coordinate is outside mesh extent: |fcoord| = "
                               + std::to_string(std::abs(fcoord)) + " > Ly = "
                               + std::to_string(Ly) + ".");

    std::cout << " Suggested fracture coordinate: " << fcoord << std::endl;

    std::cout << " min distance from source: " << min_dist << std::endl;

    return fcoord;
      
  }

template <int dim>
  std::vector<dii::Point<dim>> Mesh_operations<dim>::receiver_coordinates(const double & h,
                                                                          const double & yfrac, 
                                                                          const double & bottom_dg_y,
                                                                          const unsigned int & ncells,
                                                                          const std::vector<double> & receivers_x
                                                                          ) const
  { 
    std::vector<dii::Point<dim>> receivers_pos;
    if (receivers_x.empty())
      return receivers_pos;


    const bool is_boundary_source =(param.source == 4 );
    double source_y=0.0;

    if(!is_boundary_source)
      source_y = Source<dim>::get_source_position(param)[1];


    // Round x coordinates to 1 decimal
    std::vector<double> xcoords(receivers_x.size());
    for (unsigned int i = 0; i < receivers_x.size(); ++i)
      xcoords[i] = std::round(receivers_x[i] * 10.0) / 10.0;

    // y-coordinates for each receiver line, rounded to 1 decimal
    const double c = get_velocity();
    const double lambda=c/fo;
    double yup   = source_y- lambda;
    double ydown = bottom_dg_y - (ncells + 0.5) * h;
    double h_dg  = param.refine_dg_rows ? 0.5 * h : h;
    double ydg   = yfrac + 0.5 * h_dg;
    yup   = std::round(yup   * 10.0) / 10.0;
    ydown = std::round(ydown * 10.0) / 10.0;
    ydg   = std::round(ydg   * 10.0) / 10.0;

    auto point_at = [](double x, double y) -> dii::Point<dim>
    {
      if constexpr (dim == 2)
        return dii::Point<dim>(x, y);
      else
        return dii::Point<dim>(x, y, 0.0);
    };

    // CG receivers at yup
    for (const double x : xcoords)
      receivers_pos.push_back(point_at(x, yup));

    // DG receivers at ydg
    for (const double x : xcoords)
      receivers_pos.push_back(point_at(x, ydg));

    // CG receivers at ydown
    for (const double x : xcoords)
        receivers_pos.push_back(point_at(x, ydown));

    return receivers_pos;
  }

template <int dim>
  void Mw_helper<dim>::face_flux_worker (const typename dii::DoFHandler<dim>::cell_iterator & cell,
                                         const unsigned int & f,
                                         const unsigned int & sf,
                                         const typename dii::DoFHandler<dim>::cell_iterator & ncell,
                                         const unsigned int & nf,
                                         const unsigned int & nsf,
                                         Mw_helper<dim>::ScratchData & scratch_data,
                                         Mw_helper<dim>::CopyData & copy_data)
  { 
     // Only apply flux if:
    bool cell_is_dg1=(cell->material_id() == scratch_data.dg1_id);
    bool cell_is_dg2=(cell->material_id() == scratch_data.dg2_id);
    bool ncell_is_dg1=(ncell->material_id() == scratch_data.dg1_id);
    bool ncell_is_dg2=(ncell->material_id() == scratch_data.dg2_id);
// - One cell is CG (cg_id) and neighbor is DG (dg1_id or dg2_id), OR
// - Both cells are DG (any combination of dg1_id and dg2_id)

  bool cell_is_cg = (cell->material_id() == scratch_data.cg_id);
  bool ncell_is_cg = (ncell->material_id() == scratch_data.cg_id);
  bool cell_is_dg = (cell->material_id() == scratch_data.dg1_id || cell->material_id() == scratch_data.dg2_id);
  bool ncell_is_dg = (ncell->material_id() == scratch_data.dg1_id || ncell->material_id() == scratch_data.dg2_id);

// Skip if both are CG (internal CG faces don't need DG flux terms)
  if (cell_is_cg && ncell_is_cg)
       return;
  
    // DG1-DG2 face (fracture face) - only use fracture flux if enabled
  const bool fracture_face = (scratch_data.fracture && ((cell_is_dg1 && ncell_is_dg2) || (cell_is_dg2 && ncell_is_dg1)));
  const bool dg_dg_face = (cell_is_dg && ncell_is_dg);
  const bool cg_dg_face =  (cell_is_cg && ncell_is_dg) || (cell_is_dg && ncell_is_cg) ;

  if (cg_dg_face || dg_dg_face)
    {
      
      //penalty parameter eta
      /*******************/
      const double vp = scratch_data.vp;
      const double rho = scratch_data.rho;
      const double R = cg_dg_face ? scratch_data.R_mix : scratch_data.R;
      const int p_degree = scratch_data.p_degree;
      const double face_size =  cell->face(f)->measure();
      const double face_size_minus = ncell->face(nf)->measure();
      double face_size_avg = 0.5 * (face_size + face_size_minus);

      int S;
      if (cg_dg_face)
         S = -1;   // symmetry term at CG/DG interface
      else
         S = -1; //symmetry term flag
      
      const double eta =  R*vp*vp*rho* (p_degree+ 1)*(p_degree+1) /face_size_avg;
     const dii::SymmetricTensor<2,dim> fracture_stiffness = scratch_data.fracture_stiffness;
  
      //Flux asembly
      /*******************/

      //get FEInterfaceValues object and reinit for this face pair

    dii::FEInterfaceValues<dim> &fe_iv = scratch_data.fe_interface_values;

     const unsigned int fe_idx_cell  = cell->active_fe_index();
     const unsigned int fe_idx_ncell = ncell->active_fe_index();
     const unsigned int q_idx = fe_idx_cell; 
     const unsigned int mapping_idx = 0;
     fe_iv.reinit(cell, f, sf, ncell, nf, nsf, q_idx, mapping_idx, fe_idx_cell, fe_idx_ncell);

     // fe_iv.reinit(cell, f, sf, ncell, nf, nsf);

      // Add new element and get reference
      Mw_helper<dim>::CopyDataFace & copy_data_face = copy_data.face_data.emplace_back();  

      const unsigned int n_dofs_face = fe_iv.n_current_interface_dofs();
      copy_data_face.joint_dof_indices  = fe_iv.get_interface_dof_indices();
      copy_data_face.cell_matrix.reinit(n_dofs_face, n_dofs_face);
      
      const std::vector<double>  &JxW    = fe_iv.get_JxW_values();
      const std::vector<dii::Tensor<1, dim>> &normals = fe_iv.get_normal_vectors();
      const auto &q_points = fe_iv.get_quadrature_points();
      const unsigned int nq_points = q_points.size();

      // static int call_count = 0;
      // ++call_count;
      // std::ostringstream oss;
      // if (dg_dg_face)
      // {
      //     oss << "\n========== Call #" << call_count
      //         << " | Cell, face " << cell->active_cell_index() << ", " << f
      //         << " -> NCell, nface " << ncell->active_cell_index() << ", " << nf
      //         << " | FE id: " << fe_idx_cell << ", " << fe_idx_ncell
      //         << " | nq_points: " << nq_points << "\n";
      // }

      struct ShapeComponentInfo
      {
        unsigned int component;
        bool has_here;
        bool has_there;
      };

      //lambda function to get comp of interface dof
      auto get_shape_component =
        [&](unsigned int interface_dof) -> ShapeComponentInfo
        {
         // Test if this dof is owned by 'here = 0 = true' (cell) or 'there=1=false' (ncell)
         const auto dof_pair = fe_iv.interface_dof_to_dof_indices(interface_dof);
         const bool has_here = (dof_pair[0] != dii::numbers::invalid_unsigned_int);
         const bool has_there = (dof_pair[1] != dii::numbers::invalid_unsigned_int);
         const bool owned_here = has_here;
         const unsigned int cell_dof = owned_here ? dof_pair[0] : dof_pair[1];
        
         // Get component of the vector fe  from the owning side
         const unsigned int component = owned_here ?
           cell->get_fe().system_to_component_index(cell_dof).first :
           ncell->get_fe().system_to_component_index(cell_dof).first;

         return {component, has_here, has_there};
      };

      // Loops
      for (unsigned int point = 0; point < nq_points; ++point)
        {

          const dii::Tensor<1, dim> &n = normals[point];

          // if (cg_dg_face)
          //     oss << "\n========== Quadrature Point " << point << " nq_points=" << nq_points << " ==========\n";

          // Precompute tractions for all DOFs at this quadrature point
          std::vector<dii::Tensor<1, dim>> traction_avg(n_dofs_face);
    
          for (unsigned int j = 0; j < n_dofs_face; ++j)
           {
            
            // Get component of the vector fe  from the owning side
            const auto shape_info = get_shape_component(j);
            const unsigned int comp = shape_info.component;

            // Get strain (only on the side where this interface dof exists)
            dii::SymmetricTensor<2,dim> eps_here;
            dii::SymmetricTensor<2,dim> eps_there;
            if (shape_info.has_here)
              eps_here = scratch_data.fstrain_function(fe_iv, true, j, point, comp);
            if (shape_info.has_there)
              eps_there = scratch_data.fstrain_function(fe_iv, false, j, point, comp);

            dii::Tensor<1, dim> t_here = scratch_data.stiffness_tensor * eps_here * n;  
            dii::Tensor<1, dim> t_there = scratch_data.stiffness_tensor * eps_there * n;

             // Average traction: 

            double weight;
        
            if (cg_dg_face)
               weight = (comp < 2) ? 0.0 : 1.0;   // CG comps 0,1 -> 0; DG 2,3 -> 1
            else //dg_dg_face
               weight = 0.5;

            // Average traction: weight* (t_here + t_there)
            traction_avg[j] = weight* (t_here + t_there);

            // if (cg_dg_face)
            //     oss << "  DOF j=" << j << " comp=" << comp
            //         << " | has_here=" << shape_info.has_here
            //         << " | has_there=" << shape_info.has_there
            //         << " | t_here=[" << t_here[0] << "," << t_here[1] << "]"
            //         << " | t_there=[" << t_there[0] << "," << t_there[1] << "]\n";
           }

            // Nested loops: compute cell matrix
            for (unsigned int i = 0; i < n_dofs_face; ++i)
            { 
              const unsigned int comp_i = get_shape_component(i).component;
              const unsigned int mcomp_i = comp_i % dim;
              const double jump_i = fe_iv.jump_in_shape_values(i, point, comp_i);
   

              for (unsigned int j = 0; j < n_dofs_face; ++j)
                  {                    
                    const unsigned int comp_j = get_shape_component(j).component;
                    const unsigned int mcomp_j = comp_j % dim;

                    const double jump_j = fe_iv.jump_in_shape_values(j, point, comp_j);
                    
                    // Penalty term
                    double penalty_term = 0.0;
                    if (fracture_face)
                      penalty_term = fracture_stiffness[mcomp_i][mcomp_j] * jump_i * jump_j;  // Z^{-1}_{ij} [u_j][v_i]
                    else if (mcomp_i == mcomp_j)
                      penalty_term = eta * jump_i * jump_j;
                        
                    double consistency_term = 0.0;
                    double symmetry_term = 0.0;
                       
                  if(!fracture_face)
                    {
                    consistency_term = -traction_avg[j][mcomp_i] * jump_i;
                     symmetry_term = S*traction_avg[i][mcomp_j] * jump_j;
                    }

                    // if (cg_dg_face && (std::abs(penalty_term) > 1e-10 || std::abs(consistency_term) > 1e-10))
                    // {
                    //     oss << "    (i=" << i << ",j=" << j << ") comp_i=" << comp_i << " comp_j=" << comp_j
                    //         << " | jump_i=" << jump_i << " jump_j=" << jump_j;
                    //     if (std::abs(penalty_term) > 1e-10)
                    //         oss << " | penalty=" << penalty_term;
                    //     if (std::abs(consistency_term) > 1e-10)
                    //         oss << " | consistency=" << consistency_term;
                    //     if (std::abs(symmetry_term) > 1e-10)
                    //         oss << " | symmetry=" << symmetry_term;
                    //     oss << "\n";
                    // }

                    copy_data_face.cell_matrix(i, j) += 
                    (
                     // Term 1: {t_j} · [v_i]
                     consistency_term                 
                    // Term 2: penalty term
                     + penalty_term            
                    // Term 3: symmetry term (if needed)
                     + symmetry_term
                  
                    ) * JxW[point];
                  }
            }

        } // end loop over quadrature points

      // if (cg_dg_face)
      //     std::cout << oss.str();

      } // end if cell is CG/DG or DG/DG
  }

template <int dim>
  void Mw_helper<dim>::boundary_worker (const typename dii::DoFHandler<dim>::cell_iterator & cell,
                                        const unsigned int & face_n,
                                        Mw_helper<dim>::ScratchData & scratch_data,
                                        Mw_helper<dim>::CopyData & copy_data)
  {
    
     bool cell_is_dg = (cell->material_id() == scratch_data.dg1_id || 
                       cell->material_id() == scratch_data.dg2_id);

    if (!cell_is_dg)
         return;  // CG boundary cell; skip DG flux assembly

// Assemble DG boundary flux terms...

    //penalty parameter
    /*******************/
    const double vp = scratch_data.vp;
    const double rho = scratch_data.rho;
    const double R = scratch_data.R;
    const int p_degree = scratch_data.p_degree;
    const double face_size =  cell->face(face_n)->measure();

    const double eta=R*vp*vp*rho* (p_degree+ 1)*(p_degree+1)/face_size;
    const int S = -1; //symmetry term flag
   
    //Flux asembly
    /*******************/

    //get FEFaceValues object and reinit for this face
    scratch_data.fe_interface_values.reinit(cell, face_n);
    const dii::FEFaceValuesBase<dim> &fe_face = scratch_data.fe_interface_values.get_fe_face_values(0);
  
    const unsigned int dofs_per_cell = fe_face.get_fe().n_dofs_per_cell();
    copy_data.reinit(cell, dofs_per_cell);
    const std::vector<double>         &JxW     = fe_face.get_JxW_values();
    const std::vector<dii::Tensor<1, dim>> &normals = fe_face.get_normal_vectors();

    const auto &q_points = fe_face.get_quadrature_points();

    for (unsigned int point = 0; point < q_points.size(); ++point)
    {
        const dii::Tensor<1, dim> &n = normals[point];

        // Precompute tractions for all DOFs at this quadrature point
        std::vector<dii::Tensor<1, dim>> traction(dofs_per_cell);

        for (unsigned int j = 0; j < dofs_per_cell; ++j)
           {
            
             // Get strain on minus side (false = cell)
             auto eps = scratch_data.bstrain_function(fe_face, j, point);
        
             // Traction on minus side: C:ε⁻ · n (outward from cell)
             dii::Tensor<1, dim> t_minus = scratch_data.stiffness_tensor * eps * n;
        
             // Average traction: {t} = 0.5 * (t⁻ + t⁺)
             traction[j] = t_minus;
           }

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int comp_i = fe_face.get_fe().system_to_component_index(i).first;
              const unsigned int mcomp_i = comp_i % dim;
              const double shape_i = fe_face.shape_value_component(i, point, comp_i);
              
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                  {
                    const unsigned int comp_j = fe_face.get_fe().system_to_component_index(j).first;
                    const unsigned int mcomp_j = comp_j % dim;
                    const double shape_j = fe_face.shape_value_component(j, point, comp_j);

                    double penalty_term = 0.0;
                    if (comp_i == comp_j)
                        penalty_term = eta * shape_i * shape_j;

                    copy_data.cell_matrix(i, j) += 
                    (
                     // Term 1: {t_j} · [v_i]
                     -traction[j][mcomp_i] * shape_i
                  
                    // Term 2: penalty term
                     + penalty_term

                     + S*traction[i][mcomp_j] * shape_j
                                
                    ) * JxW[point];
                  }

            }
    }
     
  }



template <int dim>
 struct Mw_helper<dim>::ScratchData
 {
  using InterfaceStrainFunction = 
        std::function<dii::SymmetricTensor<2,dim>(const dii::FEInterfaceValues<dim>&, 
                                                  const bool &,
                                                  const unsigned int &, 
                                                  const unsigned int &,
                                                  const unsigned int & )>;
  using BoundaryStrainFunction = 
        std::function<dii::SymmetricTensor<2,dim>(const dii::FEFaceValuesBase<dim>&, 
                                                  const unsigned int &, 
                                                  const unsigned int & )>;

  ScratchData(   const dii::hp::FECollection<dim>  &fe_collection,
            const dii::hp::QCollection<dim - 1> &face_qcollection,           
            const dii::UpdateFlags  &face_flags,
            double vp_, double rho_, double R_, double R_mix_, int p_degree_,
            unsigned int cg_id_,unsigned int dg1_id_, unsigned int dg2_id_,
            bool fracture_,
            dii::SymmetricTensor<4,dim> stiffness_tensor_,
            dii::SymmetricTensor<2,dim> fracture_stiffness_,
            InterfaceStrainFunction fstrain_function_,
            BoundaryStrainFunction bstrain_function_
            );

   ScratchData(const ScratchData &scratch_data);
   
   dii::FEInterfaceValues<dim> fe_interface_values;
   const double vp;
   const double rho;
   const double R;
   const double R_mix;
   const int p_degree;
   const unsigned int cg_id, dg1_id, dg2_id;
   const bool fracture;
   dii::SymmetricTensor<4,dim> stiffness_tensor;
   dii::SymmetricTensor<2,dim> fracture_stiffness;
   InterfaceStrainFunction fstrain_function;
   BoundaryStrainFunction bstrain_function;
 };
 
template <int dim>
 Mw_helper<dim>::ScratchData::ScratchData(const dii::hp::FECollection<dim>  &fe_collection,
                                          const dii::hp::QCollection<dim - 1> &face_qcollection,           
                                          const dii::UpdateFlags  &face_flags,
                                          double vp_, double rho_, double R_, double R_mix_, int p_degree_,
                                          unsigned int cg_id_,unsigned int dg1_id_, unsigned int dg2_id_,
                                          bool fracture_,
                                          dii::SymmetricTensor<4,dim> stiffness_tensor_,
                                          dii::SymmetricTensor<2,dim> fracture_stiffness_,
                                          typename Mw_helper<dim>::ScratchData::InterfaceStrainFunction fstrain_function_,
                                          typename Mw_helper<dim>::ScratchData::BoundaryStrainFunction bstrain_function_
                                         ):

    fe_interface_values( fe_collection, face_qcollection, face_flags),
    vp(vp_),
    rho(rho_),
    R(R_),
    R_mix(R_mix_),
    p_degree(p_degree_),
    cg_id(cg_id_),
    dg1_id(dg1_id_),
    dg2_id(dg2_id_),
    fracture(fracture_),
    stiffness_tensor(stiffness_tensor_),
    fracture_stiffness(fracture_stiffness_),
    fstrain_function(fstrain_function_),
    bstrain_function(bstrain_function_)
    {}

template <int dim>
 Mw_helper<dim>::ScratchData::ScratchData(const  Mw_helper<dim>::ScratchData &scratch_data):
    fe_interface_values(scratch_data.fe_interface_values.get_fe_collection(),
                        scratch_data.fe_interface_values.get_quadrature_collection(),
                        scratch_data.fe_interface_values.get_update_flags()),
    vp(scratch_data.vp),
    rho(scratch_data.rho),
    R(scratch_data.R),
    R_mix(scratch_data.R_mix),
    p_degree(scratch_data.p_degree),
    cg_id(scratch_data.cg_id),
    dg1_id(scratch_data.dg1_id),
    dg2_id(scratch_data.dg2_id),
    fracture(scratch_data.fracture),
    stiffness_tensor(scratch_data.stiffness_tensor),
    fracture_stiffness(scratch_data.fracture_stiffness),
    fstrain_function(scratch_data.fstrain_function),
    bstrain_function(scratch_data.bstrain_function)

    {}


template <int dim>
struct Mw_helper<dim>::CopyDataFace
  {
    dii::FullMatrix<double> cell_matrix;
    std::vector<dii::types::global_dof_index> joint_dof_indices;
    std::array<double, 2>  values;
    std::array<unsigned int, 2>  cell_indices;
  };

template <int dim>
struct Mw_helper<dim>::CopyData
  {
    std::vector<dii::types::global_dof_index> local_dof_indices;
    std::vector<CopyDataFace>  face_data;
    dii::FullMatrix<double> cell_matrix;
    double value;
    unsigned int cell_index;
    template <class Iterator>
    void reinit(const Iterator &cell, unsigned int dofs_per_cell);
  };

template <int dim>
  template <class Iterator>
  void Mw_helper<dim>::CopyData::reinit(const Iterator &cell, unsigned int dofs_per_cell)
  {
    cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
    local_dof_indices.resize(dofs_per_cell);
    cell->get_dof_indices(local_dof_indices);
  }

  template <int dim>
  InitialValuesU<dim>::InitialValuesU():
  dii::Function<dim>(2*dim) // considering FE_Nothing components
  {}

template <int dim>
  void InitialValuesU<dim>::vector_value(const dii::Point<dim> & /*p*/,
                                        dii::Vector<double> &values) const 
  {
    AssertDimension(values.size(), 2*dim);
    values = 0.0;

  }

template <int dim>
  InitialValuesV<dim>::InitialValuesV():
  dii::Function<dim>(2*dim)
  {}

template <int dim>
  void InitialValuesV<dim>::vector_value(const dii::Point<dim> & /*p*/,
                                        dii::Vector<double> &values) const 
  {
    AssertDimension(values.size(), 2*dim);
    values = 0.0;

  }


template <int dim>
dii::SymmetricTensor<4, dim> get_stiffness_tensor(const double & vp,
                                                  const double & vs,
                                                  const double & rho)
 {
   const double mu =vs*vs*rho;
   const double lambda =rho*vp*vp - 2*mu;

   dii::SymmetricTensor<4, dim> tmp;
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          for (unsigned int k = 0; k < dim; ++k)
            for (unsigned int l = 0; l < dim; ++l)
              tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) +
                                 ((i == l) && (j == k) ? mu : 0.0) +
                                 ((i == j) && (k == l) ? lambda : 0.0));
      return tmp;
 }

template <int dim>
dii::SymmetricTensor<2, dim> get_strain(const dii::FEValues<dim> & fe_values,
                                        const unsigned int & shape_func,
                                        const unsigned int  & q_point)
  {
    dii::SymmetricTensor<2, dim> tmp;
    
    // Get the actual component this shape function belongs to
    const unsigned int comp = fe_values.get_fe().system_to_component_index(shape_func).first;
    // Map to physical space [0, dim) for strain tensor indexing
    const unsigned int mcomp = comp % dim;
    
    // Get the full gradient (all spatial derivatives)
    const dii::Tensor<1, dim> grad_phi = fe_values.shape_grad(shape_func, q_point);
    
    // Diagonal strain: ε_ii = ∂u_i/∂x_i
    tmp[mcomp][mcomp] = grad_phi[mcomp];
    
    // Off-diagonal strains: ε_ij = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)
    for (unsigned int i = 0; i < dim; ++i)
    {
      if (i != mcomp)
        tmp[mcomp][i] = tmp[i][mcomp] = 0.5 * grad_phi[i];
    }
           
    return tmp;
  }

  template <int dim>
  dii::SymmetricTensor<2, dim> get_strain2(const dii::FEValues<dim> &fe_values,
                                     const unsigned int & shape_func,
                                     const unsigned int  & q_point)
  {
     dii::SymmetricTensor<2, dim> tmp;
     for (unsigned int i = 0; i < dim; ++i)
        tmp[i][i] = fe_values.shape_grad_component(shape_func, q_point, i)[i];
     
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = i + 1; j < dim; ++j)
          tmp[i][j] =
            (fe_values.shape_grad_component(shape_func, q_point, i)[j] +
             fe_values.shape_grad_component(shape_func, q_point, j)[i]) /
            2;
  
      return tmp;
  }

  template <int dim>
  dii::SymmetricTensor<2, dim> get_strain_comp(const dii::FEInterfaceValues<dim> &face_values,
                                               const bool & here_there,
                                               const unsigned int & shape_index,
                                               const unsigned int  & q_point,
                                               const unsigned int comp)
  {
    dii::SymmetricTensor<2,dim> eps;

    // Gradient of the comp component of the shape function
    const dii::Tensor<1,dim> grad_phi =
        face_values.shape_grad(here_there, shape_index, q_point, comp);

   //The following maps any comp component  to 0...dim-1
   // This is relevant for the dg_fe_system, where DG components are numbered after dim-1:, dim,...2*dim-1
   const unsigned int mcomp = comp % dim;

    // Fill eps tensor
   eps[mcomp][mcomp] = grad_phi[mcomp];

   for (unsigned int i = 0; i < dim; ++i)
    {
     if (i != mcomp)
        eps[mcomp][i] = eps[i][mcomp] = 0.5 * grad_phi[i];
   }

  return eps;
}


template <int dim>
dii::SymmetricTensor<2, dim> get_strain_comp(const dii::FEFaceValuesBase<dim> & face_values,
                                             const unsigned int & shape_index,
                                             const unsigned int  & q_point)
  {
    dii::SymmetricTensor<2,dim> eps;
    const unsigned int comp = face_values.get_fe().system_to_component_index(shape_index).first;
    const dii::Tensor<1,dim> grad_phi =face_values.shape_grad(shape_index, q_point);

   //The following maps any comp component  to 0...dim-1
   // This is relevant for the dg_fe_system, where DG components are numbered after dim-1:, dim,...2*dim-1
   const unsigned int mcomp = comp % dim;

   // Fill eps tensor
   eps[mcomp][mcomp] = grad_phi[mcomp];

   for (unsigned int i = 0; i < dim; ++i)
    {
     if (i != mcomp)
        eps[mcomp][i] = eps[i][mcomp] = 0.5 * grad_phi[i];
   }


    return eps;
  }

template <int dim>
  dii::SymmetricTensor<2, dim> get_fracture_stiffness(const double & Zn,
                                                      const double & Zt)
{
  const double eps = 1e-30;
  const double inv_Zn = 1.0 / (std::abs(Zn) + eps);
  const double inv_Zt = 1.0 / (std::abs(Zt) + eps);

  dii::SymmetricTensor<2, dim> tmp;
  if constexpr (dim == 2)
  {
    tmp[0][0] = inv_Zt;
    tmp[0][1] = 0.0;
    tmp[1][1] = inv_Zn;
  }
  else
  {
    tmp[0][0] = inv_Zt;
    tmp[0][1] = tmp[0][2] = tmp[1][2] = 0.0;
    tmp[1][1] = inv_Zt;
    tmp[2][2] = inv_Zn;
  }
  return tmp;
}


template<int dim>
void cg_dg_cells(dii::DoFHandler<dim> &dof_handler, 
                 double & fcoord,
                 const int & rows_above,
                 const int & cg_id,
                 const int & dg1_id, 
                 const int & dg2_id,
                 double & top_dg_y,
                 double & bottom_dg_y)
{
    // Step 1: Find cell closest to fcoord to get h in the uniform region
    double min_dist = std::numeric_limits<double>::max();
    typename dii::DoFHandler<dim>::active_cell_iterator fracture_cell;

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      const double dist = std::abs(cell->center()[1] - fcoord);
      if (dist < min_dist)
      {
        min_dist = dist;
        fracture_cell = cell;
      }
    }

    double h = fracture_cell->extent_in_direction(1);

    // Step 2: Find the horizontal interior face closest to fcoord
    double best_dist = std::numeric_limits<double>::max();
    double face_y = fcoord;

    for (const auto &cell : dof_handler.active_cell_iterators())
      for (unsigned int f = 0; f < dii::GeometryInfo<dim>::faces_per_cell; ++f)
      {
        if (cell->face(f)->at_boundary())
          continue;
        // Skip vertical faces: only keep faces whose normal is predominantly in y
        const dii::Tensor<1, dim> normal = cell->face(f)->center() - cell->center();
        if (std::abs(normal[1]) < std::abs(normal[0]))
          continue;
        const double fy = cell->face(f)->center()[1];
        const double dist = std::abs(fy - fcoord);
        if (dist < best_dist)
        {
          best_dist = dist;
          face_y = fy;
        }
      }

    fcoord = face_y;

    // Step 3: DG zone symmetric around fcoord
    top_dg_y    = fcoord + rows_above * h;
    bottom_dg_y = fcoord - rows_above * h;
    const double tol = 1e-10;

    // Step 4: Classify cells
    int n_cg = 0, n_dg1 = 0, n_dg2 = 0;

    for (auto &cell : dof_handler.active_cell_iterators())
    {
      const double cy = cell->center()[1];

      if (cy >= bottom_dg_y - tol && cy <= top_dg_y + tol)
      {
        cell->set_active_fe_index(1);
        if (cy >= fcoord)
          { cell->set_material_id(dg1_id); ++n_dg1; }
        else
          { cell->set_material_id(dg2_id); ++n_dg2; }
      }
      else
      {
        cell->set_active_fe_index(0);
        cell->set_material_id(cg_id);
        ++n_cg;
      }
    }

    std::cout << "cg_dg_cells: fcoord = " << fcoord
              << ", top_dg_y = " << top_dg_y
              << ", bottom_dg_y = " << bottom_dg_y
              << ", CG: " << n_cg 
              << ", DG1: " << n_dg1 
              << ", DG2: " << n_dg2 << std::endl;
}

template <typename SparseMatrixType>
 void check_matrix_symmetry(const SparseMatrixType& matrix, 
                            double tol, 
                            bool verbose)
{
  double max_diff = 0.0;
  unsigned int max_i = 0, max_j = 0;
  double val_ij_at_max = 0.0;
  double val_ji_at_max = 0.0;


    for (const auto &entry : matrix)
    {
        const unsigned int i = entry.row();
        const unsigned int j = entry.column();
        const double val_ij = entry.value();
        const double val_ji = matrix(j, i);
        const double diff = std::abs(val_ij - val_ji);
        if (diff > max_diff)
        {
            max_diff = diff;
            max_i = i;
            max_j = j;
        }
    }
    if (max_diff <= tol)
       std::cout << "Symmetric" << std::endl;
    else
       std::cout << "Not symmetric" << std::endl;
    if (verbose)
       {
           std::cout << "Max asymmetry: " << max_diff << " at (" << max_i << "," << max_j << ")" << std::endl;
           const double scale = std::max(std::abs(val_ij_at_max), std::abs(val_ji_at_max));
           const double rel_asym = (scale > 1e-100) ? max_diff / scale : 0.0;
           std::cout << "  A(" << max_i << "," << max_j << ") = " << val_ij_at_max
                     << ",  A(" << max_j << "," << max_i << ") = " << val_ji_at_max
                     << ",  relative asymmetry = " << rel_asym << std::endl;
   
           // Count entries exceeding threshold
           unsigned int n_bad = 0;
           double threshold = 1e-15;
           for (const auto &entry : matrix)
           {
               if (entry.row() <= entry.column()) continue;
               double diff = std::abs(entry.value() - matrix(entry.column(), entry.row()));
               if (diff > threshold) ++n_bad;
           }
           std::cout << "Entries with |A(i,j)-A(j,i)| > " << threshold << ": " << n_bad << std::endl;
       }
   
}

template<int dim>
void check_conditioning(const dii::DoFHandler<dim> &dof_handler,
                        const dii::SparseMatrix<double> &A,
                        const std::string &name,
                        double &lambda_max,
                        double &lambda_min,
                        double &cond_est)
{
  std::cout << "Estimating conditioning for " << name << " (power iterations)..." << std::endl;

  dii::Vector<double> v(dof_handler.n_dofs());
  for (unsigned int i = 0; i < v.size(); ++i)
    v[i] = std::sin(0.37 * (i + 1));
  v /= v.l2_norm();

  // Power iteration for lambda_max (5 iterations)
  lambda_max = 0.0;
  for (int iter = 0; iter < 5; ++iter) {
    dii::Vector<double> Av(v.size());
    A.vmult(Av, v);
    lambda_max = Av.l2_norm();
    v = Av;
    v /= lambda_max;
  }

  // Inverse power for lambda_min (5 iterations with diagonal preconditioner as cheap inverse)
  v = 1.0;  // reset to constant
  v /= v.l2_norm();

  dii::Vector<double> diag(dof_handler.n_dofs());
  diag = 1.0;  // simplified: use identity for cheap diagonal inverse

  double lambda_min_inv = 0.0;
  for (int iter = 0; iter < 5; ++iter) {
    dii::Vector<double> tmp(v.size());
    A.vmult(tmp, v);
    dii::Vector<double> Dinv_Av(v.size());
    for (unsigned int i = 0; i < v.size(); ++i)
      Dinv_Av[i] = diag[i] * tmp[i];
    lambda_min_inv = Dinv_Av.l2_norm();
    v = Dinv_Av;
    v /= lambda_min_inv;
  }
  lambda_min = 1.0 / (lambda_min_inv + 1e-16);

  cond_est = lambda_max / (lambda_min + 1e-16);

  std::cout << "  λ_max~ " << lambda_max << ", λ_min~ " << lambda_min
            << ", κ~ " << cond_est << std::endl;
}

// Returns the active_fe_index (0=CG, 1=DG) for the given global DOF.
// For shared DOFs (CG), returns the FE index of the first cell found.
template <int dim>
unsigned int get_fe_index_for_dof(const dii::DoFHandler<dim> &dof_handler,
                                  unsigned int global_dof)
{
  for (const auto &cell : dof_handler.active_cell_iterators())
  {
    std::vector<dii::types::global_dof_index> local_dof_indices(cell->get_fe().n_dofs_per_cell());
    cell->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < local_dof_indices.size(); ++i)
      if (local_dof_indices[i] == global_dof)
        return cell->active_fe_index();  // 0 = CG, 1 = DG
  }
  return dii::numbers::invalid_unsigned_int;  // not found
}