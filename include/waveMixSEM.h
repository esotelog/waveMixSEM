/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2006 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2006
 */

#ifndef WAVE_MIX_SEM_H       
#define WAVE_MIX_SEM_H


#include <fstream>
#include <iostream>
#include <iomanip>
#include <exception>
#include <string>
#include <cmath>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/exceptions.h>       // for Assert and ExcMessage
#include <deal.II/base/point.h>  

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h> 

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_creator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_geometry.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/refinement.h>
#include <deal.II/hp/q_collection.h> 
#include <deal.II/fe/fe_series.h>
 
#include "functions.h"
//#include "Param_handler.h"


namespace dii = dealii;
namespace diiM=dii::MeshWorker;


template <int dim>
  class waveMixSEM
  {
  public:
    waveMixSEM(const Parameters & prm);
    void run();
  private:
    enum{
        cg_id,
        dg1_id,
        dg2_id
      };

    void make_mesh();
    void setup_system();
    void assemble_MKmatrices();
    void include_interface_fluxes();
    void assemble_boundary_matrix();
    void assemble_point_source();
    void Neumann_plane_wave();
    void timestep_loop();
    void output_results() const;

    void write_receiver_header(std::ofstream &out_yup,
                               std::ofstream &out_ydg,
                               std::ofstream &out_ydown);
    void write_receiver_data(const double time,
                             std::ofstream &out_yup,
                             std::ofstream &out_ydg,
                             std::ofstream &out_ydown);

    dii::Triangulation<dim> mesh; 
    dii::FESystem<dim> cg_fe;
    dii::FESystem<dim> dg_fe;
    dii::DoFHandler<dim>  dof_handler;
    dii::AffineConstraints<double> constraints;
    dii::hp::FECollection<dim> fe_collection;
    dii::hp::QCollection<dim>  qcollection;
    dii::hp::QCollection<dim - 1> face_qcollection;
    dii::MappingQ1<dim> mapping;

    dii::SparsityPattern sparsity_pattern;
    dii::SparsityPattern sparsity_pattern_volume;
    dii::SparseMatrix<double> stiffness_matrix;
    dii::SparseMatrix<double> boundary_matrix;
    dii::Vector<double> solution_u, solution_v;
    dii::Vector<double> mass_vector, inverse_mass_vector;
    dii::Vector<double> old_solution_u, old_solution_v;
    dii::Vector<double> source_term;
  
    Parameters param;
    dii::SymmetricTensor<4,dim> stiffness_tensor;
    Source<dim> source_function;
   
    int timestep_number;
    std::vector<dii::Point<dim>> receivers_pos;
    typename dii::DoFHandler<dim>::active_cell_iterator source_cell;
    dii::Point<dim> source_pos_ref;
    Mesh_operations<dim> mesh_op;
    const unsigned int bid_min=31, bid_max=32;
    bool source_initialized = false;
    std::ofstream out;

    const unsigned int boundary_id_left = 0;
    const unsigned int boundary_id_right = 1;

  };


//#include "template_waveMix.tcc"
#include "template_waveMixSEM.tcc"
#endif