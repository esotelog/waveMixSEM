#ifndef FUNCTIONS_H       
#define FUNCTIONS_H

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <exception>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>

#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/exceptions.h> 
#include <deal.II/base/exception_macros.h>        // for Assert and ExcMessage
#include <deal.II/base/point.h>  
#include <deal.II/base/numbers.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/scratch_data.h>

#include "Param_handler.h"
namespace dii = dealii;
namespace diiM = dii::MeshWorker;


template <int dim>
  class Source: public dii::Function <dim>
  {
   private:
    Parameters param;
    const double fo;
    const double source_rad;
    const int stf_type;
    const double M11, M12, M22, M21;
    const double loc_x, loc_y, loc_z;
    
    double source_time (const double & t) const;
    double source_space (const dii::Point<dim> & p) const;

  public:
   Source(const Parameters & param_);
   virtual double value(const dii::Point<dim> & p=dii::Point<dim>(),
                         const unsigned int component = 0) const override;
    
   double source_duration() const;

   std::vector<typename dii::DoFHandler<dim>::active_cell_iterator>
   find_source_cells(const dii::DoFHandler<dim> &dof_handler) const;

    dii::SymmetricTensor<2, dim> get_MT() const;
    dii::Tensor<1, dim> get_force_vector() const;
    bool is_isotropic() const;
    static bool is_isotropic(const Parameters & prm);
    static dii::Point<dim> get_source_position(const Parameters & prm);
    dii::Point<dim> source_pos;
   
  };

template <int dim>
 class Mesh_operations
 {
    private:
    std::vector<dii::Point<dim>> get_boundary_points() const;
    double get_velocity() const;
    const Parameters param;
    const double vs,vp;
    const int wavesamples, sourcesampling;
    const double fo, source_rad;
    const std::vector<dii::Point<dim>> boundary_points;
    double hbg, outer_radius;

    public:
    Mesh_operations(const Parameters & param_);
    void mesh_create(dii::Triangulation<dim> & mesh) const ;
    void label_extreme_boundary_faces(dii::DoFHandler<dim> &dof_handler,
                                      const unsigned int boundary_id, const unsigned int bid_min,
                                      const unsigned int bid_max, const unsigned int axis, // 0 for x, 1 for y
                                      const double fraction_nf, double &min, double &Lmax,
                                      std::vector<std::pair<typename dii::DoFHandler<dim>::cell_iterator,
                                      unsigned int>> &faces) const;

    void source_refine(dii::Triangulation<dim> & mesh);
    void refine_dg_rows(dii::Triangulation<dim> & mesh,
                       const double & top_dg_y,
                       const double & bottom_dg_y);

    void print_mesh(dii::Triangulation<dim> & mesh) const;
    double get_fcoord() const; 
    std::vector<dii::Point<dim>> receiver_coordinates(const double & h,
                                                      const double & yfrac, 
                                                      const double & bottom_dg_y,
                                                      const unsigned int & ncells,
                                                      const std::vector<double> & receivers_x
                                                    ) const;
    const double Lx, Ly;
 };

 template <int dim>
  class Mw_helper
  {
    public:
    
    struct ScratchData;
    struct CopyDataFace;
    struct CopyData;
    
    static void face_flux_worker (const typename dii::DoFHandler<dim>::cell_iterator & cell,
                                  const unsigned int & f,
                                  const unsigned int & sf,
                                  const typename dii::DoFHandler<dim>::cell_iterator & ncell,
                                  const unsigned int & nf,
                                  const unsigned int & nsf,
                                  ScratchData & scratch_data,
                                  CopyData & copy_data);
    
    static void boundary_worker (const typename dii::DoFHandler<dim>::cell_iterator & cell,
                                 const unsigned int & face_n,
                                 ScratchData & scratch_data,
                                 CopyData & copy_data);
  };

 
template <int dim>
  class InitialValuesU : public dii::Function<dim>
  {
    
  public:
    InitialValuesU();
    virtual void vector_value(const dii::Point<dim> & /*p*/,
                              dii::Vector<double> & value) const override;
  };

template <int dim>
  class InitialValuesV : public dii::Function<dim>
  {
  public:
    InitialValuesV();
    virtual void vector_value(const dii::Point<dim> & /*p*/,
                              dii::Vector<double> & value) const override;
  };

template <int dim>
dii::SymmetricTensor<4, dim> get_stiffness_tensor(const double & vp,
                                                  const double & vs,
                                                  const double & rho);


template <int dim>
dii::SymmetricTensor<2, dim> get_strain(const dii::FEValues<dim> & fe_values,
                                   const unsigned int & shape_func,
                                   const unsigned int  & q_point);
template <int dim>
dii::SymmetricTensor<2, dim> get_strain2(const dii::FEValues<dim> & fe_values,
                                   const unsigned int & shape_func,
                                   const unsigned int  & q_point);

template <int dim>
dii::SymmetricTensor<2, dim> get_strain_comp(const dii::FEInterfaceValues<dim> & face_values,
                                             const bool & here_there,
                                             const unsigned int & shape_index,
                                             const unsigned int  & q_point);
template <int dim>
dii::SymmetricTensor<2, dim> get_strain_comp(const dii::FEFaceValuesBase<dim> & face_values,
                                             const unsigned int & shape_index,
                                             const unsigned int  & q_point);

template<int dim>
void cg_dg_cells(dii::DoFHandler<dim> &dof_handler, 
                 double & fcoord,
                 const int & rows_above,
                 const int & cg_id,
                 const int & dg1_id, 
                 const int & dg2_id,
                 double & top_dg_y,
                 double & bottom_dg_y);     
                 
template <int dim>
dii::SymmetricTensor<2, dim> get_fracture_stiffness(const double & Zn,
                                                    const double & Zt);

template <typename SparseMatrixType>
void check_matrix_symmetry(const SparseMatrixType& matrix, 
                           double tol = 1e-12, 
                           bool verbose = false);

template<int dim>
void check_conditioning(const dii::DoFHandler<dim> &dof_handler,
                        const dii::SparseMatrix<double> &A,
                        const std::string &name,
                        double &lambda_max,
                        double &lambda_min,
                        double &cond_est);

template <int dim>
unsigned int get_fe_index_for_dof(const dii::DoFHandler<dim> &dof_handler,
                                  unsigned int global_dof);

double calc_time_interval(const double & space,
                         const double & vel,
                         const unsigned int & p_degree);

void rounding_digits(double & value, 
                     const int & digits, 
                     bool roundown_flag=true);

void echo_parameters(const Parameters &params);


#include "template_functions.tcc"

#endif