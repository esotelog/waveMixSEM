#include "Param_handler.h"

Parameters::Parameters():

    vp(3000.0),
    vs(2000.0),
    rho(1600.0),
    fracture(false),
    Zn(1.0e-10),
    Zt(1.0e-10),


    //source parameters
    stf_type(1),
    fo(30.0),
    source_rad(2.5),
    source(2),
    M11(1.0),
    M12(0.0),
    M22(1.0),
    M21(0.0),
    phi(0.0),
    taper(0.1),

    loc_x(0.5),
    loc_y(0.5),
    loc_z(0.5),// in 2D this is not used

    // time discretization
   
    total_time(0.5), // (total_time>0; total time = source duration)
    dt_factor(0.5),
    theta(0.5),

    //mesh paramenters
    nsamples(12), /*number of samples per min wavelength*/
    Lx(100.0), /*length of the domain in x direction*/
    Ly(100.0), /*length of the domain in y direction*/
    sourcesampling(10),/* diameter of source will be sampled 6 times*/

    //fem
    dim(2),
    p_degree(2),
    penalty(20),
    mix_penalty(20),
    dg_rows(2),
    refine_dg_rows(false),

    //output file name (vtk and trices)
    receivers_x{0.4, 0.5, 0.6},
    vtk(false),
    vtk_step(50),
    outputfile("waveMix")

    {}

    
/*ParameterReader class
 * ********************
 */

ParameterReader::ParameterReader(const std::string & parameter_file_):
  parameter_file(parameter_file_)
{
  declare_parameters();
  read_parameters();
  pass_parameters();
}

/*Declare Parameters
 * *****************
 */
void ParameterReader::declare_parameters()
{
   prm.enter_subsection("Medium");
   {
    prm.declare_entry("vp",
                      std::to_string(param.vp),
                      dii::Patterns::Double(),
                      "compressional velocity");

    prm.declare_entry("vs",
                      std::to_string(param.vs),
                      dii::Patterns::Double(),
                      "shear velocity");

    prm.declare_entry("rho",
                      std::to_string(param.rho),
                      dii::Patterns::Double(),
                      "density");

    prm.declare_entry("fracture",
                      (param.fracture ? "true" : "false"),
                      dii::Patterns::Bool(),
                      "enable fracture (linear slip)");

    prm.declare_entry("Zn",
                      std::to_string(param.Zn),
                      dii::Patterns::Double(),
                      "normal impedance");

    prm.declare_entry("Zt",
                      std::to_string(param.Zt),
                      dii::Patterns::Double(),
                      "tangential impedance");

  }
  prm.leave_subsection();

  prm.enter_subsection("Source");
   {

    prm.declare_entry("stf_type",
                      std::to_string(param.stf_type),
                      dii::Patterns::Integer(),
                      "central frequency");

    prm.declare_entry("fo",
                      std::to_string(param.fo),
                      dii::Patterns::Double(),
                      "central frequency");

    prm.declare_entry("source_rad",
                      std::to_string(param.source_rad),
                      dii::Patterns::Double(),
                      "source_radius");

    prm.declare_entry("source",
                      std::to_string(param.source),
                      dii::Patterns::Integer(1),
                      "source type: 1 = moment tensor, 2 = vector force");

    prm.declare_entry("M11",
                      std::to_string(param.M11),
                      dii::Patterns::Double(),
                      "M11");

    prm.declare_entry("M12",
                      std::to_string(param.M12),
                      dii::Patterns::Double(),
                      "M12");

    prm.declare_entry("M21",
                      std::to_string(param.M21),
                      dii::Patterns::Double(),
                      "M21");

    prm.declare_entry("M22",
                      std::to_string(param.M22),
                      dii::Patterns::Double(),
                      "M22");

    prm.declare_entry("phi",
                      std::to_string(param.phi),
                      dii::Patterns::Double(),
                      "vector source angle phi (rad), phi = 0 along +x");

    prm.declare_entry("taper",
                      std::to_string(param.taper),
                      dii::Patterns::Double(0.1),
                      "taper (>= 0.1)");

    prm.declare_entry("loc_x",
                      std::to_string(param.loc_x),
                      dii::Patterns::Double(),
                      "source x coordinate");

    prm.declare_entry("loc_y",
                      std::to_string(param.loc_y),
                      dii::Patterns::Double(),
                      "source y coordinate");

    prm.declare_entry("loc_z",
                      std::to_string(param.loc_z),
                      dii::Patterns::Double(),
                      "source z coordinate");          
  }
  prm.leave_subsection();

  prm.enter_subsection("Discretization");
   {
    prm.declare_entry("total_time",
              std::to_string(param.total_time),
              dii::Patterns::Double(),
              "total time factor");

    prm.declare_entry("dt_factor",
              std::to_string(param.dt_factor),
              dii::Patterns::Double(),
              "dt factor for implicit time step");

    prm.declare_entry("nsamples",
                      std::to_string(param.nsamples),
                      dii::Patterns::Integer(1),
                      "samples per wavelength");

    prm.declare_entry("Lx",
              std::to_string(param.Lx),
              dii::Patterns::Double(),
              "length of the domain in x direction");

    prm.declare_entry("Ly",
              std::to_string(param.Ly),
              dii::Patterns::Double(),
              "length of the domain in y direction");

    prm.declare_entry("sourcesampling",
                      std::to_string(param.sourcesampling),
                      dii::Patterns::Integer(1),
                      "samples per source radius");

   prm.declare_entry("dg_rows",
                      std::to_string(param.dg_rows),
                      dii::Patterns::Integer(2),
                      "DG rows (multiple of 2)");

   prm.declare_entry("refine_dg_rows",
                      (param.refine_dg_rows ? "true" : "false"),
                      dii::Patterns::Bool(),
                      "Refine DG rows (doubles initial DG rows)");
    
  }
  prm.leave_subsection();

  prm.enter_subsection("FEM");
  {
    prm.declare_entry("dim",
                      std::to_string(param.dim),
                      dii::Patterns::Integer(1),
                      "dimensionality");

    prm.declare_entry("p_degree",
                      std::to_string(param.p_degree),
                      dii::Patterns::Integer(1),
                      "shape function polynomial degree");
    
   prm.declare_entry("penalty",
                      std::to_string(param.penalty),
                      dii::Patterns::Double(),
                      "DG/DG penalty parameter");

   prm.declare_entry("mix_penalty",
                      std::to_string(param.mix_penalty),
                      dii::Patterns::Double(),
                      "CG/DG penalty parameter");
  }
   prm.leave_subsection();

   prm.enter_subsection("Output");

    {

      prm.declare_entry("receivers_x",
                        "0.5",
                        dii::Patterns::List(dii::Patterns::Double()),
                        "receivers x coordinates (comma-separated)");

      prm.declare_entry("vtk",
                        (param.vtk ? "true" : "false"),
                        dii::Patterns::Bool(),
                        "write VTK output");

      prm.declare_entry("vtk_step",
                        std::to_string(param.vtk_step),
                        dii::Patterns::Integer(1),
                        "number of time steps between VTK outputs");

      prm.declare_entry("outputfile",
                        param.outputfile,
                        dii::Patterns::Anything(),
                        "output file name (vtk and trices)");
    }

     prm.leave_subsection();
}

/*Read Parameters from file
 * ***********************
 */

void ParameterReader::read_parameters ()
{
  std::ifstream in(parameter_file);
  if (!in.is_open())
    throw std::runtime_error(std::string("Cannot open parameter file: ") + parameter_file);
  prm.parse_input (in, parameter_file);
}


/*Pass parameters to structure
 * ***********************
 */

void ParameterReader::pass_parameters()

{
  prm.enter_subsection("Medium");
     param.vp =prm.get_double("vp");
     param.vs =prm.get_double("vs");
     param.rho =prm.get_double("rho");
     param.fracture =prm.get_bool("fracture");
     param.Zn =prm.get_double("Zn");
     param.Zt =prm.get_double("Zt");
  prm.leave_subsection();

  prm.enter_subsection("Source");
     param.stf_type =prm.get_integer("stf_type");
     param.fo =prm.get_double("fo");
     param.source_rad =prm.get_double("source_rad");
     param.source =prm.get_integer("source");
     param.M11 =prm.get_double("M11");
     param.M12 =prm.get_double("M12");
     param.M21 =prm.get_double("M21");
     param.M22 =prm.get_double("M22");
     param.phi =prm.get_double("phi");
     param.taper =prm.get_double("taper");
     param.loc_x =prm.get_double("loc_x");
     param.loc_y =prm.get_double("loc_y");
     param.loc_z =prm.get_double("loc_z");
  prm.leave_subsection();

  prm.enter_subsection("Discretization");
     param.total_time    = prm.get_double("total_time");
     param.dt_factor     = prm.get_double("dt_factor");
     param.nsamples      = prm.get_integer("nsamples");
     param.Lx            = prm.get_double("Lx");
     param.Ly            = prm.get_double("Ly");
     param.sourcesampling = prm.get_integer("sourcesampling");
     param.dg_rows       = prm.get_integer("dg_rows");
     param.refine_dg_rows = prm.get_bool("refine_dg_rows");
  prm.leave_subsection();

  prm.enter_subsection("FEM");
     param.dim =prm.get_integer("dim");
     param.p_degree =prm.get_integer("p_degree");
      param.penalty =prm.get_double("penalty");
      param.mix_penalty =prm.get_double("mix_penalty");
  prm.leave_subsection();

  prm.enter_subsection("Output");

    //get receiver x coordinates
    const auto rx_list = dii::Utilities::split_string_list(prm.get("receivers_x"));
    param.receivers_x.clear();
    for (const auto &s : rx_list)
      param.receivers_x.push_back(dii::Utilities::string_to_double(s));

    param.vtk =prm.get_bool("vtk");
    param.vtk_step = prm.get_integer("vtk_step");
    param.outputfile =prm.get("outputfile");
  prm.leave_subsection();

}

Parameters ParameterReader::get_parameters()
{
  return param;
}
