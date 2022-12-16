
#include "Parameters.h"

#include <iostream>

//////////////////////////////////////////////////////////
//               X M L   U t i l s                      //
//////////////////////////////////////////////////////////

//-------------------
// xml_util_get_bool
//-------------------
// Get a Boolean value from the XML document.
//
// A Boolean is 0/1 or true/false
//
//
bool xml_util_get_bool(tinyxml2::XMLElement* element, std::string& name, bool& value, bool required=false)
{
  bool bool_value;
  auto result = element->QueryBoolAttribute(name.c_str(), &bool_value);
  if (!required && (result == tinyxml2::XML_NO_ATTRIBUTE)) {
      return false;
  } else if (result != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("No '" + name + "' element found.");
  }

  if (bool_value) {  
    value = 1;
  } else {
    value = 0;
  }

  return true;
}

//---------------------
// xml_util_get_double
//---------------------
//
bool xml_util_get_double(tinyxml2::XMLElement* element, std::string& name, double& value, bool required=false)
{
  if (element->QueryDoubleAttribute(name.c_str(), &value) != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("No '" + name + "' element found.");
  }
}

//--------------------
// xml_util_get_float
//--------------------
//
bool xml_util_get_float(tinyxml2::XMLElement* element, std::string& name, float& value, bool required=false)
{
  if (element->QueryFloatAttribute(name.c_str(), &value) != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("No '" + name + "' element found.");
  }
}

//------------------
// xml_util_get_int
//------------------
//
bool xml_util_get_int(tinyxml2::XMLElement* element, std::string& name, int& value, bool required=false)
{
  auto result = element->QueryIntAttribute(name.c_str(), &value);
  if (!required && (result == tinyxml2::XML_NO_ATTRIBUTE)) {
      return false;
  } else if (result != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("No '" + name + "' element found.");
  }
  return true;
}

//---------------------
// xml_util_get_string
//---------------------
//
bool xml_util_get_string(tinyxml2::XMLElement* element, std::string& name, std::string& value, bool required=true)
{
  const char* svalue;
  auto result = element->QueryStringAttribute(name.c_str(), &svalue);

  if (!required && (result == tinyxml2::XML_NO_ATTRIBUTE)) {
      return false;
  } else if (result != tinyxml2::XML_SUCCESS) {
      throw std::runtime_error("No '" + name + "' element found.");
  }

  value = std::string(svalue);
  return true;
}

//////////////////////////////////////////////////////////
//                   Parameters                         //
//////////////////////////////////////////////////////////
 
Parameters::Parameters() 
{
}

Parameters::~Parameters() 
{
}

void Parameters::get_logging_levels(int& verbose, int& warning, int& debug)
{ 
  verbose = general_parameters_.verbose_;
  warning = general_parameters_.warning_;
  debug = general_parameters_.debug_;
}

//----------------
// set_parameters
//----------------
//
void Parameters::set_parameters(std::string file_name)
{ 
  tinyxml2::XMLDocument doc;
  
  std::cout << "========== Parameters::set_parameters ==========" << std::endl;
  std::cout << "[Parameters::set_parameters] file_name: " << file_name << std::endl;

  doc.LoadFile(file_name.c_str());
  
  auto fsi_file = doc.FirstChildElement(FSI_FILE.c_str());
  if (fsi_file == nullptr) {
      throw std::runtime_error("No '" + FSI_FILE + "' element found.");
  }

  // Get general parameters.
  general_parameters_.set_parameters(fsi_file);
  
  // Get solution parameters.
  solution_parameters_.set_parameters(fsi_file);
  
  // Get mesh parameters.
  std::cout << "Meshes ..." << std::endl;
  auto xml_mesh_params = fsi_file->FirstChildElement(MESH_PARAMS.c_str());
  while (xml_mesh_params != nullptr) {
      MeshParameters mesh_params;
      mesh_params.set_parameters(xml_mesh_params);
      mesh_parameters_.emplace_back(mesh_params);
      xml_mesh_params = xml_mesh_params->NextSiblingElement(MESH_PARAMS.c_str());
  }

  // Get equation parameters.
  std::cout << "Equations ..." << std::endl;
  auto xml_eq_params = fsi_file->FirstChildElement(EQUATION_PARAMS);
  while (xml_eq_params != nullptr) {
      EquationParameters equation_params;
      equation_params.set_parameters(xml_eq_params);
      equation_parameters_.emplace_back(equation_params);
      xml_eq_params = xml_eq_params->NextSiblingElement(EQUATION_PARAMS);
  }

}


//////////////////////////////////////////////////////////
//            BoundaryConditionParameters               //
//////////////////////////////////////////////////////////

//----------------
// set_parameters
//----------------
//
void BoundaryConditionParameters::set_parameters(tinyxml2::XMLElement* bc_params)
{ 
  std::cout << std::endl; 
  std::cout << "Boundary condition parameters element ..." << std::endl;
  
  xml_util_get_string(bc_params, TYPE, type_);
  std::cout << TYPE + ": " << type_ << std::endl;
  
  xml_util_get_string(bc_params, TIME_DEP, time_dependence_);
  std::cout <<  TIME_DEP + ": " << time_dependence_ << std::endl;
  
  xml_util_get_int(bc_params, VALUE, value_);
  std::cout <<  VALUE + ": " << value_ << std::endl;
  
  bool required = false;
  int ivalue;
  if (xml_util_get_int(bc_params, ZERO_OUT_PERIM, ivalue, required)) {
      zero_out_perimeter_ = ivalue;
      std::cout <<  ZERO_OUT_PERIM + ": " << zero_out_perimeter_ << std::endl;
  }

}

//////////////////////////////////////////////////////////
//                  EquationParameters                  //
//////////////////////////////////////////////////////////

//----------------
// set_parameters
//----------------
//
void EquationParameters::set_parameters(tinyxml2::XMLElement* equation_params)
{
  using namespace tinyxml2;
  std::cout << std::endl;
  std::cout << "Equation parameters element ..." << std::endl;

  xml_util_get_string(equation_params, PHYSICS, physics_);
  std::cout << PHYSICS + ": " << physics_ << std::endl;

  xml_util_get_float(equation_params, BACKFLOW_COEF, backflow_stab_coef_);
  std::cout << BACKFLOW_COEF + ": " << backflow_stab_coef_ << std::endl;

  xml_util_get_int(equation_params, COUPLED, coupled_);
  std::cout << COUPLED + ": " << coupled_ << std::endl;

  auto bc_params = equation_params->FirstChildElement(BOUNDARY_COND);
  while (bc_params != nullptr) {
      BoundaryConditionParameters boundary_condition;
      boundary_condition.set_parameters(bc_params);
      bc_params = bc_params->NextSiblingElement(BOUNDARY_COND);
      boundary_conditions_.emplace_back(boundary_condition);
  }
}

//////////////////////////////////////////////////////////
//                  GeneralParameters                   //
//////////////////////////////////////////////////////////

GeneralParameters::GeneralParameters()
{
  debug_ = 0;
  verbose_ = 0;
  warning_ = 0;
  save_results_using_vtk_format_ = 1;
}

//----------------
// set_parameters
//----------------
//
void GeneralParameters::set_parameters(tinyxml2::XMLElement* fsi_file)
{
  using namespace tinyxml2;
  std::cout << std::endl;
  std::cout << "GeneralParameters element ..." << std::endl;
  auto general_params = fsi_file->FirstChildElement(GENERAL_PARAMS.c_str());
  if (general_params == nullptr) {
      throw std::runtime_error("No '" + GENERAL_PARAMS + "' element found.");
  }

  xml_util_get_bool(general_params,   DEBUG, debug_);
  xml_util_get_int(general_params,    INCREMENT_IN_SAVING_RESTART_FILES, restart_save_increment_);
  xml_util_get_string(general_params, NAME_PREFIX_OF_SAVED_VTK_FILES, vtk_file_name_prefix_);
  xml_util_get_bool(general_params,   SAVE_AVERAGED_RESULTS, save_averaged_results_);
  xml_util_get_string(general_params, SAVE_RESULTS_IN_FOLDER, save_results_dir_);
  xml_util_get_bool(general_params,   SAVE_RESULTS_TO_VTK_FORMAT, save_results_using_vtk_format_);
  xml_util_get_string(general_params, SEARCHED_FILE_NAME_TO_TRIGGER_STOP, trigger_stop_file_name_);
  xml_util_get_bool(general_params,   SIMULATION_REQUIRES_REMESHING, simulation_requires_remeshing_);
  xml_util_get_int(general_params,    START_SAVING_AFTER_TIME_STEP, start_saving_time_step_);
  xml_util_get_bool(general_params,   VERBOSE, verbose_);
  xml_util_get_bool(general_params,   WARNING, warning_);

}

//////////////////////////////////////////////////////////
//             M e s h P a r a m e t e r s              //
//////////////////////////////////////////////////////////

//----------------
// set_parameters
//----------------
//
void MeshParameters::set_parameters(tinyxml2::XMLElement* mesh_params)
{
  using namespace tinyxml2;
  std::cout << std::endl;
  std::cout << "Mesh parameters element ..." << std::endl;

  xml_util_get_string(mesh_params, NAME, name_);
  std::cout << NAME + ": " << name_ << std::endl;

  //xml_util_get_string(mesh_params, PHYSICS, physics_);
  //std::cout << PHYSICS + ": " << physics_ << std::endl;

  xml_util_get_int(mesh_params, DOMAIN, domain_);
  std::cout << DOMAIN + ": " << domain_ << std::endl;

  xml_util_get_string(mesh_params, MESH_FILE, mesh_file_);
  std::cout << MESH_FILE + ": " << mesh_file_ << std::endl;

  std::cout << "  Faces: " << std::endl;
  auto faces = mesh_params->FirstChildElement(FACE);
  while (faces != nullptr) {
      const char* type;
      faces->QueryStringAttribute(FACE_TYPE, &type);
      std::cout << "    Type: " << type << std::endl;
      const char* file;
      faces->QueryStringAttribute(FACE_FILE, &file);
      std::cout << "    File: " << file << std::endl;
      faces = faces->NextSiblingElement(FACE);
      faces_.push_back({std::string(type),  std::string(file)});
  }

}

//////////////////////////////////////////////////////////
//                 SolutionParameters                   //
//////////////////////////////////////////////////////////

SolutionParameters::SolutionParameters()
{
  std::cout << "+ + + SolutionParameters() ctor + + + " << std::endl;
  cont_prev_sim_ = 0;
  num_dims_ = 0;
  num_time_steps_ = 0;
  num_init_time_steps_ = 10;
  spectral_radius_inf_time_step_ = 0.0;
  time_step_ = 0.0;
}

SolutionParameters::~SolutionParameters()
{
  std::cout << "- - - ~SolutionParameters() dtor - - - " << std::endl;
}

//----------------
// set_parameters
//----------------
//
void SolutionParameters::set_parameters(tinyxml2::XMLElement* fsi_file)
{
  using namespace tinyxml2;
  std::cout << std::endl;
  std::cout << "SolutionParameters element ..." << std::endl;
  auto simulation_params = fsi_file->FirstChildElement(SOLUTION_PARAMS.c_str());
  if (simulation_params == nullptr) {
      throw std::runtime_error("No '" + SOLUTION_PARAMS + "' element found.");
  }

  xml_util_get_int(simulation_params, CONT_PREV_SIM, cont_prev_sim_);
  std::cout << CONT_PREV_SIM + ": " << cont_prev_sim_ << std::endl;

  xml_util_get_int(simulation_params, NUM_DIMS, num_dims_);
  std::cout << NUM_DIMS + ": " << num_dims_ << std::endl;

  xml_util_get_int(simulation_params, NUM_INIT_TIME_STEPS, num_init_time_steps_);
  std::cout << NUM_INIT_TIME_STEPS + ": " << num_init_time_steps_ << std::endl;

  xml_util_get_int(simulation_params, NUM_TIME_STEPS, num_time_steps_);
  std::cout << NUM_TIME_STEPS + ": " << num_time_steps_ << std::endl;

  xml_util_get_double(simulation_params, SPECTRAL_RADIUS_OF_INFINITE_TIME_STEP, spectral_radius_inf_time_step_);
  std::cout << SPECTRAL_RADIUS_OF_INFINITE_TIME_STEP + ": " << spectral_radius_inf_time_step_ << std::endl;

  xml_util_get_double(simulation_params, TIME_STEP, time_step_);
  std::cout << TIME_STEP + ": " << time_step_ << std::endl;

}

