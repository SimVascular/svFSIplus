
// The functions defined here replicate the Fortran functions defined in LOADMSH.f.

#include "load_msh.h"

#include "consts.h"
#include "nn.h"
#include "read_msh.h"
#include "vtk_xml.h"

#include <iostream>
#include <fstream>
#include <sstream>

namespace load_msh {

#define ndbg_load_msh

//-----------
// read_ccne
//-----------
// Read mesh position coordinates and connectivity.
//
// [NOTE] Not imiplemented.
//
void read_ccne(Simulation* simulation, mshType& mesh, const MeshParameters* mesh_param)
{
  auto mesh_path = mesh_param->mesh_file_path();
  auto mesh_name = mesh_param->name();
  throw std::runtime_error("[read_ccne] read_ccne() is not implemented."); 
}

//-------------
// read_ndnlff
//-------------
// Read list of end nodes from a file into the face data structure.
//
// This sets the face data: 
//   face.nNo 
//   face.gN 
//   face.nEl 
//   face.eNoN 
//
void read_ndnlff(const std::string& file_name, faceType& face)
{
  std::ifstream end_nodes_file;
  end_nodes_file.open (file_name);
  if (!end_nodes_file.is_open()) {
    throw std::runtime_error("Failed to open the end nodes face file '" + file_name + "'.");
  }

  // Read the integer node IDs.
  //
  int node_id;
  std::string str_value;
  std::vector<int> node_ids;

  while (end_nodes_file >> str_value) {
    std::istringstream str_stream(str_value);
    if (!(str_stream >> node_id)) {
      throw std::runtime_error("Incorrect integer value '" + str_value + "' found in end nodes face file '" + file_name + "'."); 
    }
    node_ids.push_back(node_id);
  }

  if (node_ids.size() == 0) { 
    throw std::runtime_error("Failed to read the end nodes face file '" + file_name + "'.");
  }

  // Set the face node IDs.
  //
  face.nNo = node_ids.size();
  face.gN = Vector<int>(face.nNo);
  for (int i = 0; i < face.nNo; i++) { 
    face.gN[i] = node_ids[i] - 1; 
    //std::cout << "[read_ndnlff] a Ac: " << i+1 << " " << node_ids[i] << std::endl;
  }

  // Set the face element properties.
  face.nEl = 0;
  face.eNoN = 1;
}

//---------
// read_sv
//---------
// Create data for a mesh.
//
// Replicates Fortran READSV subroutine defined in LOADMSH.f.
//
//   SUBROUTINE READSV(list, lM)
//
void read_sv(Simulation* simulation, mshType& mesh, const MeshParameters* mesh_param)
{
  auto mesh_path = mesh_param->mesh_file_path();
  auto mesh_name = mesh_param->get_name();
  #define n_dbg_read_sv
  #ifdef dbg_read_sv
  DebugMsg dmsg(__func__, simulation->com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "Mesh name: " << mesh_name;
  dmsg << "Mesh path: " << mesh_path;
  #endif

  // Read in volume mesh.
  vtk_xml::read_vtu(mesh_path, mesh);

  // Set mesh element properites for the input element type.
  nn::select_ele(simulation->com_mod, mesh);

  // Check the mesh element node ordering.
  //
  // Note: This may change element node ordering.
  //
  auto& com_mod = simulation->get_com_mod();
  if (com_mod.ichckIEN) {
    read_msh_ns::check_ien(simulation, mesh);
  }

  // Read face meshes. 
  //
  // Creates nodal coordinates and element connectivity data.
  //
  mesh.nFa = mesh_param->face_parameters.size();
  mesh.fa.resize(mesh.nFa);

  if (mesh.lFib && (mesh.nFa > 1)) {
    throw std::runtime_error("There are " + std::to_string(mesh.nFa) + " faces defined for the '" +
        mesh.name + "' mesh. Only one face is allowed for a 1D fiber-based mesh.");
  }

  for (int i = 0; i < mesh.nFa; i++) {
    auto face_param = mesh_param->face_parameters[i];
    auto& face = mesh.fa[i];
    face.name = face_param->name();
    #ifdef dbg_read_sv
    dmsg << "face.name: " << face.name;
    #endif

    face.qmTRI3 = face_param->quadrature_modifier_TRI3();
    if (face.qmTRI3 < (1.0/3.0) || face.qmTRI3 > 1.0) {
      throw std::runtime_error("Quadrature_modifier_TRI3 must be in the range [1/3, 1].");
    }

    if (mesh.lFib) {
      auto face_path = face_param->end_nodes_face_file_path();
      #ifdef dbg_read_sv
      dmsg << "Read end nodes face file ... " << " ";
      dmsg << "face_path: " << face_path;
      #endif
      if (face_path == "") { 
        throw std::runtime_error("No end nodes face file path provided.");
      }

      read_ndnlff(face_path, face);

    } else {
      auto face_path = face_param->face_file_path();
      vtk_xml::read_vtp(face_path, face);

      // If node IDs were not read then create them.
      if (face.gN.size() == 0) {
        read_msh_ns::calc_nbc(mesh, face);

        // Reset the connecttivity with the new node IDs?
        for (int e = 0; e < face.nEl; e++) {
          for (int a = 0; a < face.eNoN; a++) {
            int Ac = face.IEN(a,e);
            Ac = face.gN(Ac);
            face.IEN(a,e) = Ac;
          }
        }
      }
    }

    nn::select_eleb(simulation, mesh, face);
  }

}

};

