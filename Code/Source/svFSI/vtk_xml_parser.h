
#include "ComMod.h"

#ifndef VTK_XML_PARSER
#define VTK_XML_PARSER

namespace vtk_xml_parser {

enum class VtkFileFormat { VTP, VTU };

class VtkFileExtentions {
 public:
  const static std::string VTK_VTU_EXTENSION;
  const static std::string VTK_VTP_EXTENSION;
};

void load_fiber_direction_vtu(const std::string& file_name,
                              const std::string& data_name, const int idx,
                              const int nsd, mshType& mesh);

void load_vtp(const std::string& file_name, faceType& face);

void load_vtp(const std::string& file_name, mshType& mesh);

void load_vtu(const std::string& file_name, mshType& mesh);

};  // namespace vtk_xml_parser

#endif
