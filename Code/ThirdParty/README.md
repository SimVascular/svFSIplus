This directory contains the source for external applications used by svFSIplus. The applications are compiled and liked with svMultiPhysics. 

Each application has its own license agreement.

------------
Applications
------------

eigen - A header-only application for matrix objects.

gklib_internal - A library used by METIS and ParMETIS applications.

metis_internal - The METIS mesh partitioning application used by ParMETIS.

parmetis_internal - The ParMETIS parallel mesh partitioning application.

tetgen - A mesh generatin application.

tinyxml - A header-only application used to read and write XML files.

---------
IMPORTANT
---------

The *_internal directory names must agree with the names given in 

 Code/CMake/SimVascularInternals.cmake
 Code/CMake/SimVascularThirdParty.cmake 

They are also referenced in ThirdParty and the solver/CMakeLists.txt file.


