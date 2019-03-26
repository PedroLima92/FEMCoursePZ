# @file FEMCoursePZ/tinyxml/CMakeLists.txt  -- CMake file for the TinyXML module

file(GLOB sources *.cpp)

file(GLOB headers *.h)

set(FEMCoursePZ_tinyxml ${sources} ${headers} PARENT_SCOPE)

