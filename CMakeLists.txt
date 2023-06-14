cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(Point)

SET(VTK_DIR /home/angel/Documentos/Instaladores/VTK/build)

find_package(VTK COMPONENTS
  CommonColor
  CommonCore
  CommonDataModel
  FiltersSources
  IOLegacy
  InteractionStyle
  RenderingContextOpenGL2
  RenderingCore
  RenderingFreeType
  RenderingGL2PSOpenGL2
  RenderingOpenGL2
)

find_package(Eigen3)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Point: Unable to find the VTK build folder.")
endif()

# Prevent a "command line is too long" failure in Windows.
set(CMAKE_NINJA_FORCE_RESPONSE_FILE "ON" CACHE BOOL "Force Ninja to use response files.")

add_executable(Point MACOSX_BUNDLE Point.cpp )
target_link_libraries(Point PRIVATE ${VTK_LIBRARIES} Eigen3::Eigen)

# vtk_module_autoinit is needed
vtk_module_autoinit(
  TARGETS Point
  MODULES ${VTK_LIBRARIES}
)
