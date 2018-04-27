set(ExampleName "laplace_equation")

cmake_minimum_required(VERSION 3.4 FATAL_ERROR)

project(${ExampleName} VERSION 1.0 LANGUAGES C Fortran)

find_package(OpenCMISSLibs 1.3.0 REQUIRED COMPONENTS Iron CONFIG)

add_subdirectory(src/fortran)
add_subdirectory(src/c)
add_subdirectory(src/number)
