set(ExampleName "laplace_equation")

add_executable(${ExampleName} ${ExampleName}.F90)

target_link_libraries(${ExampleName} PRIVATE opencmisslibs mpi)

install(TARGETS ${ExampleName} DESTINATION bin)

configure_file(./visualise.com ../../../${ExampleName}-build/src/fortran/visualise.com)
