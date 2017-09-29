add_executable(laplace_equation laplace_equation.F90)

target_link_libraries(laplace_equation PRIVATE opencmisslibs mpi)

install(TARGETS laplace_equation DESTINATION bin)
