add_executable(${ExampleName}_Fortran_Number ${ExampleName}_number.F90)

target_link_libraries(${ExampleName}_Fortran_Number PRIVATE opencmisslibs)

install(TARGETS ${ExampleName}_Fortran_Number DESTINATION bin)

configure_file(visualise.cmgui visualise.cmgui)
