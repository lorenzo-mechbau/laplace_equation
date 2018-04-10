add_executable(${ExampleName}_C ${ExampleName}.c)

target_link_libraries(${ExampleName}_C PRIVATE opencmisslibs)

install(TARGETS ${ExampleName}_C DESTINATION bin)

configure_file(visualise.cmgui visualise.cmgui)
