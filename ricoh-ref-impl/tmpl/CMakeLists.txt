#---------- templates ----------
add_custom_target(timerC.tmpl DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/timerC.tmpl)
add_custom_command(
    TARGET timerC.tmpl
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/timerC.tmpl timerC.tmpl
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "${CMAKE_CURRENT_SOURCE_DIR}/timerC.tmpl -> ${CMAKE_CURRENT_BINARY_DIR}/timerC.tmpl"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/timerC.tmpl
)

add_custom_target(timerH.tmpl DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/timerH.tmpl)
add_custom_command(
    TARGET timerH.tmpl
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/timerH.tmpl timerH.tmpl
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "${CMAKE_CURRENT_SOURCE_DIR}/timerH.tmpl -> ${CMAKE_CURRENT_BINARY_DIR}/timerH.tmpl"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/timerH.tmpl
)

#---------- clean ----------
set(clean_files "${clean_files};${CMAKE_CURRENT_BINARY_DIR}/timerC.tmpl")
set(clean_files "${clean_files};${CMAKE_CURRENT_BINARY_DIR}/timerH.tmpl")
set_directory_properties(PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${clean_files}")
