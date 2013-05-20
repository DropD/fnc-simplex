#---------- generator scripts ----------
add_custom_target(gen-timer.py DEPENDS timerC.py timerH.py)
add_custom_command(
    TARGET gen-timer.py
    PRE_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/gen-timer.py gen-timer.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "${CMAKE_CURRENT_SOURCE_DIR}/gen-timer.py -> ${CMAKE_CURRENT_BINARY_DIR}/gen-timer.py"
    DEPENDS timerC.py timerH.py
)

#---------- template modules ----------
set(ch_dir_f --idir=${CMAKE_BINARY_DIR}/tmpl)
set(ch_setings --settings="directiveStartToken='%',directiveEndToken='%'")
set(ch_flags ${ch_dir_f} ${ch_setings})

add_custom_target(timerC.py DEPENDS timerC.tmpl)
add_custom_command(
    TARGET timerC.py
    PRE_BUILD
    COMMAND cheetah compile ${ch_flags} timerC.tmpl
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "${CMAKE_BINARY_DIR}/tmpl/timerC.tmpl -> ${CMAKE_CURRENT_BINARY_DIR}/timerC.py"
    DEPENDS timerC_tmpl
)

add_custom_target(timerH.py DEPENDS timerH.tmpl)
add_custom_command(
    TARGET timerH.py
    PRE_BUILD
    COMMAND cheetah compile ${ch_flags} timerH.tmpl
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "${CMAKE_BINARY_DIR}/tmpl/timerH.tmpl -> ${CMAKE_CURRENT_BINARY_DIR}/timerH.py"
    DEPENDS timerH_tmpl
)

#---------- clean ----------
set(clean_files "${clean_files};${CMAKE_CURRENT_BINARY_DIR}/timerC.py")
set(clean_files "${clean_files};${CMAKE_CURRENT_BINARY_DIR}/timerH.py")
set(clean_files "${clean_files};${CMAKE_CURRENT_BINARY_DIR}/gen-timer.py")
set_directory_properties(PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${clean_files}")
