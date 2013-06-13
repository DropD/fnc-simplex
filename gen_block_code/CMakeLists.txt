cmake_minimum_required(VERSION 2.8)

#---------- blocking code ----------
add_custom_target(block ALL DEPENDS gen blockdir)

add_custom_command(
    TARGET block
    COMMAND python gen_block.py
    COMMENT "generate blocking code"
    )

#---------- blockdir ----------

add_custom_target(blockdir)
add_custom_command(
    TARGET blockdir
    COMMAND ${CMAKE_COMMAND} -E make_directory code
    COMMENT "mkdir code"
    )

#---------- gen target ----------

add_custom_target(gen DEPENDS gen_block.py block.py block_avx.py)

#---------- generator script ----------

add_custom_command(
    OUTPUT gen_block.py
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/gen_block.py gen_block.py
    COMMENT "copy lp_block.py"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gen_block.py
    )

#---------- template modules ----------
set(gen_flags --settings="directiveStartToken='%',directiveEndToken='%'" --idir ${CMAKE_BINARY_DIR} --odir ${CMAKE_CURRENT_BINARY_DIR})
set(ch_flags "${gen_flags}")

add_custom_command(
    OUTPUT block.py
    COMMAND cheetah compile ${ch_flags} block.tmpl
    COMMENT "cheetah compile block.tmpl --> block.py"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/block.tmpl
    )

add_custom_command(
    OUTPUT block_avx.py
    COMMAND cheetah compile ${ch_flags} block_avx.tmpl
    COMMENT "cheetah compile block_avx.tmpl --> block_avx.py"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/block_avx.tmpl
    )

#---------- templates ----------

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/block.tmpl
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/block.tmpl block.tmpl
    COMMENT "copy block.tmpl"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/block.tmpl
    )

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/block_avx.tmpl
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/block_avx.tmpl block_avx.tmpl
    COMMENT "copy block_avx.tmpl"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/block_avx.tmpl
    )

set(cleanfiles "${cleanfiles};block.pyc")
set(cleanfiles "${cleanfiles};block_avx.pyc")
set_directory_properties(PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${cleanfiles}")
