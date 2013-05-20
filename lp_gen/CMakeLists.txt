cmake_minimum_required(VERSION 2.8)

#---------- lp benchmark ----------
add_custom_target(lpbench ALL DEPENDS gen lpbdir)

add_custom_command(
    TARGET lpbench
    COMMAND for i in {1..100}\; do python gen_lp.py 10 lpbench/lpb_1e1_\$\$i\; done
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "generate 100 lps of size 10"
)

add_custom_command(
    TARGET lpbench
    COMMAND for i in {1..100}\; do python gen_lp.py 20 lpbench/lpb_2e1_\$\$i\; done
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "generate 100 lps of size 20"
)

add_custom_command(
    TARGET lpbench
    COMMAND for i in {1..100}\; do python gen_lp.py 50 lpbench/lpb_5e1_\$\$i\; done
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "generate 100 lps of size 50"
)

add_custom_command(
    TARGET lpbench
    COMMAND for i in {1..100}\; do python gen_lp.py 100 lpbench/lpb_1e2_\$\$i\; done
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "generate 100 lps of size 100"
)

add_custom_command(
    TARGET lpbench
    COMMAND for i in {1..10}\; do python gen_lp.py 1000 lpbench/lpb_1e3_\$\$i\; done
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "generate 10 lps of size 1000"
)
#---------- benchmark dir ----------

add_custom_target(lpbdir)
add_custom_command(
    TARGET lpbdir
    COMMAND ${CMAKE_COMMAND} -E make_directory lpbench
    COMMENT "creating lpbench dir"
)

#---------- gen target ----------
add_custom_target(gen DEPENDS gen_lp.py donjlp.py cplexlp.py)

#---------- generator scripts ----------
add_custom_command(
    OUTPUT gen_lp.py
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/gen_lp.py gen_lp.py
    COMMENT "copy lp_gen.py"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/gen_lp.py
)
#---------- template modules ----------
add_custom_command(
    OUTPUT cplexlp.py
    COMMAND cheetah compile cplexlp.tmpl
    COMMENT "cheetah compile cplexlp.tmpl --> cplexlp.py"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cplexlp.tmpl
)

add_custom_command(
    OUTPUT donjlp.py
    COMMAND cheetah compile donjlp.tmpl
    COMMENT "cheetah compile donjlp.tmpl --> donjlp.py"
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/donjlp.tmpl
)

#---------- templates ----------
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cplexlp.tmpl
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/cplexlp.tmpl cplexlp.tmpl
    COMMENT "copy cplexlp.tmpl"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/cplexlp.tmpl
)

add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/donjlp.tmpl
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/donjlp.tmpl donjlp.tmpl
    COMMENT "copy donjlp.tmpl"
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/donjlp.tmpl
)

set(cleanfiles "${cleanfiles};cplexlp.pyc")
set(cleanfiles "${cleanfiles};donjlp.pyc")
set_directory_properties(PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "${cleanfiles}")
