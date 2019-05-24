set(LIBSRC
    "${SCF_SOURCE_DIR}/src/basis.c"
    "${SCF_SOURCE_DIR}/src/THO.c"
    "${SCF_SOURCE_DIR}/src/utils.c"
    )

add_library(intlib STATIC ${LIBSRC})

set(TESTS
    overlap
    )

# set(TEST_SRCS "")
foreach(TEST_NAME ${TESTS})
    # list(APPEND TEST_SRCS test_${TEST_NAME}.c)
    add_executable(test_${TEST_NAME} test_${TEST_NAME}.c)
    target_link_libraries(test_${TEST_NAME} intlib ${EXTERNAL_LIBS})
    add_test(NAME ${TEST_NAME} COMMAND test_${TEST_NAME})
endforeach()
