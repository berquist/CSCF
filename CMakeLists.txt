cmake_minimum_required(VERSION 2.8.9)
project(SCF)

if(MSVC)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4 /WX")
else()
    # Use this to treat warnings as errors
    # set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Werror")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
endif()

include_directories(include)

file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")
add_executable(scf ${SOURCES})
set(EXTERNAL_LIBS m cblas lapack)
target_link_libraries(scf ${EXTERNAL_LIBS})

enable_testing()
add_subdirectory(tests)

option(BUILD_DOC "Build documentation" ON)

find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
