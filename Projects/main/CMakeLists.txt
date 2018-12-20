find_package(Eigen3 REQUIRED)
add_library(eigen INTERFACE IMPORTED)

if(CREATE_LOCAL_MAKEFILES)
    file(RELATIVE_PATH CISPBA_RELATIVE_PATH ${CISPBA_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})
    configure_file(${CISPBA_SOURCE_DIR}/Scripts/Makefile.in ${CMAKE_CURRENT_SOURCE_DIR}/Makefile @ONLY)
endif()

file(GLOB MPM_SRC
	"*.h"
	"*.cpp")

message(${MPM_SRC})

add_cispba_executable(main ${MPM_SRC})
target_include_directories(main PUBLIC ${CISPBA_SOURCE_DIR}/tetgen/)
target_include_directories(main PUBLIC ${CISPBA_SOURCE_DIR}/voro++/include)
target_include_directories(main SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
target_compile_options(main PUBLIC -DTETLIBRARY)
target_link_libraries(main partio eigen ${CISPBA_SOURCE_DIR}/tetgen/libtet.a ${CISPBA_SOURCE_DIR}/voro++/libvoro++.a)