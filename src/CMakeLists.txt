cmake_minimum_required(VERSION 3.6)

set(CMAKE_CONFIGURATION_TYPES Debug Release)

project(arap)

# Options
set(LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../libs CACHE PATH "Path to lib folder")
set(Eigen3_INCLUDE_DIR ${LIBRARY_DIR}/Eigen/ CACHE PATH "Path to Eigen source folder")

set(IGL_INCLUDE_DIR ${LIBRARY_DIR}/libigl-master/include/ CACHE PATH "Path to igl include folder")

set(GLFW_INCLUDE_DIR ${LIBRARY_DIR}/GLFW/include CACHE PATH "Path to GLFW include folder")
set(GLFW_LIBRARY_DIR ${GLFW_INCLUDE_DIR}/../lib-vc2019 CACHE PATH "Path to GLFW lib folder")

set(GLEW_INCLUDE_DIR ${LIBRARY_DIR}/GLEW/include CACHE PATH "Path to GLEW include folder")
set(GLEW_LIBRARY_DIR ${GLEW_INCLUDE_DIR}/../lib/Release/x64 CACHE PATH "Path to GLEW lib folder")

set(glad_INCLUDE_DIR ${IGL_INCLUDE_DIR}/../external/glad/include CACHE PATH "Path to glad include folder")
set(glog_INCLUDE_DIR ${LIBRARY_DIR}/glog-lib/include CACHE PATH "Path to glog source folder")
set(Ceres_DIR ${LIBRARY_DIR}/Ceres/CMake CACHE PATH "Path to Ceres build folder")
set(Ceres_INCLUDE_DIR ${Ceres_DIR}/../include)
set(Flann_INCLUDE_DIR ${LIBRARY_DIR}/Flann/ CACHE PATH "Path to Flann source folder")

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} ${Ceres_DIR})

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

# Use pre-compiled libraries on Windows
if(WIN32)
    set(FreeImage_INCLUDE_DIR ${LIBRARY_DIR}/FreeImage/x64/ CACHE PATH "Path to FreeImage header file")
    set(FreeImage_LIBRARY_DIR ${LIBRARY_DIR}/FreeImage/x64/ CACHE PATH "Path to FreeImage .lib/.dll folder")
endif(WIN32)

# Set C++ flags
set(CMAKE_CXX_STANDARD 14)

find_package(Ceres REQUIRED)
find_package(glog REQUIRED ${glog_DIR_HINT})

add_definitions("-D_DISABLE_EXTENDED_ALIGNED_STORAGE")

if (MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /bigobj")
endif()

message(STATUS "ceres library name: ${Ceres_LIBRARIES}")
message(STATUS "path to ceres cmake information: ${Ceres_DIR}")
message(STATUS "path to ceres include directory: ${Ceres_INCLUDE_DIR}")
message(STATUS "glog library name: ${glog_LIBRARIES}")
message(STATUS "path to glog include directory = ${glog_INCLUDE_DIR}")


# Set files to be compiled
set(HEADER_FILES 
    Eigen.h 
    SimpleMesh.h 
    FreeImageHelper.h
    GUI.h
    ArapDeformer.h
)
set(SOURCE_FILES 
    FreeImageHelper.cpp
    glad.c
    ArapDeformer.cpp
)

link_directories(${FreeImage_LIBRARY_DIR} ${GLFW_LIBRARY_DIR} ${GLEW_LIBRARY_DIR})
add_executable(arap main.cpp ${HEADER_FILES} ${SOURCE_FILES})
target_link_libraries(arap ceres freeimage opengl32 glew32s glfw3)
target_include_directories(arap PUBLIC ${Eigen3_INCLUDE_DIR} ${IGL_INCLUDE_DIR} ${GLFW_INCLUDE_DIR} ${GLEW_INCLUDE_DIR} ${glad_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} ${Ceres_INCLUDE_DIR} ${glog_INCLUDE_DIR} ${FreeImage_INCLUDE_DIR} ${Flann_INCLUDE_DIR})

# # Visual Studio properties
# set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT arap)
# set_property(TARGET arap PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/)

# On Windows copy dll to executable folder
if(WIN32)
    add_custom_command(TARGET arap POST_BUILD        # Adds a post-build event to MyTest
        COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
            "${LIBRARY_DIR}/FreeImage/x64/FreeImage.dll"      # <--this is in-file
            $<TARGET_FILE_DIR:arap>)                 # <--this is out-file path
endif(WIN32)
