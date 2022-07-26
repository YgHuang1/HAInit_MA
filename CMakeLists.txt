project(HAInit_MA)
cmake_minimum_required(VERSION 3.5)
find_package(OpenMP)
set( CMAKE_BUILD_TYPE "Release" )

set( CMAKE_CXX_FLAGS  " ${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++14" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -DNDEBUG -DIL_STD -pthread -m64 ${LINUX_CXX_FLAGS}")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${HAInit_MA_SOURCE_DIR}/bin)

# set(SRC_list ./textFunc.cpp)

include_directories(HAInit)
find_package(OpenMVG REQUIRED)
include_directories( ${PCL_INCLUDE_DIRS})

aux_source_directory(HAInit SRC_list1)

if(NOT DEFINED CMAKE_SUPPRESS_DEVELOPER_WARNINGS)
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE INTERNAL "No dev warnings")
endif()
find_package(PCL 1.11 REQUIRED)
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

find_package(Boost REQUIRED)
message("boost find: ${Boost_FOUND}")
include_directories( ${Boost_INCLUDE_DIRS} )
link_directories( ${Boost_LIBRARY_DIRS} )

add_executable(HAPipeline HAPipeline.cpp ${SRC_list1})
target_link_libraries(HAPipeline
PRIVATE
    OpenMVG::openMVG_stlplus
    ${PCL_LIBRARIES}
    ${STLPLUS_LIBRARY}
    boost_thread)

