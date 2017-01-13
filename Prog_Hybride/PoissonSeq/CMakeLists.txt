set(PROJECT_NAME PoissonSeq)
set(EXEC_NAME ${PROJECT_NAME})

cmake_minimum_required(VERSION 2.8.8)
project(${PROJECT_NAME})

message(STATUS "compiler: ${CMAKE_CXX_COMPILER_ID}")
if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    message(STATUS "Clang detected, adding compile flags")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "GCC detected, adding compile flags")
    add_definitions(-Wall -Wextra)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    message(STATUS "Intel detected, adding compile flags")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    message(STATUS "Visual C/C++ detected, adding compile flags")
    add_definitions(/D_CRT_SECURE_NO_WARNINGS)
endif()

if(DEFINED ENV{POISSON_DIR})
  set(CMAKE_INSTALL_PREFIX $ENV{POISSON_DIR}/install)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "DEBUG")
  #set(CMAKE_BUILD_TYPE "RELEASE")
  #set(CMAKE_BUILD_TYPE "RELWITHDEBINFO")
  #set(CMAKE_BUILD_TYPE "MINSIZEREL")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_COMPILER_IS_GNUCXX)
  add_definitions("-Wall")
  add_definitions("-Wno-comment")
  add_definitions("-std=c++11")
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  add_definitions("-D_CRT_SECURE_NO_WARNINGS")
endif()
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  add_definitions("-std=c++11")
endif()

if (DEFINED ENV{VTK_DIR})
  set(VTK_DIR ENV{VTK_DIR} CACHE PATH "VTK directory override" FORCE)
endif()

find_package(VTK REQUIRED NO_MODULE)
if (${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} VERSION_GREATER 6.0)
  message(STATUS "Use VTK version ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
  include(${VTK_USE_FILE})
else()
  message(STATUS "VTK not found or too old, disable graphics")
  add_definitions("-DNO_GRAPHICS")
endif()

include_directories(
  .
  )

add_executable(
  ${EXEC_NAME}
  main.cpp
  Poisson_Parameters.cpp
  values.cpp
  Poisson.cpp
  timer.cpp
  plot.cpp
  osutils.cpp
  )

target_link_libraries(${EXEC_NAME} ${VTK_LIBRARIES})

install(TARGETS ${EXEC_NAME} DESTINATION .)

