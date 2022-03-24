# Install script for directory: /home/chenyidong/hpc_ot/Common

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/chenyidong/hpc_ot/../bin/libCommon.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/chenyidong/hpc_ot/../bin" TYPE STATIC_LIBRARY FILES "/home/chenyidong/hpc_ot/Common/libCommon.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/chenyidong/hpc_ot/../include/Common/eps_handler.h;/home/chenyidong/hpc_ot/../include/Common/family.h;/home/chenyidong/hpc_ot/../include/Common/ErrorCodes.h;/home/chenyidong/hpc_ot/../include/Common/GridTools.h;/home/chenyidong/hpc_ot/../include/Common/MultiScaleTools.h;/home/chenyidong/hpc_ot/../include/Common/PythonTypes.h;/home/chenyidong/hpc_ot/../include/Common/TCostFunctionProvider.h;/home/chenyidong/hpc_ot/../include/Common/TCouplingHandler.h;/home/chenyidong/hpc_ot/../include/Common/TEpsScaling.h;/home/chenyidong/hpc_ot/../include/Common/THierarchicalCostFunctionProvider.h;/home/chenyidong/hpc_ot/../include/Common/THierarchicalPartition.h;/home/chenyidong/hpc_ot/../include/Common/THierarchyBuilder.h;/home/chenyidong/hpc_ot/../include/Common/Tools.h;/home/chenyidong/hpc_ot/../include/Common/TVarListHandler.h;/home/chenyidong/hpc_ot/../include/Common/Verbose.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/chenyidong/hpc_ot/../include/Common" TYPE FILE FILES
    "/home/chenyidong/hpc_ot/Common/eps_handler.h"
    "/home/chenyidong/hpc_ot/Common/family.h"
    "/home/chenyidong/hpc_ot/Common/ErrorCodes.h"
    "/home/chenyidong/hpc_ot/Common/GridTools.h"
    "/home/chenyidong/hpc_ot/Common/MultiScaleTools.h"
    "/home/chenyidong/hpc_ot/Common/PythonTypes.h"
    "/home/chenyidong/hpc_ot/Common/TCostFunctionProvider.h"
    "/home/chenyidong/hpc_ot/Common/TCouplingHandler.h"
    "/home/chenyidong/hpc_ot/Common/TEpsScaling.h"
    "/home/chenyidong/hpc_ot/Common/THierarchicalCostFunctionProvider.h"
    "/home/chenyidong/hpc_ot/Common/THierarchicalPartition.h"
    "/home/chenyidong/hpc_ot/Common/THierarchyBuilder.h"
    "/home/chenyidong/hpc_ot/Common/Tools.h"
    "/home/chenyidong/hpc_ot/Common/TVarListHandler.h"
    "/home/chenyidong/hpc_ot/Common/Verbose.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/chenyidong/hpc_ot/Common/Models/cmake_install.cmake")

endif()

