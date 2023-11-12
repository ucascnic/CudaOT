# Install script for directory: /home/chenyidong/CudaOT/Common

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
  set(CMAKE_OBJDUMP "/home/spack/spack/opt/spack/linux-debian12-zen2/gcc-11.3.0/binutils-2.40-sjvgy243b7j35dqnfqtlqltlp6lldgef/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/chenyidong/CudaOT/../bin/libCommon.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/chenyidong/CudaOT/../bin" TYPE STATIC_LIBRARY FILES "/home/chenyidong/CudaOT/Common/libCommon.a")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/chenyidong/CudaOT/../include/Common/eps_handler.h;/home/chenyidong/CudaOT/../include/Common/family.h;/home/chenyidong/CudaOT/../include/Common/ErrorCodes.h;/home/chenyidong/CudaOT/../include/Common/GridTools.h;/home/chenyidong/CudaOT/../include/Common/MultiScaleTools.h;/home/chenyidong/CudaOT/../include/Common/PythonTypes.h;/home/chenyidong/CudaOT/../include/Common/TCostFunctionProvider.h;/home/chenyidong/CudaOT/../include/Common/TCouplingHandler.h;/home/chenyidong/CudaOT/../include/Common/TEpsScaling.h;/home/chenyidong/CudaOT/../include/Common/THierarchicalCostFunctionProvider.h;/home/chenyidong/CudaOT/../include/Common/THierarchicalPartition.h;/home/chenyidong/CudaOT/../include/Common/THierarchyBuilder.h;/home/chenyidong/CudaOT/../include/Common/Tools.h;/home/chenyidong/CudaOT/../include/Common/TVarListHandler.h;/home/chenyidong/CudaOT/../include/Common/Verbose.h")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/home/chenyidong/CudaOT/../include/Common" TYPE FILE FILES
    "/home/chenyidong/CudaOT/Common/eps_handler.h"
    "/home/chenyidong/CudaOT/Common/family.h"
    "/home/chenyidong/CudaOT/Common/ErrorCodes.h"
    "/home/chenyidong/CudaOT/Common/GridTools.h"
    "/home/chenyidong/CudaOT/Common/MultiScaleTools.h"
    "/home/chenyidong/CudaOT/Common/PythonTypes.h"
    "/home/chenyidong/CudaOT/Common/TCostFunctionProvider.h"
    "/home/chenyidong/CudaOT/Common/TCouplingHandler.h"
    "/home/chenyidong/CudaOT/Common/TEpsScaling.h"
    "/home/chenyidong/CudaOT/Common/THierarchicalCostFunctionProvider.h"
    "/home/chenyidong/CudaOT/Common/THierarchicalPartition.h"
    "/home/chenyidong/CudaOT/Common/THierarchyBuilder.h"
    "/home/chenyidong/CudaOT/Common/Tools.h"
    "/home/chenyidong/CudaOT/Common/TVarListHandler.h"
    "/home/chenyidong/CudaOT/Common/Verbose.h"
    )
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/chenyidong/CudaOT/Common/Models/cmake_install.cmake")

endif()

