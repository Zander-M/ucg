# Install script for directory: C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Users/sherl/Source/Repos/Zander-M/ucg/ext/eigen/out/install/x64-Debug")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/AdolcForward"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/AlignedVector3"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/ArpackSupport"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/AutoDiff"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/BVH"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/EulerAngles"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/FFT"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/IterativeSolvers"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/KroneckerProduct"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/LevenbergMarquardt"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/MatrixFunctions"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/MoreVectorization"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/MPRealSupport"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/NonLinearOptimization"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/NumericalDiff"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/OpenGLSupport"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/Polynomials"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/Skyline"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/SparseExtra"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/SpecialFunctions"
    "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/sherl/source/repos/Zander-M/ucg/ext/eigen/out/build/x64-Debug/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

