# Install script for directory: C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "C:/Program Files (x86)/SeismicSimulation")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "C:/mingw64/bin/objdump.exe")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/AdolcForward"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/AlignedVector3"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/ArpackSupport"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/AutoDiff"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/BVH"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/EulerAngles"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/FFT"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/IterativeSolvers"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/KroneckerProduct"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/LevenbergMarquardt"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/MatrixFunctions"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/MPRealSupport"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/NNLS"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/NonLinearOptimization"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/NumericalDiff"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/OpenGLSupport"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/Polynomials"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/SparseExtra"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/SpecialFunctions"
    "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/Splines"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Devel" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/extern/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/buildE/extern/eigen/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "C:/Users/izisr/OneDrive/Documentos/FACULDADE/2025-1/ProgramacaoPratica/Simulador/src/buildE/extern/eigen/unsupported/Eigen/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
