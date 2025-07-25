# Powershell script to set up MSVC CUDA cmake builds that mirror the CI.  Useful for reproducing issues.

param ($EIGEN_CI_ROOTDIR,
       $EIGEN_CI_BUILDDIR,
       $EIGEN_CI_BUILD_TARGET,
       $EIGEN_CI_ADDITIONAL_ARGS,
       $EIGEN_CI_BEFORE_SCRIPT,
       $EIGEN_CI_CMAKE_GENERATOR,
       $EIGEN_CI_MSVC_ARCH,
       $EIGEN_CI_MSVC_VER,
       $EIGEN_CI_TEST_CUSTOM_CXX_FLAGS,

       $EIGEN_CI_CUDA_CXX_FLAGS,
       $EIGEN_CI_CUDA_COMPUTE_ARCH
       )

# Set defaults if not already set.
IF (!$EIGEN_CI_CUDA_CXX_FLAGS)    { $EIGEN_CI_CUDA_CXX_FLAGS    = "" }
IF (!$EIGEN_CI_CUDA_COMPUTE_ARCH) { $EIGEN_CI_CUDA_COMPUTE_ARCH = "50;70" }
IF (!$EIGEN_CI_BUILD_TARGET)      { $EIGEN_CI_BUILD_TARGET      = "buildtests_gpu" }
IF (!$EIGEN_CI_ADDITIONAL_ARGS)   { $EIGEN_CI_ADDITIONAL_ARGS   = '-DCMAKE_CUDA_COMPILER=nvcc.exe -DCMAKE_CUDA_SEPARABLE_COMPILATION=OFF -DEIGEN_TEST_CUDA=on -DEIGEN_CUDA_CXX_FLAGS='+${EIGEN_CI_CUDA_CXX_FLAGS}+' -DEIGEN_CUDA_COMPUTE_ARCH='+${EIGEN_CI_CUDA_COMPUTE_ARCH} }


# Export variables into the global scope
$global:EIGEN_CI_CUDA_CXX_FLAGS    = $EIGEN_CI_CUDA_CXX_FLAGS
$global:EIGEN_CI_CUDA_COMPUTE_ARCH = $EIGEN_CI_CUDA_COMPUTE_ARCH

# Call the generic msvc setup scripts.
function Get-ScriptDirectory { Split-Path $MyInvocation.ScriptName }
$script = Join-Path (Get-ScriptDirectory) 'ci_cmake_msvc.ps1'
& $script $EIGEN_CI_ROOTDIR $EIGEN_CI_BUILDDIR  $EIGEN_CI_BUILD_TARGET $EIGEN_CI_ADDITIONAL_ARGS $EIGEN_CI_BEFORE_SCRIPT $EIGEN_CI_CMAKE_GENERATOR $EIGEN_CI_MSVC_ARCH $EIGEN_CI_MSVC_VER $EIGEN_CI_TEST_CUSTOM_CXX_FLAGS
