
macro(ei_add_property prop value)
  get_property(previous GLOBAL PROPERTY ${prop})
  if ((NOT previous) OR (previous STREQUAL ""))
    set_property(GLOBAL PROPERTY ${prop} "${value}")
  else()
    set_property(GLOBAL PROPERTY ${prop} "${previous} ${value}")
  endif()
endmacro()

#internal. See documentation of ei_add_test for details.
macro(ei_add_test_internal testname testname_with_suffix)
  set(targetname ${testname_with_suffix})

  if(EIGEN_ADD_TEST_FILENAME_EXTENSION)
    set(filename ${testname}.${EIGEN_ADD_TEST_FILENAME_EXTENSION})
  else()
    set(filename ${testname}.cpp)
  endif()

  # Add the current target to the list of subtest targets
  get_property(EIGEN_SUBTESTS_LIST GLOBAL PROPERTY EIGEN_SUBTESTS_LIST)
  set(EIGEN_SUBTESTS_LIST "${EIGEN_SUBTESTS_LIST}${targetname}\n")
  set_property(GLOBAL PROPERTY EIGEN_SUBTESTS_LIST "${EIGEN_SUBTESTS_LIST}")

  set(is_gpu_test OFF)
  if(EIGEN_ADD_TEST_FILENAME_EXTENSION STREQUAL cu)
    set(is_gpu_test ON)
    if(EIGEN_TEST_HIP)
      hip_reset_flags()
      hip_add_executable(${targetname} ${filename} HIPCC_OPTIONS -std=c++14)
      target_compile_definitions(${targetname} PRIVATE -DEIGEN_USE_HIP)
      set_property(TARGET ${targetname} PROPERTY HIP_ARCHITECTURES gfx900 gfx906 gfx908 gfx90a gfx940 gfx941 gfx942 gfx1030)
    elseif(EIGEN_TEST_CUDA_CLANG)
      set_source_files_properties(${filename} PROPERTIES LANGUAGE CXX)
      
      if(CUDA_64_BIT_DEVICE_CODE AND (EXISTS "${CUDA_TOOLKIT_ROOT_DIR}/lib64"))
        link_directories("${CUDA_TOOLKIT_ROOT_DIR}/lib64")
      else()
        link_directories("${CUDA_TOOLKIT_ROOT_DIR}/lib")
      endif()

      add_executable(${targetname} ${filename})
      set(CUDA_CLANG_LINK_LIBRARIES "cudart_static" "cuda" "dl" "pthread")
      if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set(CUDA_CLANG_LINK_LIBRARIES ${CUDA_CLANG_LINK_LIBRARIES} "rt")
      endif()
      target_link_libraries(${targetname} ${CUDA_CLANG_LINK_LIBRARIES})
    else()
      cuda_add_executable(${targetname} ${filename})
    endif()
  else()
    add_executable(${targetname} ${filename})
  endif()

  add_dependencies(buildtests ${targetname})
  
  if (is_gpu_test)
    add_dependencies(buildtests_gpu ${targetname})
  endif()

  if(EIGEN_NO_ASSERTION_CHECKING)
    target_compile_definitions(${targetname} PRIVATE EIGEN_NO_ASSERTION_CHECKING=1)
  else()
    if(EIGEN_DEBUG_ASSERTS)
      target_compile_definitions(${targetname} PRIVATE EIGEN_DEBUG_ASSERTS=1)
    endif()
  endif()

  target_compile_definitions(${targetname} PRIVATE EIGEN_TEST_MAX_SIZE=${EIGEN_TEST_MAX_SIZE})

  if(MSVC)
    target_compile_options(${targetname} PRIVATE "/bigobj")
  endif()

  # let the user pass flags.
  if(${ARGC} GREATER 2)
    separate_arguments(compile_options NATIVE_COMMAND ${ARGV2})
    target_compile_options(${targetname} PRIVATE ${compile_options})
  endif()

  if(EIGEN_TEST_CUSTOM_CXX_FLAGS)
    target_compile_options(${targetname} PRIVATE ${EIGEN_TEST_CUSTOM_CXX_FLAGS})
  endif()

  if(EIGEN_STANDARD_LIBRARIES_TO_LINK_TO)
    target_link_libraries(${targetname} ${EIGEN_STANDARD_LIBRARIES_TO_LINK_TO})
  endif()
  if(EXTERNAL_LIBS)
    target_link_libraries(${targetname} ${EXTERNAL_LIBS})
  endif()
  if(EIGEN_TEST_CUSTOM_LINKER_FLAGS)
    target_link_libraries(${targetname} ${EIGEN_TEST_CUSTOM_LINKER_FLAGS})
  endif()
  target_link_libraries(${targetname} Eigen3::Eigen)

  if(${ARGC} GREATER 3)
    set(libs_to_link ${ARGV3})
    # it could be that some cmake module provides a bad library string " "  (just spaces),
    # and that severely breaks target_link_libraries ("can't link to -l-lstdc++" errors).
    # so we check for strings containing only spaces.
    string(STRIP "${libs_to_link}" libs_to_link_stripped)
    string(LENGTH "${libs_to_link_stripped}" libs_to_link_stripped_length)
    if(${libs_to_link_stripped_length} GREATER 0)
      # notice: no double quotes around ${libs_to_link} here. It may be a list.
      target_link_libraries(${targetname} ${libs_to_link})
    endif()
  endif()

  add_test(${testname_with_suffix} "${targetname}")

  # Specify target and test labels according to EIGEN_CURRENT_SUBPROJECT
  get_property(current_subproject GLOBAL PROPERTY EIGEN_CURRENT_SUBPROJECT)
  if ((current_subproject) AND (NOT (current_subproject STREQUAL "")))
    set_property(TARGET ${targetname} PROPERTY LABELS "Build${current_subproject}")
    add_dependencies("Build${current_subproject}" ${targetname})
    set_property(TEST ${testname_with_suffix} PROPERTY LABELS "${current_subproject}")
  endif()
  if (is_gpu_test)
    # Add gpu tag for testing only GPU tests.
    set_property(TEST ${testname_with_suffix} APPEND PROPERTY LABELS "gpu")
  endif()
  
  if(EIGEN_SYCL)
    # Force include of the SYCL file at the end to avoid errors.
    set_property(TARGET ${targetname} PROPERTY COMPUTECPP_INCLUDE_AFTER 1)
    # Link against pthread and add sycl to target
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    find_package(Threads REQUIRED)
    target_link_libraries(${targetname} Threads::Threads)
    add_sycl_to_target(TARGET ${targetname} SOURCES ${filename})
  endif(EIGEN_SYCL)
endmacro(ei_add_test_internal)
# Macro to add a test
#
# the unique mandatory parameter testname must correspond to a file
# <testname>.cpp which follows this pattern:
#
# #include "main.h"
# void test_<testname>() { ... }
#
# Depending on the contents of that file, this macro can have 2 behaviors,
# see below.
#
# The optional 2nd parameter is libraries to link to.
#
# A. Default behavior
#
# this macro adds an executable <testname> as well as a ctest test
# named <testname> too.
#
# On platforms with bash simply run:
#   "ctest -V" or "ctest -V -R <testname>"
# On other platform use ctest as usual
#
# B. Multi-part behavior
#
# If the source file matches the regexp
#    CALL_SUBTEST_[0-9]+|EIGEN_TEST_PART_[0-9]+
# then it is interpreted as a multi-part test. The behavior then depends on the
# CMake option EIGEN_SPLIT_LARGE_TESTS, which is ON by default.
#
# If EIGEN_SPLIT_LARGE_TESTS is OFF, the behavior is the same as in A (the multi-part
# aspect is ignored).
#
# If EIGEN_SPLIT_LARGE_TESTS is ON, the test is split into multiple executables
#   test_<testname>_<N>
# where N runs from 1 to the greatest occurrence found in the source file. Each of these
# executables is built passing -DEIGEN_TEST_PART_N. This allows to split large tests
# into smaller executables.
#
# Moreover, targets <testname> are still generated, they
# have the effect of building all the parts of the test.
#
# Again, ctest -R allows to run all matching tests.
macro(ei_add_test testname)
  get_property(EIGEN_TESTS_LIST GLOBAL PROPERTY EIGEN_TESTS_LIST)
  set(EIGEN_TESTS_LIST "${EIGEN_TESTS_LIST}${testname}\n")
  set_property(GLOBAL PROPERTY EIGEN_TESTS_LIST "${EIGEN_TESTS_LIST}")

  if(EIGEN_ADD_TEST_FILENAME_EXTENSION)
    set(filename ${testname}.${EIGEN_ADD_TEST_FILENAME_EXTENSION})
  else()
    set(filename ${testname}.cpp)
  endif()

  file(READ "${filename}" test_source)
  string(REGEX MATCHALL "CALL_SUBTEST_[0-9]+|EIGEN_TEST_PART_[0-9]+|EIGEN_SUFFIXES(;[0-9]+)+"
         occurrences "${test_source}")
  string(REGEX REPLACE "CALL_SUBTEST_|EIGEN_TEST_PART_|EIGEN_SUFFIXES" "" suffixes "${occurrences}")
  list(REMOVE_DUPLICATES suffixes)
  set(explicit_suffixes "")
  if( (NOT EIGEN_SPLIT_LARGE_TESTS) AND suffixes)
    # Check whether we have EIGEN_TEST_PART_* statements, in which case we likely must enforce splitting.
    # For instance, indexed_view activate a different c++ version for each part.
    string(REGEX MATCHALL "EIGEN_TEST_PART_[0-9]+" occurrences "${test_source}")
    string(REGEX REPLACE "EIGEN_TEST_PART_" "" explicit_suffixes "${occurrences}")
    list(REMOVE_DUPLICATES explicit_suffixes)
  endif()
  if( (EIGEN_SPLIT_LARGE_TESTS AND suffixes) OR explicit_suffixes)
    add_custom_target(${testname})
    foreach(suffix ${suffixes})
      ei_add_test_internal(${testname} ${testname}_${suffix} "${ARGV1}" "${ARGV2}")
      add_dependencies(${testname} ${testname}_${suffix})
      target_compile_definitions(${testname}_${suffix} PRIVATE -DEIGEN_TEST_PART_${suffix}=1)
    endforeach()
  else()
    ei_add_test_internal(${testname} ${testname} "${ARGV1}" "${ARGV2}")
    target_compile_definitions(${testname} PRIVATE -DEIGEN_TEST_PART_ALL=1)
  endif()
endmacro()

# adds a failtest, i.e. a test that succeed if the program fails to compile
# note that the test runner for these is CMake itself, when passed -DEIGEN_FAILTEST=ON
# so here we're just running CMake commands immediately, we're not adding any targets.
macro(ei_add_failtest testname)

  set(test_target_ok ${testname}_ok)
  set(test_target_ko ${testname}_ko)

  # Add executables
  add_executable(${test_target_ok} ${testname}.cpp)
  add_executable(${test_target_ko} ${testname}.cpp)

  # Remove them from the normal build process
  set_target_properties(${test_target_ok} ${test_target_ko} PROPERTIES
                        EXCLUDE_FROM_ALL TRUE
                        EXCLUDE_FROM_DEFAULT_BUILD TRUE)

  # Configure the failing test
  target_compile_definitions(${test_target_ko} PRIVATE EIGEN_SHOULD_FAIL_TO_BUILD)

  # Add the tests to ctest.
  add_test(NAME ${test_target_ok}
          COMMAND ${CMAKE_COMMAND} --build . --target ${test_target_ok} --config $<CONFIG>
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  add_test(NAME ${test_target_ko}
          COMMAND ${CMAKE_COMMAND} --build . --target ${test_target_ko} --config $<CONFIG>
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # Expect the second test to fail
  set_tests_properties(${test_target_ko} PROPERTIES WILL_FAIL TRUE)
endmacro()

# print a summary of the different options
macro(ei_testing_print_summary)
  message(STATUS "************************************************************")
  message(STATUS "***    Eigen's unit tests configuration summary          ***")
  message(STATUS "************************************************************")
  message(STATUS "")
  message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
  message(STATUS "Build site:        ${SITE}")
  message(STATUS "Build string:      ${BUILDNAME}")
  get_property(EIGEN_TESTING_SUMMARY GLOBAL PROPERTY EIGEN_TESTING_SUMMARY)
  get_property(EIGEN_TESTED_BACKENDS GLOBAL PROPERTY EIGEN_TESTED_BACKENDS)
  get_property(EIGEN_MISSING_BACKENDS GLOBAL PROPERTY EIGEN_MISSING_BACKENDS)
  message(STATUS "Enabled backends:  ${EIGEN_TESTED_BACKENDS}")
  message(STATUS "Disabled backends: ${EIGEN_MISSING_BACKENDS}")

  if(EIGEN_DEFAULT_TO_ROW_MAJOR)
    message(STATUS "Default order:     Row-major")
  else()
    message(STATUS "Default order:     Column-major")
  endif()

  if(EIGEN_TEST_NO_EXPLICIT_ALIGNMENT)
    message(STATUS "Explicit alignment (hence vectorization) disabled")
  elseif(EIGEN_TEST_NO_EXPLICIT_VECTORIZATION)
    message(STATUS "Explicit vectorization disabled (alignment kept enabled)")
  else()

  message(STATUS "Maximal matrix/vector size: ${EIGEN_TEST_MAX_SIZE}")

    if(EIGEN_TEST_SSE2)
      message(STATUS "SSE2:              ON")
    else()
      message(STATUS "SSE2:              Using architecture defaults")
    endif()

    if(EIGEN_TEST_SSE3)
      message(STATUS "SSE3:              ON")
    else()
      message(STATUS "SSE3:              Using architecture defaults")
    endif()

    if(EIGEN_TEST_SSSE3)
      message(STATUS "SSSE3:             ON")
    else()
      message(STATUS "SSSE3:             Using architecture defaults")
    endif()

    if(EIGEN_TEST_SSE4_1)
      message(STATUS "SSE4.1:            ON")
    else()
      message(STATUS "SSE4.1:            Using architecture defaults")
    endif()

    if(EIGEN_TEST_SSE4_2)
      message(STATUS "SSE4.2:            ON")
    else()
      message(STATUS "SSE4.2:            Using architecture defaults")
    endif()

    if(EIGEN_TEST_AVX)
      message(STATUS "AVX:               ON")
    else()
      message(STATUS "AVX:               Using architecture defaults")
    endif()

    if(EIGEN_TEST_AVX2)
      message(STATUS "AVX2:              ON")
    else()
      message(STATUS "AVX2:              Using architecture defaults")
    endif()

    if(EIGEN_TEST_FMA)
      message(STATUS "FMA:               ON")
    else()
      message(STATUS "FMA:               Using architecture defaults")
    endif()

    if(EIGEN_TEST_AVX512)
      message(STATUS "AVX512:            ON")
    else()
      message(STATUS "AVX512:            Using architecture defaults")
    endif()

    if(EIGEN_TEST_AVX512DQ)
      message(STATUS "AVX512DQ:          ON")
    else()
      message(STATUS "AVX512DQ:          Using architecture defaults")
    endif()

    if(EIGEN_TEST_ALTIVEC)
      message(STATUS "Altivec:           ON")
    else()
      message(STATUS "Altivec:           Using architecture defaults")
    endif()

    if(EIGEN_TEST_VSX)
      message(STATUS "VSX:               ON")
    else()
      message(STATUS "VSX:               Using architecture defaults")
    endif()

    if(EIGEN_TEST_MSA)
      message(STATUS "MIPS MSA:          ON")
    else()
      message(STATUS "MIPS MSA:          Using architecture defaults")
    endif()

    if(EIGEN_TEST_NEON)
      message(STATUS "ARM NEON:          ON")
    else()
      message(STATUS "ARM NEON:          Using architecture defaults")
    endif()

    if(EIGEN_TEST_NEON64)
      message(STATUS "ARMv8 NEON:        ON")
    else()
      message(STATUS "ARMv8 NEON:        Using architecture defaults")
    endif()

    if(EIGEN_TEST_ZVECTOR)
      message(STATUS "S390X ZVECTOR:     ON")
    else()
      message(STATUS "S390X ZVECTOR:     Using architecture defaults")
    endif()

    if(EIGEN_TEST_SYCL)
      if(EIGEN_SYCL_TRISYCL)
        message(STATUS "SYCL:              ON (using triSYCL)")
      elseif(EIGEN_SYCL_ComputeCpp)
        message(STATUS "SYCL:              ON (using computeCPP)")
      elseif(EIGEN_SYCL_DPCPP)
        message(STATUS "SYCL:              ON (using DPCPP)")
      endif()
    else()
      message(STATUS "SYCL:              OFF")
    endif()
    if(EIGEN_TEST_CUDA)
      if(EIGEN_TEST_CUDA_CLANG)
        message(STATUS "CUDA:              ON (using clang)")
      else()
        message(STATUS "CUDA:              ON (using nvcc)")
      endif()
    else()
      message(STATUS "CUDA:              OFF")
    endif()
    if(EIGEN_TEST_HIP)
      message(STATUS "HIP:               ON (using hipcc)")
    else()
      message(STATUS "HIP:               OFF")
    endif()

  endif() # vectorization / alignment options

  message(STATUS "\n${EIGEN_TESTING_SUMMARY}")

  message(STATUS "************************************************************")
endmacro()

macro(ei_init_testing)
  define_property(GLOBAL PROPERTY EIGEN_CURRENT_SUBPROJECT BRIEF_DOCS " " FULL_DOCS " ")
  define_property(GLOBAL PROPERTY EIGEN_TESTED_BACKENDS BRIEF_DOCS " " FULL_DOCS " ")
  define_property(GLOBAL PROPERTY EIGEN_MISSING_BACKENDS BRIEF_DOCS " " FULL_DOCS " ")
  define_property(GLOBAL PROPERTY EIGEN_TESTING_SUMMARY BRIEF_DOCS " " FULL_DOCS " ")
  define_property(GLOBAL PROPERTY EIGEN_TESTS_LIST BRIEF_DOCS " " FULL_DOCS " ")
  define_property(GLOBAL PROPERTY EIGEN_SUBTESTS_LIST BRIEF_DOCS " " FULL_DOCS " ")

  set_property(GLOBAL PROPERTY EIGEN_TESTED_BACKENDS "")
  set_property(GLOBAL PROPERTY EIGEN_MISSING_BACKENDS "")
  set_property(GLOBAL PROPERTY EIGEN_TESTING_SUMMARY "")
  set_property(GLOBAL PROPERTY EIGEN_TESTS_LIST "")
  set_property(GLOBAL PROPERTY EIGEN_SUBTESTS_LIST "")

  define_property(GLOBAL PROPERTY EIGEN_FAILTEST_FAILURE_COUNT BRIEF_DOCS " " FULL_DOCS " ")
  define_property(GLOBAL PROPERTY EIGEN_FAILTEST_COUNT BRIEF_DOCS " " FULL_DOCS " ")

  set_property(GLOBAL PROPERTY EIGEN_FAILTEST_FAILURE_COUNT "0")
  set_property(GLOBAL PROPERTY EIGEN_FAILTEST_COUNT "0")

  # uncomment anytime you change the ei_get_compilerver_from_cxx_version_string macro
  # ei_test_get_compilerver_from_cxx_version_string()
endmacro()

macro(ei_set_sitename)
  # if the sitename is not yet set, try to set it
  if(NOT ${SITE} OR ${SITE} STREQUAL "")
    set(eigen_computername $ENV{COMPUTERNAME})
    set(eigen_hostname $ENV{HOSTNAME})
    if(eigen_hostname)
      set(SITE ${eigen_hostname})
    elseif(eigen_computername)
      set(SITE ${eigen_computername})
    endif()
  endif()
  # in case it is already set, enforce lower case
  if(SITE)
    string(TOLOWER ${SITE} SITE)
  endif()
endmacro()

macro(ei_get_compilerver VAR)
    if(MSVC)
      set(${VAR} "${CMAKE_CXX_COMPILER_VERSION}")
    elseif(${CMAKE_CXX_COMPILER_ID} MATCHES "PGI")
      set(${VAR} "${CMAKE_CXX_COMPILER_ID}-${CMAKE_CXX_COMPILER_VERSION}")
    else()
    # on all other system we rely on ${CMAKE_CXX_COMPILER}
    # supporting a "--version" or "/version" flag

    if(WIN32 AND ${CMAKE_CXX_COMPILER_ID} EQUAL "Intel")
      set(EIGEN_CXX_FLAG_VERSION "/version")
    else()
      set(EIGEN_CXX_FLAG_VERSION "--version")
    endif()

    execute_process(COMMAND ${CMAKE_CXX_COMPILER} ${EIGEN_CXX_FLAG_VERSION}
                    OUTPUT_VARIABLE eigen_cxx_compiler_version_string OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(REGEX REPLACE "^[ \n\r]+" "" eigen_cxx_compiler_version_string ${eigen_cxx_compiler_version_string})
    string(REGEX REPLACE "[\n\r].*"  ""  eigen_cxx_compiler_version_string  ${eigen_cxx_compiler_version_string})

    ei_get_compilerver_from_cxx_version_string("${eigen_cxx_compiler_version_string}" CNAME CVER)
    set(${VAR} "${CNAME}-${CVER}")

  endif()
endmacro()

# Extract compiler name and version from a raw version string
# WARNING: if you edit this macro, then please test it by uncommenting
# the testing macro call in ei_init_testing() of the EigenTesting.cmake file.
# See also the ei_test_get_compilerver_from_cxx_version_string macro at the end
# of the file
macro(ei_get_compilerver_from_cxx_version_string VERSTRING CNAME CVER)
  # extract possible compiler names
  string(REGEX MATCH "g\\+\\+"      ei_has_gpp    ${VERSTRING})
  string(REGEX MATCH "llvm|LLVM"    ei_has_llvm   ${VERSTRING})
  string(REGEX MATCH "gcc|GCC"      ei_has_gcc    ${VERSTRING})
  string(REGEX MATCH "icpc|ICC"     ei_has_icpc   ${VERSTRING})
  string(REGEX MATCH "clang|CLANG"  ei_has_clang  ${VERSTRING})
  string(REGEX MATCH "mingw32"      ei_has_mingw  ${VERSTRING})

  # combine them
  if((ei_has_llvm) AND (ei_has_gpp OR ei_has_gcc))
    set(${CNAME} "llvm-g++")
  elseif((ei_has_llvm) AND (ei_has_clang))
    set(${CNAME} "llvm-clang++")
  elseif(ei_has_clang)
    set(${CNAME} "clang++")
  elseif ((ei_has_mingw) AND (ei_has_gpp OR ei_has_gcc))
    set(${CNAME} "mingw32-g++")
  elseif(ei_has_icpc)
    set(${CNAME} "icpc")
  elseif(ei_has_gpp OR ei_has_gcc)
    set(${CNAME} "g++")
  else()
    set(${CNAME} "_")
  endif()

  # extract possible version numbers
  # first try to extract 3 isolated numbers:
  string(REGEX MATCH " [0-9]+\\.[0-9]+\\.[0-9]+" eicver ${VERSTRING})
  if(NOT eicver)
    # try to extract 2 isolated ones:
    string(REGEX MATCH " [0-9]+\\.[0-9]+" eicver ${VERSTRING})
    if(NOT eicver)
      # try to extract 3:
      string(REGEX MATCH "[^0-9][0-9]+\\.[0-9]+\\.[0-9]+" eicver ${VERSTRING})
      if(NOT eicver)
        # try to extract 2:
        string(REGEX MATCH "[^0-9][0-9]+\\.[0-9]+" eicver ${VERSTRING})
        if (NOT eicver AND ei_has_mingw)
          # try to extract 1 number plus suffix:
          string(REGEX MATCH "[^0-9][0-9]+-win32" eicver ${VERSTRING})          
        endif()
      endif()
    endif()
  endif()
  
  if (NOT eicver)
    set(eicver " _")
  endif()

  string(REGEX REPLACE ".(.*)" "\\1" ${CVER} ${eicver})

endmacro()

macro(ei_get_cxxflags VAR)
  set(${VAR} "")
  ei_is_64bit_env(IS_64BIT_ENV)
  if(EIGEN_TEST_NEON)
    set(${VAR} NEON)
  elseif(EIGEN_TEST_NEON64)
    set(${VAR} NEON)
  elseif(EIGEN_TEST_ZVECTOR)
    set(${VAR} ZVECTOR)
  elseif(EIGEN_TEST_VSX)
    set(${VAR} VSX)
  elseif(EIGEN_TEST_ALTIVEC)
    set(${VAR} ALVEC)
  elseif(EIGEN_TEST_FMA)
    set(${VAR} FMA)
  elseif(EIGEN_TEST_AVX)
    set(${VAR} AVX)
  elseif(EIGEN_TEST_SSE4_2)
    set(${VAR} SSE42)
  elseif(EIGEN_TEST_SSE4_1)
    set(${VAR} SSE41)
  elseif(EIGEN_TEST_SSSE3)
    set(${VAR} SSSE3)
  elseif(EIGEN_TEST_SSE3)
    set(${VAR} SSE3)
  elseif(EIGEN_TEST_SSE2 OR IS_64BIT_ENV)
    set(${VAR} SSE2)
  elseif(EIGEN_TEST_MSA)
    set(${VAR} MSA)
  endif()

  if(EIGEN_TEST_OPENMP)
    if (${VAR} STREQUAL "")
      set(${VAR} OMP)
    else()
      set(${VAR} ${${VAR}}-OMP)
    endif()
  endif()

  if(EIGEN_DEFAULT_TO_ROW_MAJOR)
    if (${VAR} STREQUAL "")
      set(${VAR} ROW)
    else()
      set(${VAR} ${${VAR}}-ROWMAJ)
    endif()
  endif()
endmacro()

macro(ei_set_build_string)
  ei_get_compilerver(LOCAL_COMPILER_VERSION)
  ei_get_cxxflags(LOCAL_COMPILER_FLAGS)

  set(TMP_BUILD_STRING ${CMAKE_SYSTEM}-${LOCAL_COMPILER_VERSION})

  if (NOT ${LOCAL_COMPILER_FLAGS} STREQUAL  "")
    set(TMP_BUILD_STRING ${TMP_BUILD_STRING}-${LOCAL_COMPILER_FLAGS})
  endif()

  if(EIGEN_TEST_EXTERNAL_BLAS)
    set(TMP_BUILD_STRING ${TMP_BUILD_STRING}-external_blas)
  endif()

  ei_is_64bit_env(IS_64BIT_ENV)
  if(NOT IS_64BIT_ENV)
    set(TMP_BUILD_STRING ${TMP_BUILD_STRING}-32bit)
  else()
    set(TMP_BUILD_STRING ${TMP_BUILD_STRING}-64bit)
  endif()

  if(EIGEN_BUILD_STRING_SUFFIX)
    set(TMP_BUILD_STRING ${TMP_BUILD_STRING}-${EIGEN_BUILD_STRING_SUFFIX})
  endif()

  string(TOLOWER ${TMP_BUILD_STRING} BUILDNAME)
endmacro()

macro(ei_is_64bit_env VAR)
  if(CMAKE_SIZEOF_VOID_P EQUAL 8)
    set(${VAR} 1)
  elseif(CMAKE_SIZEOF_VOID_P EQUAL 4)
    set(${VAR} 0)
  else()
    message(WARNING "Unsupported pointer size. Please contact the authors.")
  endif()
endmacro()


# helper macro for testing ei_get_compilerver_from_cxx_version_string
# STR: raw version string
# REFNAME: expected compiler name
# REFVER: expected compiler version
macro(ei_test1_get_compilerver_from_cxx_version_string STR REFNAME REFVER)
  ei_get_compilerver_from_cxx_version_string(${STR} CNAME CVER)
  if((NOT ${REFNAME} STREQUAL ${CNAME}) OR (NOT ${REFVER} STREQUAL ${CVER}))
    message("STATUS ei_get_compilerver_from_cxx_version_string error:")
    message("Expected \"${REFNAME}-${REFVER}\", got \"${CNAME}-${CVER}\"")
  endif()
endmacro()

# macro for testing ei_get_compilerver_from_cxx_version_string
# feel free to add more version strings
macro(ei_test_get_compilerver_from_cxx_version_string)
  ei_test1_get_compilerver_from_cxx_version_string("g++ (SUSE Linux) 4.5.3 20110428 [gcc-4_5-branch revision 173117]" "g++" "4.5.3")
  ei_test1_get_compilerver_from_cxx_version_string("c++ (GCC) 4.5.1 20100924 (Red Hat 4.5.1-4)" "g++" "4.5.1")
  ei_test1_get_compilerver_from_cxx_version_string("icpc (ICC) 11.0 20081105" "icpc" "11.0")
  ei_test1_get_compilerver_from_cxx_version_string("g++-3.4 (GCC) 3.4.6" "g++" "3.4.6")
  ei_test1_get_compilerver_from_cxx_version_string("SUSE Linux clang version 3.0 (branches/release_30 145598) (based on LLVM 3.0)" "llvm-clang++" "3.0")
  ei_test1_get_compilerver_from_cxx_version_string("icpc (ICC) 12.0.5 20110719" "icpc" "12.0.5")
  ei_test1_get_compilerver_from_cxx_version_string("Apple clang version 2.1 (tags/Apple/clang-163.7.1) (based on LLVM 3.0svn)" "llvm-clang++" "2.1")
  ei_test1_get_compilerver_from_cxx_version_string("i686-apple-darwin11-llvm-g++-4.2 (GCC) 4.2.1 (Based on Apple Inc. build 5658) (LLVM build 2335.15.00)" "llvm-g++" "4.2.1")
  ei_test1_get_compilerver_from_cxx_version_string("g++-mp-4.4 (GCC) 4.4.6" "g++" "4.4.6")
  ei_test1_get_compilerver_from_cxx_version_string("g++-mp-4.4 (GCC) 2011" "g++" "4.4")
  ei_test1_get_compilerver_from_cxx_version_string("x86_64-w64-mingw32-g++ (GCC) 10-win32 20210110" "mingw32-g++" "10-win32")
endmacro()

# Split all tests listed in EIGEN_TESTS_LIST into num_splits many targets
# named buildtestspartN with N = { 0, ..., num_splits-1}.
#
# The intention behind the existence of this macro is the size of Eigen's
# testsuite. Together with the relatively big compile-times building all tests
# can take a substantial amount of time depending on the available hardware.
# 
# The last buildtestspartN target will build possible remaining tests.
#
# An example:
#
#   EIGEN_TESTS_LIST= [ test1, test2, test3, test4, test5, test6, test7 ]
#
# A call to ei_split_testsuite(3) creates the following targets with dependencies
#
#   Target                      Dependencies
#   ------                      ------------
#   buildtestspart0             test1, test2
#   buildtestspart1             test3, test4
#   buildtestspart2             test5, test6, test7
#
macro(ei_split_testsuite num_splits)
  get_property(EIGEN_TESTS_LIST GLOBAL PROPERTY EIGEN_TESTS_LIST)

  # Translate EIGEN_TESTS_LIST into a CMake list
  string(REGEX REPLACE "\n" " " EIGEN_TESTS_LIST "${EIGEN_TESTS_LIST}")
  set(EIGEN_TESTS_LIST "${EIGEN_TESTS_LIST}")
  separate_arguments(EIGEN_TESTS_LIST)

  set(eigen_test_count "0")
  foreach(t IN ITEMS ${EIGEN_TESTS_LIST})
    math(EXPR eigen_test_count "${eigen_test_count}+1")
  endforeach()

  # Get number of tests per target
  math(EXPR num_tests_per_target "${eigen_test_count}/${num_splits} - ${eigen_test_count}/${num_splits} % 1")

  set(test_idx "0")
  math(EXPR target_bound "${num_splits}-1")
  foreach(part RANGE "0" "${target_bound}")
    # Create target
    set(current_target "buildtestspart${part}")
    add_custom_target("${current_target}")
    math(EXPR upper_bound "${test_idx} + ${num_tests_per_target} - 1")
    foreach(test_idx RANGE "${test_idx}" "${upper_bound}")
      list(GET EIGEN_TESTS_LIST "${test_idx}" curr_test)
      add_dependencies("${current_target}" "${curr_test}")
    endforeach()
    math(EXPR test_idx "${test_idx} + ${num_tests_per_target}")
  endforeach()
  
  # Handle the possibly remaining tests
  math(EXPR test_idx "${num_splits} * ${num_tests_per_target}")
  math(EXPR target_bound "${eigen_test_count} - 1")
  foreach(test_idx RANGE "${test_idx}" "${target_bound}")
    list(GET EIGEN_TESTS_LIST "${test_idx}" curr_test)
    add_dependencies("${current_target}" "${curr_test}")
  endforeach()
endmacro(ei_split_testsuite num_splits)

# Defines the custom command buildsmoketests to build a number of tests
# specified in smoke_test_list.
# 
# Test in smoke_test_list can be either test targets (e.g. packetmath) or
# subtests targets (e.g. packetmath_2). If any of the test are not available
# in the current configuration they are just skipped. 
#
# All tests added via this macro are labeled with the smoketest label. This
# allows running smoketests only using ctest.
#
# Smoke tests are intended to be run before the whole test suite is invoked,
# e.g., to smoke test patches.
macro(ei_add_smoke_tests smoke_test_list)
  # Set the build target to build smoketests
  set(buildtarget "buildsmoketests")
  add_custom_target("${buildtarget}")

  # Get list of all tests and translate it into a CMake list
  get_property(EIGEN_TESTS_LIST GLOBAL PROPERTY EIGEN_TESTS_LIST)
  string(REGEX REPLACE "\n" " " EIGEN_TESTS_LIST "${EIGEN_TESTS_LIST}")
  set(EIGEN_TESTS_LIST "${EIGEN_TESTS_LIST}")
  separate_arguments(EIGEN_TESTS_LIST)

  # Check if the test in smoke_test_list is a currently valid test target
  foreach(test IN ITEMS ${smoke_test_list})
    # Add tests in smoke_test_list to our smoke test target but only if the test
    # is currently available, i.e., is in EIGEN_SUBTESTS_LIST
    if ("${test}" IN_LIST EIGEN_TESTS_LIST)
      add_dependencies("${buildtarget}" "${test}")
      # In the case of a test we match all subtests
      set(ctest_regex "${ctest_regex}^${test}_[0-9]+$$|")
    endif()
  endforeach()

  # Get list of all subtests and translate it into a CMake list
  get_property(EIGEN_SUBTESTS_LIST GLOBAL PROPERTY EIGEN_SUBTESTS_LIST)
  string(REGEX REPLACE "\n" " " EIGEN_SUBTESTS_LIST "${EIGEN_SUBTESTS_LIST}")
  set(EIGEN_SUBTESTS_LIST "${EIGEN_SUBTESTS_LIST}")
  separate_arguments(EIGEN_SUBTESTS_LIST)

  # Check if the test in smoke_test_list is a currently valid subtest target
  foreach(test IN ITEMS ${smoke_test_list})
    # Add tests in smoke_test_list to our smoke test target but only if the test
    # is currently available, i.e., is in EIGEN_SUBTESTS_LIST
    if ("${test}" IN_LIST EIGEN_SUBTESTS_LIST)
      add_dependencies("${buildtarget}" "${test}")
      # Add label smoketest to be able to run smoketests using ctest
      set_property(TEST ${test} APPEND PROPERTY LABELS "smoketest")
    endif()
  endforeach()
endmacro(ei_add_smoke_tests)
