# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(SeismicSimulation)

# Especifica o padrão C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Lista de arquivos fonte
set(SOURCES
    main.cpp
    fft.cpp
    eigenValueSolver.cpp   
    seismicSimulator.cpp
    gnuPlotter.cpp
)

# Diretórios para a FFTW3
set(FFTW3_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/extern/")
set(FFTW3_LIBRARY_DIR "${CMAKE_SOURCE_DIR}/extern/fftw3")

# Inclui o diretório de cabeçalhos do FFTW3
include_directories(${FFTW3_INCLUDE_DIR})

add_subdirectory(extern/eigen)

# Adiciona o executável
add_executable(${PROJECT_NAME} ${SOURCES})

# Linka a biblioteca FFTW3 ao projeto
if (WIN32)
    # Usar o arquivo .dll.a para linkagem com MinGW
    set(FFTW3_LIBRARY "${FFTW3_LIBRARY_DIR}/libfftw3-3.dll.a")
    target_link_libraries(${PROJECT_NAME} ${FFTW3_LIBRARY} eigen)

    # Cópia da DLL para o diretório de saída
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${FFTW3_LIBRARY_DIR}/libfftw3-3.dll"
        $<TARGET_FILE_DIR:${PROJECT_NAME}>)
endif()

# Inclui diretórios adicionais, se necessário
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_SOURCE_DIR})
