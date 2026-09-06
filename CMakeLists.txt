cmake_minimum_required(VERSION 3.14)
project(FMO)

include(CTest)
include(CheckCXXCompilerFlag)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Gestione compatibilità vecchie versioni di Eigen con GCC >= 14
check_cxx_compiler_flag("-Wno-template-body" COMPILER_SUPPORTS_WNO_TEMPLATE_BODY)
if(COMPILER_SUPPORTS_WNO_TEMPLATE_BODY)
    add_compile_options($<$<COMPILE_LANGUAGE:CXX>:-Wno-template-body>)
endif()

# Opzione per dare precedenza ai repository git clonati in locale
option(FMO_USE_LOCAL_DEPS "Prioritizza i cloni git locali per FastFlow, Eigen e Matio" ON)

# ==============================================================================
# 1. FastFlow (Header-only)
# ==============================================================================
set(FASTFLOW_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../fastflow" CACHE PATH "Percorso della cartella FastFlow")
if(NOT EXISTS "${FASTFLOW_DIR}/ff")
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../fastflow/ff")
        set(FASTFLOW_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../fastflow")
    endif()
endif()

if(EXISTS "${FASTFLOW_DIR}/ff")
    message(STATUS "FastFlow trovato in: ${FASTFLOW_DIR}")
    add_library(fastflow INTERFACE)
    target_include_directories(fastflow INTERFACE "${FASTFLOW_DIR}")
else()
    message(WARNING "Cartella FastFlow non trovata in ${FASTFLOW_DIR}")
endif()

# ==============================================================================
# 2. Eigen3 (Header-only)
# ==============================================================================
set(EIGEN3_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../eigen" CACHE PATH "Percorso della cartella Eigen3")

if(FMO_USE_LOCAL_DEPS AND EXISTS "${EIGEN3_ROOT}/Eigen/Dense")
    message(STATUS "Utilizzo Eigen locale da: ${EIGEN3_ROOT}")
    add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
    target_include_directories(Eigen3::Eigen INTERFACE "${EIGEN3_ROOT}")
else()
    find_package(Eigen3 QUIET)
    if(NOT TARGET Eigen3::Eigen)
        if(EXISTS "${EIGEN3_ROOT}/Eigen/Dense")
            message(STATUS "Eigen3 di sistema non trovato: utilizzo Eigen locale da: ${EIGEN3_ROOT}")
            add_library(Eigen3::Eigen INTERFACE IMPORTED GLOBAL)
            target_include_directories(Eigen3::Eigen INTERFACE "${EIGEN3_ROOT}")
        else()
            find_package(Eigen3 REQUIRED)
        endif()
    else()
        message(STATUS "Utilizzo Eigen3 di sistema")
    endif()
endif()

# ==============================================================================
# 3. Matio (Libreria C compilata)
# ==============================================================================
set(MATIO_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../matio" CACHE PATH "Percorso della cartella matio")

if(FMO_USE_LOCAL_DEPS AND EXISTS "${MATIO_DIR}/CMakeLists.txt")
    set(USE_LOCAL_MATIO TRUE)
else()
    set(USE_LOCAL_MATIO FALSE)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
        pkg_check_modules(MATIO QUIET IMPORTED_TARGET matio)
    endif()
    if(NOT TARGET PkgConfig::MATIO AND EXISTS "${MATIO_DIR}/CMakeLists.txt")
        set(USE_LOCAL_MATIO TRUE)
    endif()
endif()

if(USE_LOCAL_MATIO)
    message(STATUS "Compilazione locale di matio come sottoprogetto da: ${MATIO_DIR}")
    set(MATIO_BUILD_TESTING OFF CACHE BOOL "" FORCE)

    # Se HDF5 non è disponibile sul sistema remoto, disabilita il supporto MAT73.
    # I dataset del progetto (formato MATLAB v5) necessitano solo di zlib.
    find_package(HDF5 QUIET)
    if(NOT HDF5_FOUND)
        message(STATUS "HDF5 non trovato: disabilitato supporto MAT73 in matio (formato v5 pienamente supportato)")
        set(MATIO_MAT73 OFF CACHE BOOL "" FORCE)
        set(MATIO_WITH_HDF5 OFF CACHE BOOL "" FORCE)
    endif()

    # Compilazione statica: evita problemi di LD_LIBRARY_PATH a runtime
    set(MATIO_SHARED OFF CACHE BOOL "" FORCE)
    set(MATIO_PIC ON CACHE BOOL "" FORCE)

    add_subdirectory("${MATIO_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/matio" EXCLUDE_FROM_ALL)
    set(MATIO_TARGET matio::matio)
elseif(TARGET PkgConfig::MATIO)
    set(MATIO_TARGET PkgConfig::MATIO)
    message(STATUS "Utilizzo matio di sistema (PkgConfig)")
else()
    message(FATAL_ERROR "matio non trovato né come pacchetto di sistema né nella cartella ${MATIO_DIR}")
endif()

# ==============================================================================
# 4. Altre dipendenze standard
# ==============================================================================
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)

# ==============================================================================
# 5. Definizione degli eseguibili
# ==============================================================================
add_executable(fmo tests/load-matrix.cpp)
add_executable(test-core tests/test-core.cpp)
add_executable(test-nsga2 tests/test-nsga2.cpp)
add_executable(nsga2-seq
	src/nsga2/seq/engine-seq.cpp
	src/nsga2/seq/main.cpp
	include/fmo/utilities/utimer.cpp
)
add_executable(nsga2-omp
	src/nsga2/omp/engine-omp.cpp
	src/nsga2/omp/main.cpp
	include/fmo/utilities/utimer.cpp
)
add_executable(nsga2-ff
	src/nsga2/ff/engine-ff.cpp
	src/nsga2/ff/main.cpp
	src/nsga2/seq/engine-seq.cpp
	include/fmo/utilities/utimer.cpp
)

set(FMO_TARGETS fmo test-core test-nsga2 nsga2-seq nsga2-omp nsga2-ff)
foreach(target IN LISTS FMO_TARGETS)
	target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
endforeach()

if(TARGET fastflow)
	target_link_libraries(nsga2-ff PRIVATE fastflow)
else()
	target_include_directories(nsga2-ff PRIVATE ${FASTFLOW_DIR})
endif()

# ==============================================================================
# 6. Linking delle librerie
# ==============================================================================
target_link_libraries(fmo PRIVATE ${MATIO_TARGET} Eigen3::Eigen)
target_link_libraries(test-core PRIVATE ${MATIO_TARGET} Eigen3::Eigen)
target_link_libraries(test-nsga2 PRIVATE ${MATIO_TARGET} Eigen3::Eigen)
target_link_libraries(nsga2-seq PRIVATE ${MATIO_TARGET} Eigen3::Eigen)
target_link_libraries(nsga2-omp PRIVATE
	${MATIO_TARGET}
	Eigen3::Eigen
	OpenMP::OpenMP_CXX
)
target_link_libraries(nsga2-ff PRIVATE
	${MATIO_TARGET}
	Eigen3::Eigen
	Threads::Threads
)

add_test(NAME test-core COMMAND test-core)
add_test(NAME test-nsga2 COMMAND test-nsga2)