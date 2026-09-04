cmake_minimum_required(VERSION 3.10)
project(FMO)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova i pacchetti richiesti
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(MATIO REQUIRED IMPORTED_TARGET matio)

# Aggiungi l'eseguibile
add_executable(fmo tests/load-matrix.cpp)
add_executable(test-core tests/test-core.cpp)
add_executable(test-nsga2 tests/test-nsga2.cpp)
add_executable(nsga2-seq src/nsga2/seq/nsga2-seq.cpp src/nsga2/seq/nsga2-steps.cpp)
add_library(nsga2-omp STATIC src/nsga2/omp/nsga2-steps-omp.cpp)

set(FMO_TARGETS fmo test-core test-nsga2 nsga2-seq nsga2-omp)
foreach(target IN LISTS FMO_TARGETS)
	target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
endforeach()

# Collega le librerie (i percorsi di include vengono ereditati automaticamente)
target_link_libraries(fmo PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(test-core PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(test-nsga2 PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(nsga2-seq PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(nsga2-omp PUBLIC Eigen3::Eigen OpenMP::OpenMP_CXX)

add_test(NAME test-core COMMAND test-core)
add_test(NAME test-nsga2 COMMAND test-nsga2)