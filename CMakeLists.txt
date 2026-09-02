cmake_minimum_required(VERSION 3.10)
project(FMO)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova i pacchetti richiesti
find_package(Eigen3 REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(MATIO REQUIRED IMPORTED_TARGET matio)

# Aggiungi il percorso di include per le funzioni helper
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Aggiungi l'eseguibile
add_executable(fmo tests/load-matrix.cpp)
add_executable(test-core tests/test-core.cpp)

# Collega le librerie (i percorsi di include vengono ereditati automaticamente)
target_link_libraries(fmo PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(test-core PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)