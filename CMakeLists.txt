cmake_minimum_required(VERSION 3.10)
project(FMO)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova i pacchetti richiesti
find_package(Eigen3 REQUIRED)
find_package(PkgConfig REQUIRED)
pkg_check_modules(MATIO REQUIRED IMPORTED_TARGET matio)

# Aggiungi l'eseguibile
add_executable(fmo preprocessing.cpp)

# Collega le librerie (i percorsi di include vengono ereditati automaticamente)
target_link_libraries(fmo PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)