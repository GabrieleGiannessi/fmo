cmake_minimum_required(VERSION 3.10)
project(FMO)

include(CTest)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Trova i pacchetti richiesti
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)
pkg_check_modules(MATIO REQUIRED IMPORTED_TARGET matio)

# FastFlow si trova nella directory imrt, due livelli sopra questa directory.
set(FASTFLOW_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../fastflow")

# Aggiungi l'eseguibile
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
target_include_directories(nsga2-ff PRIVATE ${FASTFLOW_DIR})

# Collega le librerie (i percorsi di include vengono ereditati automaticamente)
target_link_libraries(fmo PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(test-core PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(test-nsga2 PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(nsga2-seq PRIVATE ${MATIO_LIBRARIES} Eigen3::Eigen)
target_link_libraries(nsga2-omp PRIVATE
	${MATIO_LIBRARIES}
	Eigen3::Eigen
	OpenMP::OpenMP_CXX
)
target_link_libraries(nsga2-ff PRIVATE
	${MATIO_LIBRARIES}
	Eigen3::Eigen
	Threads::Threads
)


add_test(NAME test-core COMMAND test-core)
add_test(NAME test-nsga2 COMMAND test-nsga2)