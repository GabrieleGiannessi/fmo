#!/bin/bash
set -e

# Pulisci la build precedente
echo "Cleaning the previous build dir..."
rm -rf build/

# Crea e entra in build
echo "Creation of directory build..."
mkdir build
cd build

# Configura e compila
echo "cmake execution..."
cmake ..

echo "Compilation in progress..."
make

echo ""
echo "Build completed successfully, the files are in the "build/" directory!"
