#!/bin/bash
set -e

# Pulisci la build precedente
echo "Pulizia della build precedente..."
rm -rf build/

# Crea e entra in build
echo "Creazione della directory build..."
mkdir build
cd build

# Configura e compila
echo "Esecuzione di cmake..."
cmake ..

echo "Compilazione in corso..."
make

echo ""
echo "Build completed successfully, execute in: build/fmo !"
