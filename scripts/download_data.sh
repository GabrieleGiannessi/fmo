#!/bin/bash
set -e

# Script per scaricare i dataset di matrici di influenza

# URL dei dataset
TG119_URL="https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100110/TG119.zip"
PROSTATE_URL="https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100110/PROSTATE.zip"

# Directory di destinazione
DATA_DIR="../data"

# Colori per output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}================================${NC}"
echo -e "${YELLOW}Download Influence Matrices${NC}"
echo -e "${YELLOW}================================${NC}"
echo ""

# Crea la directory se non esiste
if [ ! -d "$DATA_DIR" ]; then
    echo "Creazione directory: $DATA_DIR"
    mkdir -p "$DATA_DIR"
fi

# Funzione per scaricare e estrarre
download_and_extract() {
    local url=$1
    local filename=$2
    local filepath="$DATA_DIR/$filename"
    
    echo -e "${YELLOW}Scaricamento: $filename${NC}"
    
    # Scarica il file
    if curl -L -o "$filepath" "$url"; then
        echo -e "${GREEN}✓ Download completato: $filename${NC}"
        
        # Estrai il file
        echo -e "${YELLOW}Estrazione: $filename${NC}"
        unzip -q "$filepath" -d "$DATA_DIR"
        
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ Estrazione completata${NC}"
            # Rimuovi il file zip dopo l'estrazione
            rm "$filepath"
            echo -e "${GREEN}✓ File ZIP rimosso${NC}"
        else
            echo -e "${RED}✗ Errore durante l'estrazione di $filename${NC}"
            exit 1
        fi
    else
        echo -e "${RED}✗ Errore durante il download di $filename${NC}"
        exit 1
    fi
    
    echo ""
}

# Controlla se curl è disponibile
if ! command -v curl &> /dev/null; then
    echo -e "${RED}✗ curl non è installato. Installalo con: sudo apt install curl${NC}"
    exit 1
fi

# Controlla se unzip è disponibile
if ! command -v unzip &> /dev/null; then
    echo -e "${RED}✗ unzip non è installato. Installalo con: sudo apt install unzip${NC}"
    exit 1
fi

# Scarica TG119
download_and_extract "$TG119_URL" "TG119.zip"

# Scarica PROSTATE
download_and_extract "$PROSTATE_URL" "PROSTATE.zip"

echo -e "${GREEN}✓ Download completato con successo!${NC}"
echo -e "${GREEN}I dati sono in: $DATA_DIR${NC}"

# Mostra il contenuto della directory
echo ""
echo "Contenuto della directory:"
ls -lh "$DATA_DIR"
