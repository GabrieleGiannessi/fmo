#!/bin/bash
set -e

# Script per scaricare i dataset di matrici di influenza

# URL dei dataset
TG119_URL="https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100110/TG119.zip"
PROSTATE_URL="https://s3.ap-northeast-1.wasabisys.com/gigadb-datasets/live/pub/10.5524/100001_101000/100110/PROSTATE.zip"

# Directory di destinazione (allo stesso livello delle librerie fastflow, eigen, matio e del progetto fmo)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"
DATA_DIR="$PROJECT_ROOT/data"

# Colori per output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${YELLOW}================================${NC}"
echo -e "${YELLOW}Download Influence Matrices${NC}"
echo -e "${YELLOW}================================${NC}"
echo ""

# Crea la directory principale e le sottocartelle se non esistono
mkdir -p "$DATA_DIR/phantom"
mkdir -p "$DATA_DIR/prostate"

# Funzione per scaricare e estrarre
download_and_extract() {
    local url=$1
    local filename=$2
    local target_dir=$3
    local filepath="$target_dir/$filename"
    local folder_name
    folder_name="$(basename "$target_dir")"
    
    echo -e "${YELLOW}Scaricamento: $filename in $folder_name/${NC}"
    
    # Scarica il file
    if curl -L -o "$filepath" "$url"; then
        echo -e "${GREEN}✓ Download completato: $filename${NC}"
        
        # Estrai il file nella sottocartella appropriata
        echo -e "${YELLOW}Estrazione: $filename in $target_dir${NC}"
        unzip -q -o "$filepath" -d "$target_dir"
        
        if [ $? -eq 0 ]; then
            echo -e "${GREEN}✓ Estrazione completata${NC}"
            # Rimuovi il file zip dopo l'estrazione
            rm -f "$filepath"
            echo -e "${GREEN}✓ File ZIP rimosso${NC}"
            
            # Se lo zip contiene una sottocartella (es. TG119 o PROSTATE), sposta i file al livello principale
            local name_without_ext="${filename%.*}"
            if [ -d "$target_dir/$name_without_ext" ]; then
                mv "$target_dir/$name_without_ext"/* "$target_dir/" 2>/dev/null || true
                rmdir "$target_dir/$name_without_ext" 2>/dev/null || true
            fi
        else
            echo -e "${RED}✗ Errore durante l'estrazione di $filename${NC}"
            rm -f "$filepath"
            exit 1
        fi
    else
        echo -e "${RED}✗ Errore durante il download di $filename${NC}"
        rm -f "$filepath"
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

# Scarica TG119 nella cartella phantom
download_and_extract "$TG119_URL" "TG119.zip" "$DATA_DIR/phantom"

# Scarica PROSTATE nella cartella prostate
download_and_extract "$PROSTATE_URL" "PROSTATE.zip" "$DATA_DIR/prostate"

echo -e "${GREEN}✓ Download completato con successo!${NC}"
echo -e "${GREEN}I dati sono in: $DATA_DIR${NC}"

# Mostra il contenuto della directory
echo ""
echo "Contenuto della directory $DATA_DIR:"
ls -lh "$DATA_DIR"
echo ""
echo "File in phantom: $(ls -1 "$DATA_DIR/phantom" 2>/dev/null | wc -l)"
echo "File in prostate: $(ls -1 "$DATA_DIR/prostate" 2>/dev/null | wc -l)"
