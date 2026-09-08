#!/usr/bin/env bash
# benchmark.sh - HPC Benchmark & Scalability Measurement for FMO (NSGA-II)
#
# Misura tempi di esecuzione, speedup, scalabilità ed efficienza delle varie
# implementazioni (Sequential, OpenMP, FastFlow parfor/farm) fino a 256 core.


set -uo pipefail

# Colori per terminale
BOLD="$(printf '\033[1m')"
RED="$(printf '\033[0;31m')"
GREEN="$(printf '\033[0;32m')"
YELLOW="$(printf '\033[0;33m')"
BLUE="$(printf '\033[0;34m')"
CYAN="$(printf '\033[0;36m')"
NC="$(printf '\033[0m')" # No Color

# Percorsi principali
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FMO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${FMO_DIR}/build"

# Configurazione di default
DEFAULT_WORKERS="1 2 4 8 16 32 64 128 256"
DEFAULT_VARIANTS="seq,omp,ff-parfor,ff-farm"
DEFAULT_REPETITIONS=3
DEFAULT_POP_SIZE=256
DEFAULT_GENERATIONS=50
DEFAULT_OUTPUT_DIR=""
BUILD_BEFORE_RUN=false
DRY_RUN=false
TIMEOUT_SEC=""

WORKERS=""
VARIANTS=""
REPETITIONS=""
POP_SIZE=""
GENERATIONS=""
OUTPUT_DIR=""

# PID del processo figlio attualmente in esecuzione (per gestione SIGINT)
CURRENT_CHILD_PID=""

# Stampa help
show_help() {
    cat <<EOF
${BOLD}USO:${NC}
  $(basename "$0") [OPZIONI]

${BOLD}DESCRIZIONE:${NC}
  Esegue benchmark automatizzati per le versioni dell'algoritmo NSGA-II in FMO:
  - nsga2-seq        (Sequenziale, baseline a 1 worker)
  - nsga2-omp        (OpenMP, variazione numero di workers)
  - nsga2-ff parfor  (FastFlow ParallelFor, modalità 0)
  - nsga2-ff farm    (FastFlow Farm, modalità 1)

  Cattura l'output TIMING di ogni fase (nsga2_total, population_evaluation, ecc.),
  salva i log raw in file separati e genera file CSV (timings_raw.csv e summary.csv)
  con calcolo automatico di Speedup assoluto (vs Seq), Speedup relativo (vs T1)
  ed Efficienza parallela.

${BOLD}OPZIONI:${NC}
  -w, --workers <LISTA>       Lista di worker/cores (default: "$DEFAULT_WORKERS")
                              Separati da spazio o virgola (es. "1 2 4 8 16 32 64 128 256")
  -v, --variants <LISTA>      Varianti da eseguire (default: "$DEFAULT_VARIANTS")
                              Valori ammessi: seq, omp, ff-parfor, ff-farm, all
  -r, --repetitions <NUM>     Numero di ripetizioni per configurazione (default: $DEFAULT_REPETITIONS)
  -p, --pop-size <NUM>        Dimensione popolazione N (default: $DEFAULT_POP_SIZE)
  -g, --generations <NUM>     Numero di generazioni G (default: $DEFAULT_GENERATIONS)
  -o, --output-dir <PATH>     Directory di output per log e CSV (default: fmo/results/run_<TIMESTAMP>)
  -b, --build                 Esegue la compilazione (scripts/build.sh) prima dei benchmark
  -t, --timeout <SEC>         Timeout in secondi per ogni singolo run (opzionale)
  -d, --dry-run               Stampa il piano dei benchmark pianificati senza eseguirli
  -h, --help                  Mostra questa schermata di aiuto

${BOLD}ESEMPI:${NC}
  # Esecuzione completa di default (seq, omp, ff-parfor, ff-farm fino a 256 core, N=256, 3 ripetizioni):
  ./scripts/benchmark.sh

  # Test rapido di verifica su OpenMP con 1 ripetizione, N=64 e 10 generazioni:
  ./scripts/benchmark.sh -v omp -w "1 2 4 8" -r 1 -p 64 -g 10

  # Test di scalabilità esteso solo per OpenMP e FastFlow Farm su macchina a 256 core con N=256:
  ./scripts/benchmark.sh -v "omp,ff-farm" -w "1 2 4 8 16 32 64 128 256" -r 5 -p 256

EOF
}

# Gestione interruzione da tastiera (Ctrl+C)
handle_sigint() {
    echo ""
    echo -e "${RED}${BOLD}[INTERRUZIONE]${NC} Ricevuto segnale SIGINT (Ctrl+C)."
    if [[ -n "$CURRENT_CHILD_PID" ]] && kill -0 "$CURRENT_CHILD_PID" 2>/dev/null; then
        echo "Terminazione del processo in esecuzione (PID $CURRENT_CHILD_PID)..."
        kill -TERM "$CURRENT_CHILD_PID" 2>/dev/null || true
        wait "$CURRENT_CHILD_PID" 2>/dev/null || true
    fi
    if [[ -d "$OUTPUT_DIR" && -f "${OUTPUT_DIR}/timings_raw.csv" ]]; then
        echo -e "${YELLOW}Tentativo di generare summary.csv con i dati parziali raccolti finora...${NC}"
        generate_summary
    fi
    echo -e "${YELLOW}Benchmark interrotto. Dati salvati in: $OUTPUT_DIR${NC}"
    exit 130
}

trap handle_sigint SIGINT SIGTERM

# Parsing argomenti CLI
while [[ $# -gt 0 ]]; do
    case "$1" in
        -w|--workers)
            WORKERS="$2"
            shift 2
            ;;
        -v|--variants)
            VARIANTS="$2"
            shift 2
            ;;
        -r|--repetitions)
            REPETITIONS="$2"
            shift 2
            ;;
        -p|--pop-size|--population)
            POP_SIZE="$2"
            shift 2
            ;;
        -g|--generations)
            GENERATIONS="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -b|--build)
            BUILD_BEFORE_RUN=true
            shift
            ;;
        -t|--timeout)
            TIMEOUT_SEC="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        *)
            echo -e "${RED}Opzione sconosciuta: $1${NC}" >&2
            show_help
            exit 1
            ;;
    esac
done

# Applicazione default
WORKERS="${WORKERS:-$DEFAULT_WORKERS}"
# Converti virgole in spazi nella lista dei workers
WORKERS=$(echo "$WORKERS" | tr ',' ' ')

VARIANTS="${VARIANTS:-$DEFAULT_VARIANTS}"
REPETITIONS="${REPETITIONS:-$DEFAULT_REPETITIONS}"
POP_SIZE="${POP_SIZE:-$DEFAULT_POP_SIZE}"
GENERATIONS="${GENERATIONS:-$DEFAULT_GENERATIONS}"

TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
if [[ -z "$OUTPUT_DIR" ]]; then
    OUTPUT_DIR="${FMO_DIR}/results/run_${TIMESTAMP}"
fi

RAW_LOGS_DIR="${OUTPUT_DIR}/raw_logs"

# Normalizza elenco varianti
SELECTED_VARIANTS=()
IFS=',' read -ra VAR_ARRAY <<< "$(echo "$VARIANTS" | tr ' ' ',')"
for v in "${VAR_ARRAY[@]}"; do
    v=$(echo "$v" | tr -d ' ' | tr '[:upper:]' '[:lower:]')
    case "$v" in
        all)
            SELECTED_VARIANTS=("seq" "omp" "ff-parfor" "ff-farm")
            break
            ;;
        seq)
            SELECTED_VARIANTS+=("seq")
            ;;
        omp)
            SELECTED_VARIANTS+=("omp")
            ;;
        ff)
            SELECTED_VARIANTS+=("ff-parfor" "ff-farm")
            ;;
        ff-parfor|parfor)
            SELECTED_VARIANTS+=("ff-parfor")
            ;;
        ff-farm|farm)
            SELECTED_VARIANTS+=("ff-farm")
            ;;
        *)
            echo -e "${RED}Variante non valida: '$v'. Varianti supportate: seq, omp, ff-parfor, ff-farm, all${NC}" >&2
            exit 1
            ;;
    esac
done

# Rimuovi duplicati dalle varianti mantenendo l'ordine
UNIQUE_VARIANTS=()
for v in "${SELECTED_VARIANTS[@]}"; do
    already=false
    for u in "${UNIQUE_VARIANTS[@]}"; do
        if [[ "$u" == "$v" ]]; then already=true; break; fi
    done
    if ! $already; then UNIQUE_VARIANTS+=("$v"); fi
done

# Compilazione su richiesta
if [[ "$BUILD_BEFORE_RUN" == true ]]; then
    echo -e "${BLUE}${BOLD}=== Compilazione del progetto (scripts/build.sh) ===${NC}"
    (cd "$FMO_DIR" && bash scripts/build.sh)
    echo ""
fi

# Verifica eseguibili
check_binary() {
    local bin="$1"
    local path="${BUILD_DIR}/${bin}"
    if [[ ! -x "$path" ]]; then
        echo -e "${RED}Errore: eseguibile non trovato o non eseguibile: ${path}${NC}" >&2
        echo -e "${YELLOW}Suggerimento: compila il progetto con 'bash scripts/build.sh' oppure usa il flag '-b'${NC}" >&2
        exit 1
    fi
}

for v in "${UNIQUE_VARIANTS[@]}"; do
    case "$v" in
        seq) check_binary "nsga2-seq" ;;
        omp) check_binary "nsga2-omp" ;;
        ff-parfor|ff-farm) check_binary "nsga2-ff" ;;
    esac
done
check_binary "compute-rmse"

# Calcolo numero totale di test previsti
TOTAL_TESTS=0
for v in "${UNIQUE_VARIANTS[@]}"; do
    if [[ "$v" == "seq" ]]; then
        TOTAL_TESTS=$((TOTAL_TESTS + REPETITIONS))
    else
        num_w=$(echo "$WORKERS" | wc -w)
        TOTAL_TESTS=$((TOTAL_TESTS + num_w * REPETITIONS))
    fi
done

# Stampa intestazione benchmark
echo -e "${CYAN}${BOLD}======================================================================${NC}"
echo -e "${CYAN}${BOLD}           FMO HPC BENCHMARK & SCALABILITY RUNNER                     ${NC}"
echo -e "${CYAN}${BOLD}======================================================================${NC}"
echo -e "${BOLD}Data di avvio:${NC}     $(date)"
echo -e "${BOLD}Varianti:${NC}          ${UNIQUE_VARIANTS[*]}"
echo -e "${BOLD}Workers testati:${NC}   $WORKERS"
echo -e "${BOLD}Ripetizioni:${NC}       $REPETITIONS"
echo -e "${BOLD}Popolazione (N):${NC}   $POP_SIZE"
echo -e "${BOLD}Generazioni (G):${NC}   $GENERATIONS"
echo -e "${BOLD}Totale esecuzioni:${NC} $TOTAL_TESTS"
echo -e "${BOLD}Output directory:${NC}  $OUTPUT_DIR"
if [[ -n "$TIMEOUT_SEC" ]]; then
    echo -e "${BOLD}Timeout per run:${NC}   ${TIMEOUT_SEC}s"
fi
echo -e "${CYAN}----------------------------------------------------------------------${NC}"

# Se dry-run, elenca solo i test e termina
if [[ "$DRY_RUN" == true ]]; then
    echo -e "${YELLOW}${BOLD}[MODALITÀ DRY-RUN ATTIVA] - I seguenti comandi verrebbero eseguiti:${NC}"
    idx=1
    for rep in $(seq 1 "$REPETITIONS"); do
        for v in "${UNIQUE_VARIANTS[@]}"; do
            if [[ "$v" == "seq" ]]; then
                echo "  [$idx/$TOTAL_TESTS] (cd $FMO_DIR && ./build/nsga2-seq $POP_SIZE $GENERATIONS) [rep=$rep]"
                idx=$((idx + 1))
            else
                for w in $WORKERS; do
                    case "$v" in
                        omp)
                            echo "  [$idx/$TOTAL_TESTS] (cd $FMO_DIR && OMP_NUM_THREADS=$w ./build/nsga2-omp $w $POP_SIZE $GENERATIONS) [rep=$rep]"
                            ;;
                        ff-parfor)
                            echo "  [$idx/$TOTAL_TESTS] (cd $FMO_DIR && ./build/nsga2-ff $w 0 $POP_SIZE $GENERATIONS) [rep=$rep]"
                            ;;
                        ff-farm)
                            echo "  [$idx/$TOTAL_TESTS] (cd $FMO_DIR && ./build/nsga2-ff $w 1 $POP_SIZE $GENERATIONS) [rep=$rep]"
                            ;;
                    esac
                    idx=$((idx + 1))
                done
            fi
        done
    done
    echo ""
    echo -e "${GREEN}Dry run completato. Nessuna operazione eseguita.${NC}"
    exit 0
fi

# Creazione cartelle di output
mkdir -p "$RAW_LOGS_DIR"
POPULATIONS_DIR="${OUTPUT_DIR}/populations"
mkdir -p "$POPULATIONS_DIR"

# Raccolta informazioni di sistema (hardware, OS, compilatore, git)
collect_system_info() {
    local info_file="${OUTPUT_DIR}/system_info.txt"
    {
        echo "======================================================================"
        echo "                    SYSTEM & ENVIRONMENT METADATA                     "
        echo "======================================================================"
        echo "Timestamp:    $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
        echo "Hostname:     $(hostname 2>/dev/null || cat /etc/hostname 2>/dev/null || echo 'unknown')"
        echo "Kernel:       $(uname -s -r -v -m)"
        echo ""
        echo "--- CPU INFORMATION ---"
        if command -v lscpu &>/dev/null; then
            lscpu
        else
            grep "model name" /proc/cpuinfo 2>/dev/null | head -n 1 || echo "CPU info unavailable"
            echo "Total logical cores: $(nproc 2>/dev/null || grep -c processor /proc/cpuinfo 2>/dev/null || echo 'unknown')"
        fi
        echo ""
        echo "--- MEMORY INFORMATION ---"
        if command -v free &>/dev/null; then
            free -h
        elif [[ -f /proc/meminfo ]]; then
            grep -E "MemTotal|MemFree|MemAvailable" /proc/meminfo
        fi
        echo ""
        echo "--- COMPILER & TOOLCHAIN ---"
        if command -v g++ &>/dev/null; then
            g++ --version | head -n 1
        fi
        if command -v cmake &>/dev/null; then
            cmake --version | head -n 1
        fi
        echo ""
        echo "--- GIT STATUS ---"
        if git rev-parse --is-inside-work-tree &>/dev/null; then
            echo "Commit: $(git rev-parse HEAD)"
            echo "Branch: $(git rev-parse --abbrev-ref HEAD)"
            echo "Status: $(git status --short)"
        else
            echo "Not a git repository"
        fi
        echo ""
        echo "--- BENCHMARK CONFIGURATION ---"
        echo "Variants:     ${UNIQUE_VARIANTS[*]}"
        echo "Workers:      $WORKERS"
        echo "Repetitions:  $REPETITIONS"
        echo "Population:   $POP_SIZE"
        echo "Generations:  $GENERATIONS"
        echo "Output dir:   $OUTPUT_DIR"
    } > "$info_file"
}

collect_system_info

# Inizializzazione del CSV dei dati grezzi di timing
TIMINGS_RAW_CSV="${OUTPUT_DIR}/timings_raw.csv"
echo "timestamp,variant,mode,workers,repetition,phase,samples,mean_ms,median_ms,wall_time_sec" > "$TIMINGS_RAW_CSV"

# Inizializzazione del CSV delle metriche di stabilità (RMSE & Spread)
STABILITY_RAW_CSV="${OUTPUT_DIR}/stability_rmse.csv"
echo "timestamp,variant,mode,workers,repetition,rmse_ptv,rmse_rectum,rmse_bladder,rmse_total_fitness,rmse_genes,spread_seq,spread_par" > "$STABILITY_RAW_CSV"

# Funzione per parsare le righe TIMING dall'output del processo
# TIMING variant=omp mode=openmp phase=population_evaluation samples=49 mean_ms=8130.47 median_ms=8308
parse_timings() {
    local log_file="$1"
    local variant="$2"
    local mode="$3"
    local workers="$4"
    local rep="$5"
    local wall_sec="$6"
    local ts
    ts="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

    while IFS= read -r line; do
        if [[ "$line" =~ ^TIMING[[:space:]]+variant=([^[:space:]]+)[[:space:]]+mode=([^[:space:]]+)[[:space:]]+phase=([^[:space:]]+)[[:space:]]+samples=([^[:space:]]+)[[:space:]]+mean_ms=([^[:space:]]+)[[:space:]]+median_ms=([^[:space:]]+) ]]; then
            local v="${BASH_REMATCH[1]}"
            local m="${BASH_REMATCH[2]}"
            local ph="${BASH_REMATCH[3]}"
            local sm="${BASH_REMATCH[4]}"
            local mn="${BASH_REMATCH[5]}"
            local md="${BASH_REMATCH[6]}"
            echo "${ts},${v},${m},${workers},${rep},${ph},${sm},${mn},${md},${wall_sec}" >> "$TIMINGS_RAW_CSV"
        fi
    done < "$log_file"
}

# Funzione per eseguire il comando con misurazione del tempo di clock
run_benchmark() {
    local variant="$1"
    local mode="$2"
    local workers="$3"
    local rep="$4"
    local current_idx="$5"

    local log_name="${variant}_${mode}_w${workers}_rep${rep}.log"
    local log_path="${RAW_LOGS_DIR}/${log_name}"

    local pop_file="${POPULATIONS_DIR}/${variant}_${mode}_w${workers}_rep${rep}.csv"
    if [[ "$variant" == "seq" ]]; then
        pop_file="${POPULATIONS_DIR}/seq_rep${rep}.csv"
    fi

    echo -ne "[${current_idx}/${TOTAL_TESTS}] ${BOLD}${variant}${NC} (mode=${mode}, workers=${workers}, rep=${rep}/${REPETITIONS})... "

    local cmd=()
    case "$variant" in
        seq)
            cmd=("./build/nsga2-seq" "$POP_SIZE" "$GENERATIONS")
            ;;
        omp)
            cmd=("./build/nsga2-omp" "$workers" "$POP_SIZE" "$GENERATIONS")
            ;;
        ff)
            local mod_flag="0"
            if [[ "$mode" == "farm" ]]; then mod_flag="1"; fi
            cmd=("./build/nsga2-ff" "$workers" "$mod_flag" "$POP_SIZE" "$GENERATIONS")
            ;;
    esac

    # Timer wall-clock preciso in nanosecondi
    local start_ns
    start_ns=$(date +%s%N)

    local exit_code=0

    # L'esecuzione DEVE avvenire nella root fmo per risolvere "../data/prostate"
    (
        cd "$FMO_DIR"
        export FMO_POPULATION_OUT="$pop_file"
        export FMO_POPULATION_SIZE="$POP_SIZE"
        export FMO_GENERATIONS="$GENERATIONS"
        if [[ "$variant" == "omp" ]]; then
            export OMP_NUM_THREADS="$workers"
        fi

        if [[ -n "$TIMEOUT_SEC" ]]; then
            timeout --kill-after=10s "${TIMEOUT_SEC}s" "${cmd[@]}"
        else
            "${cmd[@]}"
        fi
    ) > "$log_path" 2>&1 &

    CURRENT_CHILD_PID=$!
    wait "$CURRENT_CHILD_PID" || exit_code=$?
    CURRENT_CHILD_PID=""

    local end_ns
    end_ns=$(date +%s%N)
    local elapsed_ns=$((end_ns - start_ns))
    local elapsed_sec
    elapsed_sec=$(awk -v ns="$elapsed_ns" 'BEGIN { printf "%.3f", ns / 1000000000 }')

    if [[ $exit_code -eq 0 ]]; then
        parse_timings "$log_path" "$variant" "$mode" "$workers" "$rep" "$elapsed_sec"

        # Estrai tempo nsga2_total dal log se presente
        local total_ms
        total_ms=$(grep -E "^TIMING .* phase=nsga2_total " "$log_path" | awk '{ for(i=1;i<=NF;i++) if($i ~ /^mean_ms=/) { split($i,a,"="); print a[2] } }' || true)

        local ts
        ts="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

        if [[ "$variant" == "seq" ]]; then
            local spread_val
            spread_val=$(grep -E "^QUALITY variant=seq " "$log_path" | awk -F'spread_metric=' '{print $2}' | tr -d ' \r\n' || echo "")
            if [[ -z "$spread_val" ]]; then spread_val="0.0"; fi
            echo "${ts},seq,sequential,1,${rep},0.00000000,0.00000000,0.00000000,0.00000000,0.00000000,${spread_val},${spread_val}" >> "$STABILITY_RAW_CSV"

            if [[ -n "$total_ms" ]]; then
                echo -e "${GREEN}OK${NC} (${elapsed_sec}s wall, nsga2_total=${total_ms}ms, spread=${spread_val})"
            else
                echo -e "${GREEN}OK${NC} (${elapsed_sec}s wall, spread=${spread_val})"
            fi
        else
            # Calcolo RMSE vs baseline sequenziale
            local ref_seq="${POPULATIONS_DIR}/seq_rep${rep}.csv"
            if [[ ! -f "$ref_seq" ]]; then ref_seq="${POPULATIONS_DIR}/seq_rep1.csv"; fi
            if [[ ! -f "$ref_seq" && -f "${POPULATIONS_DIR}/seq_baseline.csv" ]]; then ref_seq="${POPULATIONS_DIR}/seq_baseline.csv"; fi

            local rmse_info=""
            if [[ -f "$ref_seq" && -f "$pop_file" ]]; then
                local rmse_line
                rmse_line=$("${BUILD_DIR}/compute-rmse" "$ref_seq" "$pop_file" "$variant" "$mode" "$workers" "$rep" 2>/dev/null || true)
                if [[ "$rmse_line" =~ ^RMSE ]]; then
                    local r_ptv=$(echo "$rmse_line" | awk '{for(i=1;i<=NF;i++) if($i ~ /^rmse_ptv=/) {split($i,a,"="); print a[2]}}')
                    local r_rec=$(echo "$rmse_line" | awk '{for(i=1;i<=NF;i++) if($i ~ /^rmse_rectum=/) {split($i,a,"="); print a[2]}}')
                    local r_bla=$(echo "$rmse_line" | awk '{for(i=1;i<=NF;i++) if($i ~ /^rmse_bladder=/) {split($i,a,"="); print a[2]}}')
                    local r_fit=$(echo "$rmse_line" | awk '{for(i=1;i<=NF;i++) if($i ~ /^rmse_total_fitness=/) {split($i,a,"="); print a[2]}}')
                    local r_gen=$(echo "$rmse_line" | awk '{for(i=1;i<=NF;i++) if($i ~ /^rmse_genes=/) {split($i,a,"="); print a[2]}}')
                    local s_seq=$(echo "$rmse_line" | awk '{for(i=1;i<=NF;i++) if($i ~ /^spread_seq=/) {split($i,a,"="); print a[2]}}')
                    local s_par=$(echo "$rmse_line" | awk '{for(i=1;i<=NF;i++) if($i ~ /^spread_par=/) {split($i,a,"="); print a[2]}}')
                    echo "${ts},${variant},${mode},${workers},${rep},${r_ptv},${r_rec},${r_bla},${r_fit},${r_gen},${s_seq},${s_par}" >> "$STABILITY_RAW_CSV"
                    echo "$rmse_line" >> "$log_path"
                    rmse_info=", RMSE_fit=${r_fit}, RMSE_genes=${r_gen}"
                fi
            fi

            if [[ -n "$total_ms" ]]; then
                echo -e "${GREEN}OK${NC} (${elapsed_sec}s wall, nsga2_total=${total_ms}ms${rmse_info})"
            else
                echo -e "${GREEN}OK${NC} (${elapsed_sec}s wall${rmse_info})"
            fi
        fi
    elif [[ $exit_code -eq 124 ]]; then
        echo -e "${RED}TIMEOUT (${TIMEOUT_SEC}s superati)${NC}"
    else
        echo -e "${RED}FALLITO (exit code: ${exit_code})${NC}"
        echo -e "   ${YELLOW}Dettagli errore salvati in: ${log_path}${NC}"
    fi
}

# Generazione summary aggregato (Python standard library per calcolo statistiche, speedup ed efficienza)
generate_summary() {
    local summary_file="${OUTPUT_DIR}/summary.csv"
    echo -e "\n${BLUE}${BOLD}=== Calcolo metriche di Speedup, Scalabilità ed Efficienza ===${NC}"

    python3 - <<PY_EOF
import csv
import sys
import math
import os
from collections import defaultdict

raw_csv_path = "$TIMINGS_RAW_CSV"
summary_csv_path = "$summary_file"
stability_raw_path = "$STABILITY_RAW_CSV"
output_dir = "$OUTPUT_DIR"

try:
    with open(raw_csv_path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
except Exception as e:
    print(f"Errore lettura raw CSV: {e}", file=sys.stderr)
    sys.exit(1)

if not rows:
    print("Nessun dato raccolto in timings_raw.csv")
    sys.exit(0)

# Raggruppa campioni per (variant, mode, workers, phase)
groups = defaultdict(list)
for r in rows:
    key = (r['variant'], r['mode'], int(r['workers']), r['phase'])
    try:
        mean_ms = float(r['mean_ms'])
        median_ms = float(r['median_ms'])
        wall_sec = float(r['wall_time_sec'])
        groups[key].append({
            'mean_ms': mean_ms,
            'median_ms': median_ms,
            'wall_sec': wall_sec
        })
    except ValueError:
        continue

# Calcola statistiche per gruppo
stats = {}
for key, items in groups.items():
    variant, mode, workers, phase = key
    means = [x['mean_ms'] for x in items]
    medians = [x['median_ms'] for x in items]
    count = len(means)
    avg_mean = sum(means) / count
    avg_median = sum(medians) / count
    min_val = min(means)
    max_val = max(means)
    std_mean = math.sqrt(sum((x - avg_mean) ** 2 for x in means) / count) if count > 1 else 0.0

    stats[key] = {
        'count': count,
        'avg_mean_ms': avg_mean,
        'std_mean_ms': std_mean,
        'avg_median_ms': avg_median,
        'min_ms': min_val,
        'max_ms': max_val
    }

# Cerca baseline sequenziale per ciascuna fase
seq_baselines = {}
for (variant, mode, workers, phase), data in stats.items():
    if variant == 'seq' and workers == 1:
        seq_baselines[phase] = data['avg_mean_ms']

# Cerca baseline 1-worker (T1) per ciascuna implementazione (variant, mode, phase)
t1_baselines = {}
for (variant, mode, workers, phase), data in stats.items():
    if workers == 1:
        t1_baselines[(variant, mode, phase)] = data['avg_mean_ms']

# Elaborazione metriche di stabilità (RMSE & Spread)
stab_rows = []
if os.path.exists(stability_raw_path):
    try:
        with open(stability_raw_path, 'r', newline='') as f:
            stab_rows = list(csv.DictReader(f))
    except Exception as e:
        print(f"Avviso: impossibile leggere stability_rmse.csv: {e}")

stab_groups = defaultdict(list)
for r in stab_rows:
    try:
        w_val = int(r['workers'])
        stab_groups[(r['variant'], r['mode'], w_val)].append({
            'rmse_ptv': float(r['rmse_ptv']),
            'rmse_rectum': float(r['rmse_rectum']),
            'rmse_bladder': float(r['rmse_bladder']),
            'rmse_total_fitness': float(r['rmse_total_fitness']),
            'rmse_genes': float(r['rmse_genes']),
            'spread_seq': float(r['spread_seq']),
            'spread_par': float(r['spread_par'])
        })
    except (ValueError, KeyError):
        continue

stab_stats = {}
for key, items in stab_groups.items():
    cnt = len(items)
    avg_ptv = sum(x['rmse_ptv'] for x in items) / cnt
    avg_rec = sum(x['rmse_rectum'] for x in items) / cnt
    avg_bla = sum(x['rmse_bladder'] for x in items) / cnt
    avg_fit = sum(x['rmse_total_fitness'] for x in items) / cnt
    avg_gen = sum(x['rmse_genes'] for x in items) / cnt
    avg_spread = sum(x['spread_par'] for x in items) / cnt

    if key[0] == 'seq':
        status = "BASELINE"
    elif avg_fit == 0.0 and avg_gen == 0.0:
        status = "DETERMINISTICO (0 err)"
    elif avg_fit < 1e-6:
        status = "QUASI-IDENTICO"
    else:
        status = "STABILE (Dev. GA)"

    stab_stats[key] = {
        'count': cnt,
        'avg_ptv': avg_ptv,
        'avg_rec': avg_rec,
        'avg_bla': avg_bla,
        'avg_fit': avg_fit,
        'avg_gen': avg_gen,
        'avg_spread': avg_spread,
        'status': status
    }

# Scrittura stability_summary.csv
if stab_stats:
    stab_summary_path = os.path.join(output_dir, "stability_summary.csv")
    with open(stab_summary_path, 'w', newline='') as f:
        fn = ['variant', 'mode', 'workers', 'repetitions', 'avg_rmse_ptv', 'avg_rmse_rectum', 'avg_rmse_bladder', 'avg_rmse_fitness', 'avg_rmse_genes', 'avg_spread', 'status']
        w = csv.DictWriter(f, fieldnames=fn)
        w.writeheader()
        for key in sorted(stab_stats.keys(), key=lambda k: (k[0], k[1], k[2])):
            d = stab_stats[key]
            w.writerow({
                'variant': key[0],
                'mode': key[1],
                'workers': key[2],
                'repetitions': d['count'],
                'avg_rmse_ptv': f"{d['avg_ptv']:.6f}",
                'avg_rmse_rectum': f"{d['avg_rec']:.6f}",
                'avg_rmse_bladder': f"{d['avg_bla']:.6f}",
                'avg_rmse_fitness': f"{d['avg_fit']:.6f}",
                'avg_rmse_genes': f"{d['avg_gen']:.6f}",
                'avg_spread': f"{d['avg_spread']:.6f}",
                'status': d['status']
            })

# Scrittura di summary.csv
fieldnames = [
    'variant', 'mode', 'workers', 'phase', 'repetitions',
    'avg_mean_ms', 'std_mean_ms', 'avg_median_ms', 'min_ms', 'max_ms',
    'speedup_seq', 'speedup_t1', 'efficiency_seq', 'efficiency_t1',
    'rmse_fitness', 'rmse_genes'
]

with open(summary_csv_path, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()

    # Ordina per phase, variant, mode, workers
    sorted_keys = sorted(stats.keys(), key=lambda k: (k[3], k[0], k[1], k[2]))
    for key in sorted_keys:
        variant, mode, workers, phase = key
        data = stats[key]
        tp = data['avg_mean_ms']

        # Speedup vs Seq: S_seq = T_seq / T_p
        t_seq = seq_baselines.get(phase)
        if t_seq is not None and tp > 0:
            s_seq = t_seq / tp
            eff_seq = s_seq / workers
            s_seq_str = f"{s_seq:.4f}"
            eff_seq_str = f"{eff_seq:.4f}"
        else:
            s_seq_str = ""
            eff_seq_str = ""

        # Speedup vs T1: S_t1 = T1 / T_p
        t1 = t1_baselines.get((variant, mode, phase))
        if t1 is not None and tp > 0:
            s_t1 = t1 / tp
            eff_t1 = s_t1 / workers
            s_t1_str = f"{s_t1:.4f}"
            eff_t1_str = f"{eff_t1:.4f}"
        else:
            s_t1_str = ""
            eff_t1_str = ""

        s_info = stab_stats.get((variant, mode, workers))
        rmse_fit_str = f"{s_info['avg_fit']:.6f}" if s_info else ""
        rmse_gen_str = f"{s_info['avg_gen']:.6f}" if s_info else ""

        writer.writerow({
            'variant': variant,
            'mode': mode,
            'workers': workers,
            'phase': phase,
            'repetitions': data['count'],
            'avg_mean_ms': f"{data['avg_mean_ms']:.2f}",
            'std_mean_ms': f"{data['std_mean_ms']:.2f}",
            'avg_median_ms': f"{data['avg_median_ms']:.2f}",
            'min_ms': f"{data['min_ms']:.2f}",
            'max_ms': f"{data['max_ms']:.2f}",
            'speedup_seq': s_seq_str,
            'speedup_t1': s_t1_str,
            'efficiency_seq': eff_seq_str,
            'efficiency_t1': eff_t1_str,
            'rmse_fitness': rmse_fit_str,
            'rmse_genes': rmse_gen_str
        })

print(f"File di riepilogo generato con successo: {summary_csv_path}")

# Stampa riassunto compatto a terminale per la fase principale: nsga2_total
print("\n" + "=" * 90)
print(f"{'VARIANTE':<12} {'MODE':<12} {'WORKERS':<9} {'AVG_TIME(ms)':<14} {'SPEEDUP(seq)':<14} {'SPEEDUP(T1)':<14} {'EFFICIENZA':<12}")
print("-" * 90)

for key in sorted_keys:
    variant, mode, workers, phase = key
    if phase != 'nsga2_total':
        continue
    data = stats[key]
    tp = data['avg_mean_ms']
    t_seq = seq_baselines.get(phase)
    t1 = t1_baselines.get((variant, mode, phase))
    s_seq = f"{(t_seq / tp):.2f}x" if t_seq and tp > 0 else "-"
    s_t1 = f"{(t1 / tp):.2f}x" if t1 and tp > 0 else "-"
    eff = f"{(t1 / tp / workers * 100):.1f}%" if t1 and tp > 0 else "-"

    print(f"{variant:<12} {mode:<12} {workers:<9} {tp:<14.1f} {s_seq:<14} {s_t1:<14} {eff:<12}")

print("=" * 90)

# Stampa tabella riassuntiva stabilità e determinismo delle soluzioni
if stab_stats:
    print("\n" + "=" * 105)
    print(f"  STABILITÀ DELLE SOLUZIONI (RMSE vs Baseline Sequenziale) & DETERMINISMO")
    print("=" * 105)
    print(f"{'VARIANTE':<10} {'MODE':<10} {'WORKERS':<9} {'RMSE_FITNESS':<14} {'RMSE_GENI':<14} {'SPREAD_PAR':<12} {'STATO DETERMINISMO':<28}")
    print("-" * 105)
    for key in sorted(stab_stats.keys(), key=lambda k: (k[0], k[1], k[2])):
        d = stab_stats[key]
        print(f"{key[0]:<10} {key[1]:<10} {key[2]:<9} {d['avg_fit']:<14.6f} {d['avg_gen']:<14.6f} {d['avg_spread']:<12.4f} {d['status']:<28}")
    print("=" * 105)

PY_EOF
}

# Se la variante sequenziale non è inclusa nelle varianti selezionate,
# esegui una baseline preliminare per abilitare il calcolo di RMSE
HAS_SEQ=false
for v in "${UNIQUE_VARIANTS[@]}"; do
    if [[ "$v" == "seq" ]]; then HAS_SEQ=true; break; fi
done

if [[ "$HAS_SEQ" == false && "$DRY_RUN" == false ]]; then
    echo -e "${YELLOW}${BOLD}Nota: 'seq' non inclusa in -v. Esecuzione baseline sequenziale per calcolo RMSE...${NC}"
    (
        cd "$FMO_DIR"
        export FMO_POPULATION_OUT="${POPULATIONS_DIR}/seq_baseline.csv"
        export FMO_POPULATION_SIZE="$POP_SIZE"
        export FMO_GENERATIONS="$GENERATIONS"
        ./build/nsga2-seq "$POP_SIZE" "$GENERATIONS" > "${RAW_LOGS_DIR}/seq_baseline.log" 2>&1
    )
    if [[ -f "${POPULATIONS_DIR}/seq_baseline.csv" ]]; then
        echo -e "${GREEN}✓ Baseline sequenziale salvata in: ${POPULATIONS_DIR}/seq_baseline.csv${NC}\n"
    fi
fi

# Esecuzione del loop principale
current_test=1
for rep in $(seq 1 "$REPETITIONS"); do
    echo ""
    echo -e "${YELLOW}${BOLD}>>> RIPETIZIONE ${rep}/${REPETITIONS} <<<${NC}"
    for v in "${UNIQUE_VARIANTS[@]}"; do
        case "$v" in
            seq)
                run_benchmark "seq" "sequential" 1 "$rep" "$current_test"
                current_test=$((current_test + 1))
                ;;
            omp)
                for w in $WORKERS; do
                    run_benchmark "omp" "openmp" "$w" "$rep" "$current_test"
                    current_test=$((current_test + 1))
                done
                ;;
            ff-parfor)
                for w in $WORKERS; do
                    run_benchmark "ff" "parfor" "$w" "$rep" "$current_test"
                    current_test=$((current_test + 1))
                done
                ;;
            ff-farm)
                for w in $WORKERS; do
                    run_benchmark "ff" "farm" "$w" "$rep" "$current_test"
                    current_test=$((current_test + 1))
                done
                ;;
        esac
    done
done

# Calcolo summary
generate_summary

# Avvio automatico dello script di plotting se disponibile
PLOT_SCRIPT="${SCRIPT_DIR}/plot_benchmark.py"
if [[ -f "$PLOT_SCRIPT" ]]; then
    echo -e "\n${BLUE}${BOLD}=== Generazione automatica grafici e report HTML ===${NC}"
    python3 "$PLOT_SCRIPT" "$OUTPUT_DIR" || echo -e "${YELLOW}Generazione grafici saltata (eseguibile manualmente con: python3 scripts/plot_benchmark.py $OUTPUT_DIR)${NC}"
fi

echo ""
echo -e "${GREEN}${BOLD}======================================================================${NC}"
echo -e "${GREEN}${BOLD}                BENCHMARK COMPLETATO CON SUCCESSO!                    ${NC}"
echo -e "${GREEN}${BOLD}======================================================================${NC}"
echo -e "Tutti i risultati sono stati salvati in:"
echo -e "  ${BOLD}${OUTPUT_DIR}${NC}"
echo -e "  ├── ${BOLD}timings_raw.csv${NC}       (Tutti i campioni di timing misurati per fase)"
echo -e "  ├── ${BOLD}summary.csv${NC}           (Riepilogo statistiche, speedup ed efficienza)"
echo -e "  ├── ${BOLD}stability_rmse.csv${NC}    (Campioni di RMSE e Spread vs baseline)"
echo -e "  ├── ${BOLD}stability_summary.csv${NC} (Riepilogo discrepanze e determinismo)"
echo -e "  ├── ${BOLD}system_info.txt${NC}       (Metadati CPU, NUMA, OS, toolchain)"
echo -e "  ├── ${BOLD}populations/${NC}          (Popolazioni finali esportate in CSV)"
echo -e "  └── ${BOLD}raw_logs/${NC}             (File di log completi per ogni singola esecuzione)"
echo ""
echo -e "Per visualizzare o rigenerare i grafici in qualsiasi momento:"
echo -e "  ${CYAN}python3 scripts/plot_benchmark.py ${OUTPUT_DIR}${NC}"
echo ""
