#!/usr/bin/env bash
# test-cli-params.sh - Verifica parametrizzazione CLI e help per eseguibili FMO
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FMO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
BUILD_DIR="${FMO_DIR}/build"

FAILURES=0

assert_success() {
    local cmd="$1"
    local desc="$2"
    if eval "$cmd" > /dev/null 2>&1; then
        echo "[PASS] $desc"
    else
        echo "[FAIL] $desc (comando: $cmd)"
        FAILURES=$((FAILURES + 1))
    fi
}

assert_failure() {
    local cmd="$1"
    local desc="$2"
    if ! eval "$cmd" > /dev/null 2>&1; then
        echo "[PASS] $desc"
    else
        echo "[FAIL] $desc (comando doveva fallire: $cmd)"
        FAILURES=$((FAILURES + 1))
    fi
}

echo "=== Test CLI Arguments & Help per FMO NSGA-II ==="

# 1. Test help exits with 0 without loading dataset
assert_success "${BUILD_DIR}/nsga2-seq -h" "nsga2-seq -h restituisce codice 0"
assert_success "${BUILD_DIR}/nsga2-seq --help" "nsga2-seq --help restituisce codice 0"
assert_success "${BUILD_DIR}/nsga2-omp -h" "nsga2-omp -h restituisce codice 0"
assert_success "${BUILD_DIR}/nsga2-omp --help" "nsga2-omp --help restituisce codice 0"
assert_success "${BUILD_DIR}/nsga2-ff -h" "nsga2-ff -h restituisce codice 0"
assert_success "${BUILD_DIR}/nsga2-ff --help" "nsga2-ff --help restituisce codice 0"

# 2. Test invalid arguments fail
assert_failure "${BUILD_DIR}/nsga2-omp" "nsga2-omp senza argomenti fallisce"
assert_failure "${BUILD_DIR}/nsga2-ff" "nsga2-ff senza argomenti fallisce"
assert_failure "${BUILD_DIR}/nsga2-ff 4 2" "nsga2-ff con modalità non valida (2) fallisce"
assert_failure "${BUILD_DIR}/nsga2-seq --pop-size 0" "nsga2-seq con pop-size 0 fallisce"
assert_failure "${BUILD_DIR}/nsga2-seq --generations -5" "nsga2-seq con generazioni negative fallisce"

# 3. Test benchmark.sh dry-run
dry_out=$("${FMO_DIR}/scripts/benchmark.sh" --dry-run)
if echo "$dry_out" | grep -q "256 50"; then
    echo "[PASS] benchmark.sh --dry-run usa di default N=256 e G=50"
else
    echo "[FAIL] benchmark.sh --dry-run non include N=256 e G=50 di default"
    FAILURES=$((FAILURES + 1))
fi

dry_custom=$("${FMO_DIR}/scripts/benchmark.sh" --dry-run -p 512 -g 20)
if echo "$dry_custom" | grep -q "512 20"; then
    echo "[PASS] benchmark.sh --dry-run accetta parametri custom -p 512 -g 20"
else
    echo "[FAIL] benchmark.sh --dry-run non imposta correttamente -p 512 -g 20"
    FAILURES=$((FAILURES + 1))
fi

if [[ $FAILURES -eq 0 ]]; then
    echo "Tutti i test CLI sono passati con successo!"
    exit 0
else
    echo "Rilevati $FAILURES fallimenti nei test CLI."
    exit 1
fi
