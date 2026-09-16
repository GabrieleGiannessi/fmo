#!/usr/bin/env python3
import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import time
from collections import defaultdict


def backup_directory(base_dir):
    """Crea una copia di backup dei file CSV e grafici prima di effettuare il merge."""
    ts = time.strftime("%Y%m%d_%H%M%S")
    backup_dir = os.path.join(base_dir, f"backup_{ts}")
    os.makedirs(backup_dir, exist_ok=True)

    files_to_backup = [
        "summary.csv",
        "timings_raw.csv",
        "stability_rmse.csv",
        "stability_summary.csv",
    ]
    for fn in files_to_backup:
        src = os.path.join(base_dir, fn)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(backup_dir, fn))

    print(f"Backup salvato in: {backup_dir}")
    return backup_dir


def load_csv(path):
    if not os.path.exists(path):
        return []
    with open(path, "r", newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)


def save_csv(path, rows, fieldnames):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            filtered_r = {k: r.get(k, "") for k in fieldnames}
            writer.writerow(filtered_r)


def recalculate_summary(output_dir):
    """
    Ricalcola summary.csv e stability_summary.csv da timings_raw.csv e stability_rmse.csv.
    Utilizza la stessa logica consolidata di benchmark.sh.
    """
    raw_csv_path = os.path.join(output_dir, "timings_raw.csv")
    summary_csv_path = os.path.join(output_dir, "summary.csv")
    stability_raw_path = os.path.join(output_dir, "stability_rmse.csv")

    raw_rows = load_csv(raw_csv_path)
    if not raw_rows:
        print("Nessun dato trovato in timings_raw.csv", file=sys.stderr)
        return

    # Raggruppa campioni per (variant, mode, workers, phase)
    groups = defaultdict(list)
    for r in raw_rows:
        key = (r["variant"], r["mode"], int(r["workers"]), r["phase"])
        try:
            mean_ms = float(r["mean_ms"])
            median_ms = float(r["median_ms"])
            wall_sec = float(r.get("wall_time_sec", 0.0))
            groups[key].append(
                {"mean_ms": mean_ms, "median_ms": median_ms, "wall_sec": wall_sec}
            )
        except (ValueError, KeyError):
            continue

    # Calcola statistiche aggregate
    stats = {}
    for key, items in groups.items():
        means = [x["mean_ms"] for x in items]
        medians = [x["median_ms"] for x in items]
        count = len(means)
        avg_mean = sum(means) / count
        avg_median = sum(medians) / count
        min_val = min(means)
        max_val = max(means)
        std_mean = (
            math.sqrt(sum((x - avg_mean) ** 2 for x in means) / count)
            if count > 1
            else 0.0
        )

        stats[key] = {
            "count": count,
            "avg_mean_ms": avg_mean,
            "std_mean_ms": std_mean,
            "avg_median_ms": avg_median,
            "min_ms": min_val,
            "max_ms": max_val,
        }

    # Baseline sequenziale (variant=seq, workers=1)
    seq_baselines = {}
    for (variant, mode, workers, phase), data in stats.items():
        if variant == "seq" and workers == 1:
            seq_baselines[phase] = data["avg_mean_ms"]

    # Baseline T1 per ogni implementazione
    t1_baselines = {}
    for (variant, mode, workers, phase), data in stats.items():
        if workers == 1:
            t1_baselines[(variant, mode, phase)] = data["avg_mean_ms"]

    # Stability stats
    stab_rows = load_csv(stability_raw_path)
    stab_groups = defaultdict(list)
    for r in stab_rows:
        try:
            w_val = int(r["workers"])
            stab_groups[(r["variant"], r["mode"], w_val)].append(
                {
                    "rmse_ptv": float(r.get("rmse_ptv", 0.0)),
                    "rmse_rectum": float(r.get("rmse_rectum", 0.0)),
                    "rmse_bladder": float(r.get("rmse_bladder", 0.0)),
                    "rmse_total_fitness": float(r.get("rmse_total_fitness", 0.0)),
                    "rmse_genes": float(r.get("rmse_genes", 0.0)),
                    "spread_seq": float(r.get("spread_seq", 0.0)),
                    "spread_par": float(r.get("spread_par", 0.0)),
                }
            )
        except (ValueError, KeyError):
            continue

    stab_stats = {}
    for key, items in stab_groups.items():
        cnt = len(items)
        avg_ptv = sum(x["rmse_ptv"] for x in items) / cnt
        avg_rec = sum(x["rmse_rectum"] for x in items) / cnt
        avg_bla = sum(x["rmse_bladder"] for x in items) / cnt
        avg_fit = sum(x["rmse_total_fitness"] for x in items) / cnt
        avg_gen = sum(x["rmse_genes"] for x in items) / cnt
        avg_spread = sum(x["spread_par"] for x in items) / cnt

        if key[0] == "seq":
            status = "BASELINE"
        elif avg_fit == 0.0 and avg_gen == 0.0:
            status = "DETERMINISTICO (0 err)"
        elif avg_fit < 1e-6:
            status = "QUASI-IDENTICO"
        else:
            status = "STABILE (Dev. GA)"

        stab_stats[key] = {
            "count": cnt,
            "avg_ptv": avg_ptv,
            "avg_rec": avg_rec,
            "avg_bla": avg_bla,
            "avg_fit": avg_fit,
            "avg_gen": avg_gen,
            "avg_spread": avg_spread,
            "status": status,
        }

    # Salva stability_summary.csv se presenti dati di stabilità
    if stab_stats:
        stab_summary_path = os.path.join(output_dir, "stability_summary.csv")
        fn = [
            "variant",
            "mode",
            "workers",
            "repetitions",
            "avg_rmse_ptv",
            "avg_rmse_rectum",
            "avg_rmse_bladder",
            "avg_rmse_fitness",
            "avg_rmse_genes",
            "avg_spread",
            "status",
        ]
        with open(stab_summary_path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fn)
            w.writeheader()
            for key in sorted(stab_stats.keys(), key=lambda k: (k[0], k[1], k[2])):
                d = stab_stats[key]
                w.writerow(
                    {
                        "variant": key[0],
                        "mode": key[1],
                        "workers": key[2],
                        "repetitions": d["count"],
                        "avg_rmse_ptv": f"{d['avg_ptv']:.6f}",
                        "avg_rmse_rectum": f"{d['avg_rec']:.6f}",
                        "avg_rmse_bladder": f"{d['avg_bla']:.6f}",
                        "avg_rmse_fitness": f"{d['avg_fit']:.6f}",
                        "avg_rmse_genes": f"{d['avg_gen']:.6f}",
                        "avg_spread": f"{d['avg_spread']:.6f}",
                        "status": d["status"],
                    }
                )

    # Scrivi summary.csv
    fieldnames = [
        "variant",
        "mode",
        "workers",
        "phase",
        "repetitions",
        "avg_mean_ms",
        "std_mean_ms",
        "avg_median_ms",
        "min_ms",
        "max_ms",
        "speedup_seq",
        "speedup_t1",
        "efficiency_seq",
        "efficiency_t1",
        "rmse_fitness",
        "rmse_genes",
    ]

    with open(summary_csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        sorted_keys = sorted(stats.keys(), key=lambda k: (k[3], k[0], k[1], k[2]))
        for key in sorted_keys:
            variant, mode, workers, phase = key
            data = stats[key]
            tp = data["avg_mean_ms"]

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

            writer.writerow(
                {
                    "variant": variant,
                    "mode": mode,
                    "workers": workers,
                    "phase": phase,
                    "repetitions": data["count"],
                    "avg_mean_ms": f"{tp:.2f}",
                    "std_mean_ms": f"{data['std_mean_ms']:.2f}",
                    "avg_median_ms": f"{data['avg_median_ms']:.2f}",
                    "min_ms": f"{data['min_ms']:.2f}",
                    "max_ms": f"{data['max_ms']:.2f}",
                    "speedup_seq": s_seq_str,
                    "speedup_t1": s_t1_str,
                    "efficiency_seq": eff_seq_str,
                    "efficiency_t1": eff_t1_str,
                    "rmse_fitness": rmse_fit_str,
                    "rmse_genes": rmse_gen_str,
                }
            )

    print(f"File summary.csv ricalcolato con successo in: {summary_csv_path}")


def merge_benchmark_runs(base_dir, new_dir, target_variant="ff", target_mode="parfor", keep_old_as=None):
    """
    Effettua il merge dei dati raw da new_dir dentro base_dir:
    - Se keep_old_as è None: le vecchie righe (target_variant, target_mode) in base_dir vengono SOVRASCRITTE.
    - Se keep_old_as è specificato (es. 'parfor-naive'): le vecchie righe vengono rinominate con tale modalità.
    """
    base_raw_path = os.path.join(base_dir, "timings_raw.csv")
    new_raw_path = os.path.join(new_dir, "timings_raw.csv")

    if not os.path.exists(base_raw_path):
        print(f"Errore: {base_raw_path} non trovato.", file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(new_raw_path):
        print(f"Errore: {new_raw_path} non trovato.", file=sys.stderr)
        sys.exit(1)

    backup_directory(base_dir)

    base_raw = load_csv(base_raw_path)
    new_raw = load_csv(new_raw_path)

    # Filtra solo i campioni rilevanti dalla nuova run
    new_filtered = [
        r
        for r in new_raw
        if r.get("variant") == target_variant and r.get("mode") == target_mode
    ]

    if not new_filtered:
        print(
            f"Attenzione: nessun record con variant='{target_variant}' e mode='{target_mode}' trovato in {new_raw_path}.",
            file=sys.stderr,
        )
        print("Trovate invece le seguenti combinazioni nella nuova run:")
        combos = {(r.get("variant"), r.get("mode")) for r in new_raw}
        for v, m in combos:
            print(f"  - variant={v}, mode={m}")
        sys.exit(1)

    merged_raw = []
    old_count = 0
    for r in base_raw:
        if r.get("variant") == target_variant and r.get("mode") == target_mode:
            old_count += 1
            if keep_old_as:
                r_copy = dict(r)
                r_copy["mode"] = keep_old_as
                merged_raw.append(r_copy)
        else:
            merged_raw.append(r)

    # Aggiunge i nuovi record
    merged_raw.extend(new_filtered)

    if keep_old_as:
        print(
            f"Rinominati {old_count} vecchi record come mode='{keep_old_as}'. Aggiunti {len(new_filtered)} nuovi record per mode='{target_mode}'."
        )
    else:
        print(
            f"Sovrascritti {old_count} vecchi record con {len(new_filtered)} nuovi record per {target_variant}/{target_mode}."
        )

    # Salva il nuovo timings_raw.csv unificato
    fieldnames = list(base_raw[0].keys())
    save_csv(base_raw_path, merged_raw, fieldnames)

    # Merge stability_rmse.csv se presente
    base_stab_path = os.path.join(base_dir, "stability_rmse.csv")
    new_stab_path = os.path.join(new_dir, "stability_rmse.csv")
    if os.path.exists(new_stab_path) and os.path.exists(base_stab_path):
        base_stab = load_csv(base_stab_path)
        new_stab = load_csv(new_stab_path)
        new_stab_filtered = [
            r
            for r in new_stab
            if r.get("variant") == target_variant and r.get("mode") == target_mode
        ]
        merged_stab = []
        for r in base_stab:
            if r.get("variant") == target_variant and r.get("mode") == target_mode:
                if keep_old_as:
                    r_copy = dict(r)
                    r_copy["mode"] = keep_old_as
                    merged_stab.append(r_copy)
            else:
                merged_stab.append(r)
        merged_stab.extend(new_stab_filtered)
        if base_stab:
            save_csv(base_stab_path, merged_stab, list(base_stab[0].keys()))

    # Copia i file di log grezzi da new_dir/raw_logs a base_dir/raw_logs
    new_logs_dir = os.path.join(new_dir, "raw_logs")
    base_logs_dir = os.path.join(base_dir, "raw_logs")
    if os.path.exists(new_logs_dir) and os.path.exists(base_logs_dir):
        for fname in os.listdir(new_logs_dir):
            if target_mode in fname and fname.endswith(".log"):
                shutil.copy2(
                    os.path.join(new_logs_dir, fname),
                    os.path.join(base_logs_dir, fname),
                )

    # Ricalcola summary.csv e stability_summary.csv
    recalculate_summary(base_dir)


def main():
    parser = argparse.ArgumentParser(
        description="Merge dei risultati di una nuova run di benchmark (es. parfor ottimizzato) nella cartella base.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Esempi di utilizzo:
  # 1. Sovrascrive i dati di ff parfor nella cartella results/ con quelli della nuova run:
  python3 scripts/merge_benchmark.py results/ results_parfor_new/

  # 2. Conserva la vecchia versione rinominandola 'parfor-naive' (per confrontarle nei grafici):
  python3 scripts/merge_benchmark.py results/ results_parfor_new/ --keep-old parfor-naive

  # 3. Solo ricalcolo di summary.csv da timings_raw.csv esistente:
  python3 scripts/merge_benchmark.py results/ --recalculate-only
""",
    )
    parser.add_argument(
        "base_dir",
        help="Directory principale dei risultati consolidati (es. fmo/results/)",
    )
    parser.add_argument(
        "new_dir",
        nargs="?",
        default=None,
        help="Directory contenente i dati della nuova run di benchmark",
    )
    parser.add_argument(
        "--variant",
        default="ff",
        help="Variante da aggiornare (default: 'ff')",
    )
    parser.add_argument(
        "--mode",
        default="parfor",
        help="Modalità da aggiornare (default: 'parfor')",
    )
    parser.add_argument(
        "--keep-old",
        default=None,
        metavar="LABEL",
        help="Se specificato, non cancella il vecchio parfor ma lo rinomina con LABEL (es. 'parfor-naive') per confrontarlo nei grafici.",
    )
    parser.add_argument(
        "--recalculate-only",
        action="store_true",
        help="Ricalcola solo summary.csv e stability_summary.csv per base_dir senza fare merge.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Non rigenerare automaticamente i grafici e il report HTML alla fine.",
    )

    args = parser.parse_args()

    if args.recalculate_only:
        recalculate_summary(args.base_dir)
    else:
        if not args.new_dir:
            print(
                "Errore: specificare sia la directory base sia la nuova directory, oppure usare --recalculate-only.",
                file=sys.stderr,
            )
            sys.exit(1)
        merge_benchmark_runs(
            args.base_dir,
            args.new_dir,
            target_variant=args.variant,
            target_mode=args.mode,
            keep_old_as=args.keep_old,
        )

    if not args.no_plot:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        plot_script = os.path.join(script_dir, "plot_benchmark.py")
        if os.path.exists(plot_script):
            print("\nRigenerazione grafici e report HTML con plot_benchmark.py...")
            subprocess.run([sys.executable, plot_script, args.base_dir])


if __name__ == "__main__":
    main()
