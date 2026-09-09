#!/usr/bin/env python3
"""
plot_benchmark.py - Analisi e visualizzazione delle prestazioni di FMO (NSGA-II)

Genera:
1. Tabelle riassuntive a terminale (Speedup seq, Speedup T1, Efficienza)
2. Grafici vettoriali SVG standalone (speedup.svg, efficiency.svg, execution_time.svg)
   funzionanti al 100% offline senza alcuna dipendenza esterna!
3. Report interattivo HTML (report.html) con grafici dinamici e selettore di fase
4. Se disponibile matplotlib, figure PNG/PDF ad alta risoluzione

Uso:
  python3 plot_benchmark.py logs/
"""

import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict

# Controllo presenza matplotlib (opzionale)
HAS_MATPLOTLIB = False
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    pass


class PlotConfig:
    """Configurazione per nomi dei file, percorsi e titoli dei grafici."""

    def __init__(
        self,
        output_dir=None,
        prefix="",
        suffix="",
        speedup_name="speedup",
        scalability_name="scalability",
        efficiency_name="efficiency",
        time_name="execution_time",
        stability_name="solution_stability",
        report_name="report",
        title_prefix="",
        title_suffix="",
        speedup_title=None,
        scalability_title=None,
        efficiency_title=None,
        time_title=None,
        stability_title=None,
        theme="light",
        amdahl_f=None,
        show_amdahl=True,
    ):
        self.output_dir = output_dir or "."
        self.prefix = prefix
        self.suffix = suffix
        self.speedup_name = speedup_name
        self.scalability_name = scalability_name
        self.efficiency_name = efficiency_name
        self.time_name = time_name
        self.stability_name = stability_name
        self.report_name = report_name
        self.title_prefix = title_prefix
        self.title_suffix = title_suffix
        self.speedup_title = speedup_title
        self.scalability_title = scalability_title
        self.efficiency_title = efficiency_title
        self.time_title = time_title
        self.stability_title = stability_title
        self.theme = theme
        self.amdahl_f = amdahl_f
        self.show_amdahl = show_amdahl

    def file(self, base_name, ext):
        if not ext.startswith("."):
            ext = "." + ext
        return f"{self.prefix}{base_name}{self.suffix}{ext}"

    def path(self, base_name, ext):
        return os.path.join(self.output_dir, self.file(base_name, ext))

    def get_title(self, default_title, override=None):
        if override:
            return override
        res = default_title
        if self.title_prefix:
            res = f"{self.title_prefix}{res}"
        if self.title_suffix:
            res = f"{res}{self.title_suffix}"
        return res


def find_latest_results_dir(base_dir):
    results_dir = os.path.join(base_dir, "results")
    if not os.path.exists(results_dir):
        return None
    runs = [
        os.path.join(results_dir, d)
        for d in os.listdir(results_dir)
        if d.startswith("run_")
    ]
    if not runs:
        return None
    runs.sort(key=os.path.getmtime, reverse=True)
    return runs[0]


def compute_amdahl_fraction(summary_rows, target_phase="nsga2_total"):
    """
    Calcola la frazione non parallelizzabile (f) secondo la Legge di Amdahl dai dati del baseline sequenziale.
    Restituisce None se non ci sono dati sequenziali sufficienti.
    """
    seq_map = {
        r["phase"]: float(r["avg_mean_ms"])
        for r in summary_rows
        if r.get("variant") == "seq"
    }
    if not seq_map:
        return None

    if target_phase == "nsga2_total" and "nsga2_total" in seq_map:
        t_total = seq_map["nsga2_total"]
        t_init_eval = seq_map.get("initial_evaluation", 0.0)
        t_pop_eval = seq_map.get("population_evaluation", 0.0)
        t_par = t_init_eval + 50 * t_pop_eval
        if t_total > 0 and 0 < t_par < t_total:
            return (t_total - t_par) / t_total
    elif target_phase == "generation_total" and "generation_total" in seq_map:
        t_gen = seq_map["generation_total"]
        t_pop_eval = seq_map.get("population_evaluation", 0.0)
        if t_gen > 0 and 0 < t_pop_eval < t_gen:
            return (t_gen - t_pop_eval) / t_gen

    return None


def load_data(target_path):
    """
    Carica i dati da summary.csv oppure calcola il sommario da timings_raw.csv.
    """
    if os.path.isdir(target_path):
        summary_csv = os.path.join(target_path, "summary.csv")
        raw_csv = os.path.join(target_path, "timings_raw.csv")
        results_dir = target_path
    elif os.path.isfile(target_path):
        results_dir = os.path.dirname(target_path)
        if target_path.endswith("summary.csv"):
            summary_csv = target_path
            raw_csv = os.path.join(results_dir, "timings_raw.csv")
        else:
            raw_csv = target_path
            summary_csv = os.path.join(results_dir, "summary.csv")
    else:
        print(f"Errore: Percorso '{target_path}' non valido.", file=sys.stderr)
        sys.exit(1)

    summary_rows = []
    if os.path.exists(summary_csv):
        with open(summary_csv, "r", newline="") as f:
            reader = csv.DictReader(f)
            summary_rows = list(reader)

    raw_rows = []
    if os.path.exists(raw_csv):
        with open(raw_csv, "r", newline="") as f:
            reader = csv.DictReader(f)
            raw_rows = list(reader)

    stab_summary_csv = os.path.join(results_dir, "stability_summary.csv")
    stab_raw_csv = os.path.join(results_dir, "stability_rmse.csv")

    stab_rows = []
    if os.path.exists(stab_summary_csv):
        try:
            with open(stab_summary_csv, "r", newline="") as f:
                stab_rows = list(csv.DictReader(f))
        except Exception:
            pass
    elif os.path.exists(stab_raw_csv):
        try:
            with open(stab_raw_csv, "r", newline="") as f:
                reader = csv.DictReader(f)
                groups = defaultdict(list)
                for r in reader:
                    w_val = int(r["workers"])
                    groups[(r["variant"], r["mode"], w_val)].append(
                        {
                            "rmse_ptv": float(r["rmse_ptv"]),
                            "rmse_rectum": float(r["rmse_rectum"]),
                            "rmse_bladder": float(r["rmse_bladder"]),
                            "rmse_total_fitness": float(r["rmse_total_fitness"]),
                            "rmse_genes": float(r["rmse_genes"]),
                            "spread_par": float(r["spread_par"]),
                        }
                    )
                for (v, m, w_val), items in sorted(groups.items()):
                    c = len(items)
                    avg_fit = sum(x["rmse_total_fitness"] for x in items) / c
                    avg_gen = sum(x["rmse_genes"] for x in items) / c
                    avg_ptv = sum(x["rmse_ptv"] for x in items) / c
                    avg_rec = sum(x["rmse_rectum"] for x in items) / c
                    avg_bla = sum(x["rmse_bladder"] for x in items) / c
                    avg_spread = sum(x["spread_par"] for x in items) / c
                    status = (
                        "BASELINE"
                        if v == "seq"
                        else (
                            "DETERMINISTICO (0 err)"
                            if avg_fit == 0.0 and avg_gen == 0.0
                            else "STABILE (Dev. GA)"
                        )
                    )
                    stab_rows.append(
                        {
                            "variant": v,
                            "mode": m,
                            "workers": str(w_val),
                            "repetitions": str(c),
                            "avg_rmse_ptv": f"{avg_ptv:.6f}",
                            "avg_rmse_rectum": f"{avg_rec:.6f}",
                            "avg_rmse_bladder": f"{avg_bla:.6f}",
                            "avg_rmse_fitness": f"{avg_fit:.6f}",
                            "avg_rmse_genes": f"{avg_gen:.6f}",
                            "avg_spread": f"{avg_spread:.6f}",
                            "status": status,
                        }
                    )
        except Exception:
            pass

    return summary_rows, raw_rows, stab_rows, results_dir


def print_terminal_summary(summary_rows, stab_rows=None, config=None):
    if not summary_rows and not stab_rows:
        print("Nessun dato di riepilogo disponibile.")
        return

    if summary_rows:
        phases = sorted({r["phase"] for r in summary_rows})
        target_phase = "nsga2_total" if "nsga2_total" in phases else phases[0]

        print("\n" + "=" * 107)
        print(f"  RIEPILOGO PRESTAZIONI (Fase: {target_phase})")
        print("=" * 107)
        print(
            f"{'VARIANTE':<14} {'MODALITÀ':<14} {'WORKERS':<9} {'TEMPO MEDIO (s)':<18} {'SPEEDUP (seq)':<16} {'SCALABILITÀ (T1)':<18} {'EFFICIENZA':<14}"
        )
        print("-" * 107)

        def sort_key(r):
            return (r["variant"], r["mode"], int(r["workers"]))

        phase_rows = [r for r in summary_rows if r["phase"] == target_phase]
        phase_rows.sort(key=sort_key)

        for r in phase_rows:
            variant = r["variant"]
            mode = r["mode"]
            workers = r["workers"]
            mean_s = float(r["avg_mean_ms"]) / 1000.0

            s_seq = f"{float(r['speedup_seq']):.2f}x" if r.get("speedup_seq") else "-"
            s_t1 = f"{float(r['speedup_t1']):.2f}x" if r.get("speedup_t1") else "-"
            eff = (
                f"{float(r['efficiency_t1']) * 100:.1f}%"
                if r.get("efficiency_t1")
                else "-"
            )

            print(
                f"{variant:<14} {mode:<14} {workers:<9} {mean_s:<18.3f} {s_seq:<16} {s_t1:<18} {eff:<14}"
            )

        amdahl_f = None
        if config is not None and getattr(config, "show_amdahl", True):
            if getattr(config, "amdahl_f", None) is not None:
                amdahl_f = config.amdahl_f
            else:
                amdahl_f = compute_amdahl_fraction(summary_rows, target_phase)
        elif config is None:
            amdahl_f = compute_amdahl_fraction(summary_rows, target_phase)

        if amdahl_f is not None and 0.0 < amdahl_f < 1.0:
            print("-" * 107)
            s_max_amd = 1.0 / amdahl_f
            print(
                f"  [Legge di Amdahl] Frazione seriale (f): {amdahl_f * 100:.2f}% | Limite Speedup teorico asintotico (1/f): {s_max_amd:.2f}x"
            )

        print("=" * 107 + "\n")

    if stab_rows:
        print("=" * 110)
        print(
            "  STABILITÀ DELLE SOLUZIONI (RMSE vs Baseline Sequenziale) & DETERMINISMO"
        )
        print("=" * 110)
        print(
            f"{'VARIANTE':<10} {'MODE':<10} {'WORKERS':<9} {'RMSE_FITNESS':<14} {'RMSE_GENI':<14} {'SPREAD_PAR':<12} {'STATO DETERMINISMO':<30}"
        )
        print("-" * 110)
        for r in stab_rows:
            fit_val = float(r.get("avg_rmse_fitness", 0.0))
            gen_val = float(r.get("avg_rmse_genes", 0.0))
            spread_val = float(r.get("avg_spread", 0.0))
            print(
                f"{r['variant']:<10} {r['mode']:<10} {r['workers']:<9} {fit_val:<14.6f} {gen_val:<14.6f} {spread_val:<12.4f} {r.get('status', '-'):<30}"
            )
        print("=" * 110 + "\n")


def generate_svg_chart(
    title,
    x_label,
    y_label,
    series_data,
    out_path,
    is_speedup=False,
    is_efficiency=False,
    is_time=False,
    is_stability=False,
    speedup_mode="log",
    time_mode="linear",
    theme="light",
    amdahl_f=None,
):
    """
    Genera un grafico vettoriale SVG professionale e pulito, pronto per report e tesi.
    Supporta tema chiaro (sfondo bianco puro, Ideal per report/stampe) e scuro,
    spaziatura uniforme log2 per le potenze di 2 e marcatori distinti per ogni serie.
    """
    w, h = 800, 480
    pad_left, pad_right, pad_top, pad_bottom = 75, 185, 55, 60
    plot_w = w - pad_left - pad_right
    plot_h = h - pad_top - pad_bottom

    # Palette ad alto contrasto adatta sia a schermo che a stampa
    palette = [
        {"color": "#0969da", "shape": "circle"},  # Blu
        {"color": "#1a7f37", "shape": "rect"},  # Verde
        {"color": "#d1242f", "shape": "triangle"},  # Rosso
        {"color": "#8250df", "shape": "diamond"},  # Viola
        {"color": "#e36209", "shape": "circle"},  # Arancione
        {"color": "#0891b2", "shape": "rect"},  # Ciano
    ]

    if theme == "light":
        bg_canvas = "#ffffff"
        bg_plot = "#ffffff"
        border_plot = "#334155"
        grid_major = "#e2e8f0"
        grid_minor = "#f8fafc"
        text_title = "#0f172a"
        text_axis = "#1e293b"
        text_tick = "#475569"
        card_bg = "#f8fafc"
        card_border = "#cbd5e1"
        legend_text = "#1e293b"
        ideal_color = "#64748b"
    else:
        bg_canvas = "#0d1117"
        bg_plot = "#161b22"
        border_plot = "#30363d"
        grid_major = "#21262d"
        grid_minor = "#161b22"
        text_title = "#f0f6fc"
        text_axis = "#8b949e"
        text_tick = "#8b949e"
        card_bg = "#161b22"
        card_border = "#30363d"
        legend_text = "#c9d1d9"
        ideal_color = "#8b949e"

    all_x = []
    all_y = []
    for label, points in series_data.items():
        for px, py in points:
            if px is not None and py is not None:
                all_x.append(px)
                all_y.append(py)

    if not all_x or not all_y:
        return

    unique_x = sorted(set(all_x))
    min_x, max_x = min(unique_x), max(unique_x)

    # Rilevamento scala log2 su X (evita sovrapposizioni su 1, 2, 4...)
    use_log_x = min_x > 0 and (max_x / min_x >= 4) and all(x > 0 for x in unique_x)
    log2_min_x = math.log2(min_x) if use_log_x else 0
    log2_max_x = math.log2(max_x) if use_log_x else 0

    def to_px(val_x):
        if use_log_x:
            if log2_max_x == log2_min_x:
                return pad_left + plot_w / 2
            return (
                pad_left
                + (math.log2(val_x) - log2_min_x) / (log2_max_x - log2_min_x) * plot_w
            )
        else:
            if max_x == min_x:
                return pad_left + plot_w / 2
            return pad_left + (val_x - min_x) / (max_x - min_x) * plot_w

    use_log_y = False
    if is_speedup and speedup_mode == "log":
        use_log_y = True
        min_y = 1.0
        max_y = max(2.0, float(max_x))
        log2_min_y = 0.0
        log2_max_y = math.log2(max_y)
    elif is_speedup and speedup_mode == "linear":
        min_y = 0.0
        max_actual_y = max(all_y)
        max_y = math.ceil(max_actual_y * 1.15 / 5.0) * 5.0
        max_y = max(max_y, 2.0)
    elif is_efficiency:
        min_y = 0.0
        max_y = 105.0
    elif is_time and time_mode == "log":
        use_log_y = True
        min_y = max(1.0, 10 ** math.floor(math.log10(min(all_y))))
        max_y = 10 ** math.ceil(math.log10(max(all_y)))
        if max_y <= min_y:
            max_y = min_y * 10.0
        log10_min_y = math.log10(min_y)
        log10_max_y = math.log10(max_y)
    elif is_time:
        min_y = 0.0
        max_y = max(all_y) * 1.08
    elif is_stability:
        min_y = 0.0
        max_act = max(all_y)
        max_y = math.ceil(max_act * 1.25 * 10.0) / 10.0 if max_act > 0 else 1.0
    else:
        min_y = 0.0
        max_y = max(all_y) * 1.15 if max(all_y) > 0 else 1.0

    def to_py(val_y):
        if use_log_y and is_speedup:
            safe_y = max(val_y, min_y)
            denom = log2_max_y - log2_min_y
            return (
                pad_top
                + plot_h
                - (
                    (math.log2(safe_y) - log2_min_y) / denom * plot_h
                    if denom != 0
                    else 0.0
                )
            )
        elif use_log_y and is_time:
            safe_y = max(val_y, min_y)
            denom = log10_max_y - log10_min_y
            return (
                pad_top
                + plot_h
                - (
                    (math.log10(safe_y) - log10_min_y) / denom * plot_h
                    if denom != 0
                    else 0.0
                )
            )
        else:
            if max_y == min_y:
                return pad_top + plot_h / 2
            return pad_top + plot_h - (val_y - min_y) / (max_y - min_y) * plot_h

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">',
        f'<rect width="{w}" height="{h}" fill="{bg_canvas}" />',
        f'<text x="{pad_left}" y="36" fill="{text_title}" font-family="system-ui, -apple-system, sans-serif" font-size="16" font-weight="bold">{title}</text>',
        f'<rect x="{pad_left}" y="{pad_top}" width="{plot_w}" height="{plot_h}" fill="{bg_plot}" stroke="{border_plot}" stroke-width="1.2" />',
    ]

    # Asse Y
    if use_log_y and is_speedup:
        powers = range(int(log2_min_y), int(log2_max_y) + 1)
        for p in powers:
            val_y = 2**p
            py = to_py(val_y)
            svg.append(
                f'<line x1="{pad_left}" y1="{py:.1f}" x2="{pad_left + plot_w}" y2="{py:.1f}" stroke="{grid_major}" stroke-width="1" stroke-dasharray="3,3" />'
            )
            svg.append(
                f'<text x="{pad_left - 10}" y="{py + 4:.1f}" fill="{text_tick}" font-family="system-ui, -apple-system, sans-serif" font-size="11" text-anchor="end">{int(val_y)}x</text>'
            )
    elif use_log_y and is_time:
        for p in range(int(log10_min_y), int(log10_max_y) + 1):
            val_y = 10**p
            py = to_py(val_y)
            svg.append(
                f'<line x1="{pad_left}" y1="{py:.1f}" x2="{pad_left + plot_w}" y2="{py:.1f}" stroke="{grid_major}" stroke-width="1" stroke-dasharray="3,3" />'
            )
            svg.append(
                f'<text x="{pad_left - 10}" y="{py + 4:.1f}" fill="{text_tick}" font-family="system-ui, -apple-system, sans-serif" font-size="11" text-anchor="end">{int(val_y)}s</text>'
            )
    else:
        y_ticks = 5
        for i in range(y_ticks + 1):
            y_val = min_y + (max_y - min_y) * (i / y_ticks)
            py = to_py(y_val)
            svg.append(
                f'<line x1="{pad_left}" y1="{py:.1f}" x2="{pad_left + plot_w}" y2="{py:.1f}" stroke="{grid_major}" stroke-width="1" stroke-dasharray="3,3" />'
            )
            if is_efficiency:
                y_text = f"{y_val:.0f}%"
            elif is_speedup:
                y_text = f"{y_val:.1f}x"
            elif is_time:
                y_text = f"{y_val:.0f}s" if y_val >= 10 else f"{y_val:.2f}s"
            elif is_stability:
                y_text = f"{y_val:.2f}"
            else:
                y_text = f"{y_val:.1f}"
            svg.append(
                f'<text x="{pad_left - 10}" y="{py + 4:.1f}" fill="{text_tick}" font-family="system-ui, -apple-system, sans-serif" font-size="11" text-anchor="end">{y_text}</text>'
            )

    # Asse X
    for ux in unique_x:
        px = to_px(ux)
        svg.append(
            f'<line x1="{px:.1f}" y1="{pad_top}" x2="{px:.1f}" y2="{pad_top + plot_h}" stroke="{grid_minor}" stroke-width="1" stroke-dasharray="2,2" />'
        )
        svg.append(
            f'<line x1="{px:.1f}" y1="{pad_top + plot_h}" x2="{px:.1f}" y2="{pad_top + plot_h + 5}" stroke="{border_plot}" stroke-width="1.2" />'
        )
        svg.append(
            f'<text x="{px:.1f}" y="{pad_top + plot_h + 18}" fill="{text_tick}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500" text-anchor="middle">{ux}</text>'
        )

    # Linea ideale (Speedup o Efficienza)
    if is_speedup:
        if use_log_y:
            p1_x, p1_y = to_px(min_x), to_py(min_x)
            p2_x, p2_y = to_px(max_x), to_py(max_x)
            svg.append(
                f'<line x1="{p1_x:.1f}" y1="{p1_y:.1f}" x2="{p2_x:.1f}" y2="{p2_y:.1f}" stroke="{ideal_color}" stroke-width="1.6" stroke-dasharray="5,4" />'
            )
        else:
            ideal_pts = []
            for ux in unique_x:
                if ux <= max_y:
                    ideal_pts.append((to_px(ux), to_py(ux)))
                else:
                    ideal_pts.append((to_px(ux), to_py(max_y)))
                    break
            if len(ideal_pts) > 1:
                poly = " ".join([f"{x:.1f},{y:.1f}" for x, y in ideal_pts])
                svg.append(
                    f'<polyline fill="none" stroke="{ideal_color}" stroke-width="1.6" stroke-dasharray="5,4" points="{poly}" />'
                )

        if amdahl_f is not None and 0.0 < amdahl_f < 1.0:
            amdahl_color = "#e36209" if theme == "light" else "#f0883e"
            step_count = 60
            if use_log_x and min_x > 0:
                xs = [
                    2 ** (log2_min_x + i * (log2_max_x - log2_min_x) / step_count)
                    for i in range(step_count + 1)
                ]
            else:
                xs = [
                    min_x + i * (max_x - min_x) / step_count
                    for i in range(step_count + 1)
                ]

            amdahl_pts = []
            for cur_x in xs:
                if cur_x <= 0:
                    continue
                s_amd = 1.0 / (amdahl_f + (1.0 - amdahl_f) / cur_x)
                if not use_log_y:
                    s_amd_clamped = min(s_amd, max_y)
                    amdahl_pts.append((to_px(cur_x), to_py(s_amd_clamped)))
                else:
                    if s_amd >= min_y:
                        amdahl_pts.append((to_px(cur_x), to_py(min(s_amd, max_y))))
            if len(amdahl_pts) > 1:
                poly_amd = " ".join([f"{x:.1f},{y:.1f}" for x, y in amdahl_pts])
                svg.append(
                    f'<polyline fill="none" stroke="{amdahl_color}" stroke-width="1.8" stroke-dasharray="4,3" points="{poly_amd}" />'
                )
    elif is_efficiency:
        py_100 = to_py(100.0)
        svg.append(
            f'<line x1="{pad_left}" y1="{py_100:.1f}" x2="{pad_left + plot_w}" y2="{py_100:.1f}" stroke="{ideal_color}" stroke-width="1.6" stroke-dasharray="5,4" />'
        )

    # Tracciamento serie dati
    legend_entries = []
    for idx, (label, points) in enumerate(series_data.items()):
        style = palette[idx % len(palette)]
        color = style["color"]
        shape = style["shape"]
        valid_pts = [(px, py) for px, py in points if px is not None and py is not None]
        valid_pts.sort(key=lambda p: p[0])

        if len(valid_pts) > 1:
            poly_points = " ".join(
                [f"{to_px(px):.1f},{to_py(py):.1f}" for px, py in valid_pts]
            )
            svg.append(
                f'<polyline fill="none" stroke="{color}" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" points="{poly_points}" />'
            )

        for px, py in valid_pts:
            cx, cy = to_px(px), to_py(py)
            if shape == "rect":
                svg.append(
                    f'<rect x="{cx - 4:.1f}" y="{cy - 4:.1f}" width="8" height="8" fill="{color}" stroke="#ffffff" stroke-width="1.5" />'
                )
            elif shape == "triangle":
                svg.append(
                    f'<polygon points="{cx:.1f},{cy - 5.5:.1f} {cx + 5:.1f},{cy + 4.5:.1f} {cx - 5:.1f},{cy + 4.5:.1f}" fill="{color}" stroke="#ffffff" stroke-width="1.5" />'
                )
            elif shape == "diamond":
                svg.append(
                    f'<polygon points="{cx:.1f},{cy - 5.5:.1f} {cx + 5:.1f},{cy:.1f} {cx:.1f},{cy + 5.5:.1f} {cx - 5:.1f},{cy:.1f}" fill="{color}" stroke="#ffffff" stroke-width="1.5" />'
                )
            else:
                svg.append(
                    f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.5" fill="{color}" stroke="#ffffff" stroke-width="1.5" />'
                )

        legend_entries.append((label, color, shape))

    # Box Legenda
    has_amdahl = is_speedup and (amdahl_f is not None and 0.0 < amdahl_f < 1.0)
    extra_h = 0
    if is_speedup:
        extra_h = 46 if has_amdahl else 24
    elif is_efficiency:
        extra_h = 24
    else:
        extra_h = 10
    legend_h = len(legend_entries) * 22 + extra_h + 12
    legend_x = pad_left + plot_w + 14
    legend_y = pad_top
    svg.append(
        f'<rect x="{legend_x}" y="{legend_y}" width="{pad_right - 24}" height="{legend_h}" fill="{card_bg}" stroke="{card_border}" stroke-width="1" rx="6" />'
    )

    cur_y = legend_y + 18
    for label, color, shape in legend_entries:
        icon_cx = legend_x + 16
        icon_cy = cur_y - 4
        svg.append(
            f'<line x1="{icon_cx - 8}" y1="{icon_cy}" x2="{icon_cx + 8}" y2="{icon_cy}" stroke="{color}" stroke-width="2" />'
        )
        if shape == "rect":
            svg.append(
                f'<rect x="{icon_cx - 3.5}" y="{icon_cy - 3.5}" width="7" height="7" fill="{color}" stroke="#ffffff" stroke-width="1" />'
            )
        elif shape == "triangle":
            svg.append(
                f'<polygon points="{icon_cx},{icon_cy - 4} {icon_cx + 4},{icon_cy + 3.5} {icon_cx - 4},{icon_cy + 3.5}" fill="{color}" stroke="#ffffff" stroke-width="1" />'
            )
        elif shape == "diamond":
            svg.append(
                f'<polygon points="{icon_cx},{icon_cy - 4.5} {icon_cx + 4},{icon_cy} {icon_cx},{icon_cy + 4.5} {icon_cx - 4},{icon_cy}" fill="{color}" stroke="#ffffff" stroke-width="1" />'
            )
        else:
            svg.append(
                f'<circle cx="{icon_cx}" cy="{icon_cy}" r="3.5" fill="{color}" stroke="#ffffff" stroke-width="1" />'
            )

        svg.append(
            f'<text x="{legend_x + 32}" y="{cur_y}" fill="{legend_text}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500">{label}</text>'
        )
        cur_y += 22

    if is_speedup:
        icon_cx = legend_x + 16
        icon_cy = cur_y - 4
        svg.append(
            f'<line x1="{icon_cx - 8}" y1="{icon_cy}" x2="{icon_cx + 8}" y2="{icon_cy}" stroke="{ideal_color}" stroke-dasharray="3,2" stroke-width="1.8" />'
        )
        svg.append(
            f'<text x="{legend_x + 32}" y="{cur_y}" fill="{ideal_color}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500">Ideal (S = p)</text>'
        )
        cur_y += 22
        if has_amdahl:
            amdahl_color = "#e36209" if theme == "light" else "#f0883e"
            icon_cy = cur_y - 4
            svg.append(
                f'<line x1="{icon_cx - 8}" y1="{icon_cy}" x2="{icon_cx + 8}" y2="{icon_cy}" stroke="{amdahl_color}" stroke-dasharray="4,3" stroke-width="1.8" />'
            )
            pct_label = f"{amdahl_f * 100:.2f}%"
            svg.append(
                f'<text x="{legend_x + 32}" y="{cur_y}" fill="{amdahl_color}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500">Amdahl (f={pct_label})</text>'
            )
            cur_y += 22
    elif is_efficiency:
        icon_cx = legend_x + 16
        icon_cy = cur_y - 4
        svg.append(
            f'<line x1="{icon_cx - 8}" y1="{icon_cy}" x2="{icon_cx + 8}" y2="{icon_cy}" stroke="{ideal_color}" stroke-dasharray="3,2" stroke-width="1.8" />'
        )
        svg.append(
            f'<text x="{legend_x + 32}" y="{cur_y}" fill="{ideal_color}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500">Ideal (100%)</text>'
        )

    # Etichette assi
    svg.append(
        f'<text x="{pad_left + plot_w / 2}" y="{h - 15}" fill="{text_axis}" font-family="system-ui, -apple-system, sans-serif" font-size="12" font-weight="600" text-anchor="middle">{x_label}</text>'
    )
    svg.append(
        f'<text x="22" y="{pad_top + plot_h / 2}" fill="{text_axis}" font-family="system-ui, -apple-system, sans-serif" font-size="12" font-weight="600" text-anchor="middle" transform="rotate(-90 22 {pad_top + plot_h / 2})">{y_label}</text>'
    )

    svg.append("</svg>")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(svg))

    print(f"Grafico SVG salvato: {out_path}")


def generate_all_svg_charts(summary_rows, results_dir, stab_rows=None, theme="light", config=None):
    """
    Genera i grafici SVG standalone per la fase principale e per la stabilità.
    """
    if config is None:
        config = PlotConfig(output_dir=results_dir, theme=theme)

    target_dir = config.output_dir
    os.makedirs(target_dir, exist_ok=True)
    theme = config.theme

    phases = sorted({r["phase"] for r in summary_rows})
    target_phase = "nsga2_total" if "nsga2_total" in phases else phases[0]
    rows = [r for r in summary_rows if r["phase"] == target_phase]

    speedup_series = defaultdict(list)
    scalability_series = defaultdict(list)
    eff_series = defaultdict(list)
    time_series = defaultdict(list)

    for r in rows:
        label = f"{r['variant']} ({r['mode']})"
        w = int(r["workers"])
        time_s = float(r["avg_mean_ms"]) / 1000.0
        s_seq = float(r["speedup_seq"]) if r.get("speedup_seq") else None
        s_t1 = float(r["speedup_t1"]) if r.get("speedup_t1") else None
        eff_t1 = float(r["efficiency_t1"]) * 100.0 if r.get("efficiency_t1") else None

        time_series[label].append((w, time_s))
        if s_seq is not None:
            speedup_series[label].append((w, s_seq))
        if s_t1 is not None:
            scalability_series[label].append((w, s_t1))
        if eff_t1 is not None:
            eff_series[label].append((w, eff_t1))

    amdahl_f = None
    if getattr(config, "show_amdahl", True):
        if getattr(config, "amdahl_f", None) is not None:
            amdahl_f = config.amdahl_f
        else:
            amdahl_f = compute_amdahl_fraction(summary_rows, target_phase)

    # 1. Speedup Assoluto (vs baseline sequenziale pura T_seq)
    speedup_title = config.get_title(f"Speedup Assoluto vs Cores ({target_phase})", config.speedup_title)
    generate_svg_chart(
        speedup_title,
        "Cores / Workers",
        "Speedup Assoluto (Tseq / Tp)",
        speedup_series,
        config.path(config.speedup_name, "svg"),
        is_speedup=True,
        speedup_mode="log",
        theme=theme,
        amdahl_f=amdahl_f,
    )

    speedup_linear_title = config.get_title(
        f"Speedup Assoluto vs Cores - Scala Lineare ({target_phase})",
        f"{config.speedup_title} (Scala Lineare)" if config.speedup_title else None,
    )
    generate_svg_chart(
        speedup_linear_title,
        "Cores / Workers",
        "Speedup Assoluto (Tseq / Tp)",
        speedup_series,
        config.path(f"{config.speedup_name}_linear", "svg"),
        is_speedup=True,
        speedup_mode="linear",
        theme=theme,
        amdahl_f=amdahl_f,
    )

    # 2. Scalabilità Forte (Relativa a 1 worker T1: scalab(n) = T1 / Tn)
    scalability_title = config.get_title(f"Scalabilità vs Cores ({target_phase})", config.scalability_title)
    generate_svg_chart(
        scalability_title,
        "Cores / Workers",
        "Scalabilità (T1 / Tp)",
        scalability_series,
        config.path(config.scalability_name, "svg"),
        is_speedup=True,
        speedup_mode="log",
        theme=theme,
        amdahl_f=amdahl_f,
    )

    scalability_linear_title = config.get_title(
        f"Scalabilità vs Cores - Scala Lineare ({target_phase})",
        f"{config.scalability_title} (Scala Lineare)" if config.scalability_title else None,
    )
    generate_svg_chart(
        scalability_linear_title,
        "Cores / Workers",
        "Scalabilità (T1 / Tp)",
        scalability_series,
        config.path(f"{config.scalability_name}_linear", "svg"),
        is_speedup=True,
        speedup_mode="linear",
        theme=theme,
        amdahl_f=amdahl_f,
    )

    # 3. Efficienza Parallela
    efficiency_title = config.get_title(f"Efficiency ({target_phase})", config.efficiency_title)
    generate_svg_chart(
        efficiency_title,
        "Cores / Workers",
        "Efficiency (%)",
        eff_series,
        config.path(config.efficiency_name, "svg"),
        is_efficiency=True,
        theme=theme,
    )

    # 4. Completion time lineare
    time_title = config.get_title(f"Completion time vs Cores ({target_phase})", config.time_title)
    generate_svg_chart(
        time_title,
        "Cores / Workers",
        "Tempo (secondi)",
        time_series,
        config.path(config.time_name, "svg"),
        is_time=True,
        time_mode="linear",
        theme=theme,
    )

    # 5. Completion time logaritmico
    time_log_title = config.get_title(
        f"Completion time vs Cores - Log scale ({target_phase})",
        f"{config.time_title} (Log scale)" if config.time_title else None,
    )
    generate_svg_chart(
        time_log_title,
        "Cores / Workers",
        "log_10(seconds)",
        time_series,
        config.path(f"{config.time_name}_log", "svg"),
        is_time=True,
        time_mode="log",
        theme=theme,
    )

    # 6. Stabilità delle soluzioni
    if stab_rows:
        stab_series = defaultdict(list)
        for r in stab_rows:
            if r["variant"] == "seq":
                continue
            label = f"{r['variant']} ({r['mode']})"
            w = int(r["workers"])
            fit_rmse = float(r.get("avg_rmse_fitness", 0.0))
            stab_series[label].append((w, fit_rmse))

        if stab_series:
            stability_title = config.get_title(
                "Solution stability: RMSE Fitness vs Cores", config.stability_title
            )
            generate_svg_chart(
                stability_title,
                "Cores / Workers",
                "RMSE Fitness (vs Baseline)",
                stab_series,
                config.path(config.stability_name, "svg"),
                is_stability=True,
                theme=theme,
            )



def generate_html_report(summary_rows, raw_rows, results_dir, stab_rows=None, config=None):
    """
    Genera un report HTML completo e interattivo con grafici Chart.js (online)
    e link diretti agli SVG vettoriali (offline).
    """
    if config is None:
        config = PlotConfig(output_dir=results_dir)

    target_dir = config.output_dir
    os.makedirs(target_dir, exist_ok=True)
    html_path = config.path(config.report_name, "html")

    phase_data = defaultdict(lambda: defaultdict(dict))
    for r in summary_rows:
        phase = r["phase"]
        label = f"{r['variant']} ({r['mode']})"
        w = int(r["workers"])
        time_s = float(r["avg_mean_ms"]) / 1000.0
        s_seq = float(r["speedup_seq"]) if r.get("speedup_seq") else None
        s_t1 = float(r["speedup_t1"]) if r.get("speedup_t1") else None
        eff_t1 = float(r["efficiency_t1"]) * 100.0 if r.get("efficiency_t1") else None

        phase_data[phase][label][w] = {
            "time_s": time_s,
            "speedup_seq": s_seq,
            "speedup_t1": s_t1,
            "efficiency_t1": eff_t1,
        }

    sys_info = ""
    sys_info_path = os.path.join(results_dir, "system_info.txt")
    if os.path.exists(sys_info_path):
        with open(sys_info_path, "r") as f:
            sys_info = f.read()

    chart_payload = {}
    amdahl_payload = {}
    for phase, series_dict in phase_data.items():
        chart_payload[phase] = {}
        f_val = compute_amdahl_fraction(summary_rows, phase)
        if f_val is not None and 0.0 < f_val < 1.0:
            amdahl_payload[phase] = round(f_val, 6)
        for label, w_dict in series_dict.items():
            sorted_w = sorted(w_dict.keys())
            chart_payload[phase][label] = {
                "workers": sorted_w,
                "times": [w_dict[w]["time_s"] for w in sorted_w],
                "speedup_seq": [w_dict[w]["speedup_seq"] for w in sorted_w],
                "speedup_t1": [w_dict[w]["speedup_t1"] for w in sorted_w],
                "efficiency_t1": [w_dict[w]["efficiency_t1"] for w in sorted_w],
            }

    stab_card_html = ""
    stab_link_html = ""
    if stab_rows:
        stab_svg_file = config.file(config.stability_name, "svg")
        stab_link_html = f"""
            <a href="{stab_svg_file}" target="_blank">🎯 {stab_svg_file}</a>
            <a href="stability_summary.csv" download>📄 stability_summary.csv</a>
            <a href="stability_rmse.csv" download>📄 stability_rmse.csv</a>
        """
        rows_html = []
        for r in stab_rows:
            st = r.get("status", "")
            if "DETERMINISTICO" in st:
                badge = f'<span class="badge" style="background: rgba(63, 185, 80, 0.2); color: #3fb950;">{st}</span>'
            elif "BASELINE" in st:
                badge = f'<span class="badge" style="background: rgba(88, 166, 255, 0.2); color: #58a6ff;">{st}</span>'
            else:
                badge = f'<span class="badge" style="background: rgba(210, 153, 34, 0.2); color: #d29922;">{st}</span>'

            rows_html.append(f"""
                <tr>
                    <td><strong>{r.get("variant", "")}</strong></td>
                    <td>{r.get("mode", "")}</td>
                    <td>{r.get("workers", "1")}</td>
                    <td>{r.get("repetitions", "1")}</td>
                    <td><code>{r.get("avg_rmse_fitness", "-")}</code></td>
                    <td><code>{r.get("avg_rmse_genes", "-")}</code></td>
                    <td><code>{r.get("avg_rmse_ptv", "-")}</code></td>
                    <td><code>{r.get("avg_rmse_rectum", "-")}</code></td>
                    <td><code>{r.get("avg_rmse_bladder", "-")}</code></td>
                    <td><code>{r.get("avg_spread", "-")}</code></td>
                    <td>{badge}</td>
                </tr>
            """)
        stab_card_html = f"""
        <div class="card" style="margin-bottom: 32px;">
            <h2>🎯 Stabilità delle Soluzioni & Determinismo (RMSE vs Baseline Sequenziale)</h2>
            <p style="color: var(--text-muted); font-size: 14px; margin-top: 4px; margin-bottom: 16px;">
                Misurazione della discrepanza tra la popolazione prodotta dalla baseline sequenziale e quella parallela
                tramite <code>computePopulationRMSE</code> e <code>computeSpreadMetric</code>.
                Un valore di RMSE pari a 0 attesta il perfetto determinismo e l'assenza di race condition.
            </p>
            <div style="overflow-x: auto;">
                <table>
                    <thead>
                        <tr>
                            <th>Variante</th>
                            <th>Modo</th>
                            <th>Workers</th>
                            <th>Rip.</th>
                            <th>RMSE Total Fitness</th>
                            <th>RMSE Geni</th>
                            <th>RMSE PTV</th>
                            <th>RMSE Retto</th>
                            <th>RMSE Vescica</th>
                            <th>Spread Δ</th>
                            <th>Stato Determinismo</th>
                        </tr>
                    </thead>
                    <tbody>
                        {"".join(rows_html)}
                    </tbody>
                </table>
            </div>
        </div>
        """

    speedup_svg = config.file(config.speedup_name, "svg")
    speedup_lin_svg = config.file(f"{config.speedup_name}_linear", "svg")
    scalability_svg = config.file(config.scalability_name, "svg")
    scalability_lin_svg = config.file(f"{config.scalability_name}_linear", "svg")
    eff_svg = config.file(config.efficiency_name, "svg")
    time_svg = config.file(config.time_name, "svg")
    time_log_svg = config.file(f"{config.time_name}_log", "svg")

    html_content = f"""<!DOCTYPE html>
<html lang="it">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>FMO NSGA-II HPC Benchmark Report</title>
    <!-- Chart.js (CDN opzionale per grafici dinamici interattivi) -->
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
    <style>
        :root {{
            --bg-color: #0d1117;
            --card-bg: #161b22;
            --border-color: #30363d;
            --text-color: #c9d1d9;
            --text-muted: #8b949e;
            --primary: #58a6ff;
            --success: #3fb950;
            --warning: #d29922;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: var(--bg-color);
            color: var(--text-color);
            margin: 0;
            padding: 24px;
            line-height: 1.5;
        }}
        .container {{
            max-width: 1240px;
            margin: 0 auto;
        }}
        header {{
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 16px;
            margin-bottom: 24px;
        }}
        h1 {{
            color: #f0f6fc;
            margin: 0 0 8px 0;
            font-size: 28px;
        }}
        .links-bar {{
            margin-bottom: 20px;
            padding: 12px 16px;
            background: #21262d;
            border-radius: 6px;
            font-size: 14px;
        }}
        .links-bar a {{
            color: var(--primary);
            margin-right: 16px;
            text-decoration: none;
            font-weight: 600;
        }}
        .links-bar a:hover {{ text-decoration: underline; }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(560px, 1fr));
            gap: 24px;
            margin-bottom: 32px;
        }}
        .card {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.3);
        }}
        .card h2 {{
            margin-top: 0;
            color: #f0f6fc;
            font-size: 18px;
            border-bottom: 1px solid var(--border-color);
            padding-bottom: 8px;
        }}
        .chart-container {{
            position: relative;
            height: 360px;
            width: 100%;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 12px;
            font-size: 14px;
        }}
        th, td {{
            padding: 10px 12px;
            text-align: left;
            border-bottom: 1px solid var(--border-color);
        }}
        th {{
            background-color: #21262d;
            color: #f0f6fc;
            font-weight: 600;
        }}
        tr:hover td {{
            background-color: #21262d;
        }}
        pre {{
            background: #090d13;
            border: 1px solid var(--border-color);
            padding: 12px;
            border-radius: 6px;
            overflow-x: auto;
            font-size: 12px;
            color: #7ee787;
        }}
        .badge {{
            display: inline-block;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 12px;
            font-weight: 600;
        }}
        .badge-primary {{ background: rgba(56, 139, 253, 0.15); color: var(--primary); }}
        select {{
            background: #21262d;
            color: #f0f6fc;
            border: 1px solid var(--border-color);
            padding: 6px 12px;
            border-radius: 6px;
            font-size: 14px;
            margin-bottom: 16px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>⚡ FMO NSGA-II Scalability & Benchmark Report</h1>
            <p style="color: var(--text-muted); margin: 0;">
                Analisi di Speedup, Efficienza Parallela, Scalabilità e Stabilità delle Soluzioni.
            </p>
        </header>

        <div class="links-bar">
            <strong>Grafici Vettoriali e Dati:</strong>
            <a href="{speedup_svg}" target="_blank">🚀 {speedup_svg} (log)</a>
            <a href="{speedup_lin_svg}" target="_blank">🚀 {speedup_lin_svg}</a>
            <a href="{scalability_svg}" target="_blank">📈 {scalability_svg} (log)</a>
            <a href="{scalability_lin_svg}" target="_blank">📈 {scalability_lin_svg}</a>
            <a href="{eff_svg}" target="_blank">⚡ {eff_svg}</a>
            <a href="{time_svg}" target="_blank">⏱️ {time_svg}</a>
            <a href="{time_log_svg}" target="_blank">⏱️ {time_log_svg}</a>
            {stab_link_html}
            <a href="summary.csv" download>📄 summary.csv</a>
            <a href="timings_raw.csv" download>📄 timings_raw.csv</a>
        </div>

        <label for="phaseSelect" style="font-weight: 600;">Seleziona Fase da Visualizzare: </label>
        <select id="phaseSelect" onchange="updateCharts()"></select>

        <div class="grid">
            <div class="card">
                <h2>📈 Scalabilità vs Cores (Scalab = T1 / Tp)</h2>
                <div class="chart-container">
                    <canvas id="scalabilityChart"></canvas>
                </div>
            </div>
            <div class="card">
                <h2>🚀 Speedup Assoluto vs Cores (S = Tseq / Tp)</h2>
                <div class="chart-container">
                    <canvas id="speedupChart"></canvas>
                </div>
            </div>
            <div class="card">
                <h2>⚡ Efficienza Parallela vs Cores (E = Scalab / p)</h2>
                <div class="chart-container">
                    <canvas id="efficiencyChart"></canvas>
                </div>
            </div>
            <div class="card">
                <h2>⏱️ Completion time vs Cores</h2>
                <div class="chart-container">
                    <canvas id="timeChart"></canvas>
                </div>
            </div>
        </div>

        <div class="card" style="margin-bottom: 32px;">
            <h2>📊 Tabella di Riepilogo Prestazioni</h2>
            <div style="overflow-x: auto; max-height: 380px;">
                <table id="summaryTable">
                    <thead>
                        <tr>
                            <th>Variante</th>
                            <th>Modo</th>
                            <th>Cores</th>
                            <th>Tempo (s)</th>
                            <th>Speedup Ass. (Tseq/Tp)</th>
                            <th>Scalabilità (T1/Tp)</th>
                            <th>Efficienza</th>
                        </tr>
                    </thead>
                    <tbody id="tableBody"></tbody>
                </table>
            </div>
        </div>

        {stab_card_html}

        <div class="card" style="margin-bottom: 32px;">
            <h2>🖥️ Metadati Hardware & Ambiente di Esecuzione</h2>
            <pre><code>{sys_info if sys_info else "Nessun metadato hardware rilevato in system_info.txt"}</code></pre>
        </div>
    </div>

    <script>
        const chartData = {json.dumps(chart_payload)};
        const amdahlData = {json.dumps(amdahl_payload)};
        const phases = Object.keys(chartData);

        const selectElem = document.getElementById('phaseSelect');
        phases.forEach(ph => {{
            const opt = document.createElement('option');
            opt.value = ph;
            opt.textContent = ph;
            if (ph === 'nsga2_total') opt.selected = true;
            selectElem.appendChild(opt);
        }});

        let scalabilityChart, speedupChart, efficiencyChart, timeChart;
        const colors = [
            '#58a6ff', '#3fb950', '#f0883e', '#d29922', '#bc8cff', '#79c0ff'
        ];

        function createCharts(currentPhase) {{
            const phaseObj = chartData[currentPhase];
            if (!phaseObj) return;

            let allWorkers = new Set();
            Object.values(phaseObj).forEach(s => s.workers.forEach(w => allWorkers.add(w)));
            const sortedWorkers = Array.from(allWorkers).sort((a,b) => a - b);

            // 1. Scalabilità (T1 / Tp)
            const scalabilityDatasets = [];
            let colorIdx = 0;
            for (const [label, data] of Object.entries(phaseObj)) {{
                const color = colors[colorIdx % colors.length];
                colorIdx++;
                scalabilityDatasets.push({{
                    label: label,
                    data: data.workers.map((w, i) => ({{ x: w, y: data.speedup_t1[i] }})),
                    borderColor: color,
                    backgroundColor: color,
                    tension: 0.1,
                    pointRadius: 4
                }});
            }}

            scalabilityDatasets.push({{
                label: 'Ideal (Lineare: S=p)',
                data: sortedWorkers.map(w => ({{ x: w, y: w }})),
                borderColor: '#8b949e',
                borderDash: [5, 5],
                pointRadius: 0,
                fill: false
            }});

            if (amdahlData && amdahlData[currentPhase]) {{
                const fVal = amdahlData[currentPhase];
                const pctStr = (fVal * 100).toFixed(2);
                scalabilityDatasets.push({{
                    label: `Amdahl (f=${{pctStr}}%)`,
                    data: sortedWorkers.map(w => ({{
                        x: w,
                        y: +(1.0 / (fVal + (1.0 - fVal) / w)).toFixed(2)
                    }})),
                    borderColor: '#f0883e',
                    borderDash: [3, 3],
                    pointRadius: 0,
                    fill: false
                }});
            }}

            // 2. Speedup Assoluto (Tseq / Tp)
            const speedupDatasets = [];
            colorIdx = 0;
            for (const [label, data] of Object.entries(phaseObj)) {{
                const color = colors[colorIdx % colors.length];
                colorIdx++;
                speedupDatasets.push({{
                    label: label,
                    data: data.workers.map((w, i) => ({{ x: w, y: data.speedup_seq[i] }})),
                    borderColor: color,
                    backgroundColor: color,
                    tension: 0.1,
                    pointRadius: 4
                }});
            }}

            speedupDatasets.push({{
                label: 'Ideal (Lineare: S=p)',
                data: sortedWorkers.map(w => ({{ x: w, y: w }})),
                borderColor: '#8b949e',
                borderDash: [5, 5],
                pointRadius: 0,
                fill: false
            }});

            if (amdahlData && amdahlData[currentPhase]) {{
                const fVal = amdahlData[currentPhase];
                const pctStr = (fVal * 100).toFixed(2);
                speedupDatasets.push({{
                    label: `Amdahl (f=${{pctStr}}%)`,
                    data: sortedWorkers.map(w => ({{
                        x: w,
                        y: +(1.0 / (fVal + (1.0 - fVal) / w)).toFixed(2)
                    }})),
                    borderColor: '#f0883e',
                    borderDash: [3, 3],
                    pointRadius: 0,
                    fill: false
                }});
            }}

            // 3. Efficienza Parallela
            const effDatasets = [];
            colorIdx = 0;
            for (const [label, data] of Object.entries(phaseObj)) {{
                const color = colors[colorIdx % colors.length];
                colorIdx++;
                effDatasets.push({{
                    label: label,
                    data: data.workers.map((w, i) => ({{ x: w, y: data.efficiency_t1[i] }})),
                    borderColor: color,
                    backgroundColor: color,
                    tension: 0.1,
                    pointRadius: 4
                }});
            }}
            effDatasets.push({{
                label: 'Ideal (100%)',
                data: sortedWorkers.map(w => ({{ x: w, y: 100 }})),
                borderColor: '#8b949e',
                borderDash: [5, 5],
                pointRadius: 0,
                fill: false
            }});

            // 4. Completion Time
            const timeDatasets = [];
            colorIdx = 0;
            for (const [label, data] of Object.entries(phaseObj)) {{
                const color = colors[colorIdx % colors.length];
                colorIdx++;
                timeDatasets.push({{
                    label: label,
                    data: data.workers.map((w, i) => ({{ x: w, y: data.times[i] }})),
                    borderColor: color,
                    backgroundColor: color,
                    tension: 0.1,
                    pointRadius: 4
                }});
            }}

            if (scalabilityChart) scalabilityChart.destroy();
            scalabilityChart = new Chart(document.getElementById('scalabilityChart'), {{
                type: 'line',
                data: {{ datasets: scalabilityDatasets }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        x: {{ type: 'linear', title: {{ display: true, text: 'Cores / Workers', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
                        y: {{ type: 'linear', title: {{ display: true, text: 'Scalabilità (T1 / Tp)', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }}
                    }},
                    plugins: {{ legend: {{ labels: {{ color: '#c9d1d9' }} }} }}
                }}
            }});

            if (speedupChart) speedupChart.destroy();
            speedupChart = new Chart(document.getElementById('speedupChart'), {{
                type: 'line',
                data: {{ datasets: speedupDatasets }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        x: {{ type: 'linear', title: {{ display: true, text: 'Cores / Workers', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
                        y: {{ type: 'linear', title: {{ display: true, text: 'Speedup Assoluto (Tseq / Tp)', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }}
                    }},
                    plugins: {{ legend: {{ labels: {{ color: '#c9d1d9' }} }} }}
                }}
            }});

            if (efficiencyChart) efficiencyChart.destroy();
            efficiencyChart = new Chart(document.getElementById('efficiencyChart'), {{
                type: 'line',
                data: {{ datasets: effDatasets }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        x: {{ type: 'linear', title: {{ display: true, text: 'Cores / Workers', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
                        y: {{ type: 'linear', min: 0, max: 110, title: {{ display: true, text: 'Efficienza (%)', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }}
                    }},
                    plugins: {{ legend: {{ labels: {{ color: '#c9d1d9' }} }} }}
                }}
            }});

            if (timeChart) timeChart.destroy();
            timeChart = new Chart(document.getElementById('timeChart'), {{
                type: 'line',
                data: {{ datasets: timeDatasets }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        x: {{ type: 'linear', title: {{ display: true, text: 'Cores / Workers', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
                        y: {{ type: 'linear', title: {{ display: true, text: 'Secondi (s)', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }}
                    }},
                    plugins: {{ legend: {{ labels: {{ color: '#c9d1d9' }} }} }}
                }}
            }});

            const tbody = document.getElementById('tableBody');
            tbody.innerHTML = '';
            for (const [label, data] of Object.entries(phaseObj)) {{
                data.workers.forEach((w, i) => {{
                    const tr = document.createElement('tr');
                    const parts = label.replace(')', '').split(' (');
                    const vName = parts[0];
                    const vMode = parts[1] || '-';
                    const tVal = data.times[i] !== null ? data.times[i].toFixed(3) : '-';
                    const sSeqVal = data.speedup_seq[i] !== null ? data.speedup_seq[i].toFixed(2) + 'x' : '-';
                    const sT1Val = data.speedup_t1[i] !== null ? data.speedup_t1[i].toFixed(2) + 'x' : '-';
                    const eVal = data.efficiency_t1[i] !== null ? data.efficiency_t1[i].toFixed(1) + '%' : '-';

                    tr.innerHTML = `
                        <td><strong>${{vName}}</strong></td>
                        <td>${{vMode}}</td>
                        <td><span class="badge badge-primary">${{w}}</span></td>
                        <td>${{tVal}} s</td>
                        <td>${{sSeqVal}}</td>
                        <td>${{sT1Val}}</td>
                        <td>${{eVal}}</td>
                    `;
                    tbody.appendChild(tr);
                }});
            }}
        }}

        function updateCharts() {{
            createCharts(selectElem.value);
        }}

        createCharts(selectElem.value);
    </script>
</body>
</html>
"""
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"Report interattivo HTML salvato in: {html_path}")


def generate_matplotlib_figures(summary_rows, results_dir, config=None):
    if not HAS_MATPLOTLIB:
        print(
            "Note: matplotlib non installato. Figure PNG saltate (usati i file vettoriali SVG e report.html)."
        )
        return

    if config is None:
        config = PlotConfig(output_dir=results_dir)

    target_dir = config.output_dir
    os.makedirs(target_dir, exist_ok=True)

    phases = sorted({r["phase"] for r in summary_rows})
    target_phase = "nsga2_total" if "nsga2_total" in phases else phases[0]
    rows = [r for r in summary_rows if r["phase"] == target_phase]
    if not rows:
        return

    series = defaultdict(list)
    for r in rows:
        key = f"{r['variant']} ({r['mode']})"
        series[key].append(
            {
                "w": int(r["workers"]),
                "time_s": float(r["avg_mean_ms"]) / 1000.0,
                "s_t1": float(r["speedup_t1"]) if r.get("speedup_t1") else None,
                "s_seq": float(r["speedup_seq"]) if r.get("speedup_seq") else None,
                "eff_t1": float(r["efficiency_t1"]) * 100.0
                if r.get("efficiency_t1")
                else None,
            }
        )

    for k in series:
        series[k].sort(key=lambda x: x["w"])

    all_workers = sorted({x["w"] for s in series.values() for x in s})

    amdahl_f = None
    if getattr(config, "show_amdahl", True):
        if getattr(config, "amdahl_f", None) is not None:
            amdahl_f = config.amdahl_f
        else:
            amdahl_f = compute_amdahl_fraction(summary_rows, target_phase)

    # 1. Speedup Assoluto PNG (T_seq / T_p)
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for label, items in series.items():
        ws = [x["w"] for x in items if x["s_seq"] is not None]
        ss = [x["s_seq"] for x in items if x["s_seq"] is not None]
        if ws:
            ax.plot(ws, ss, marker="o", label=label, linewidth=2)

    if all_workers:
        ax.plot(
            all_workers,
            all_workers,
            "--",
            color="gray",
            label="Ideale (S = p)",
            linewidth=1.5,
        )
        if amdahl_f is not None and 0.0 < amdahl_f < 1.0:
            step_count = 100
            min_w, max_w = min(all_workers), max(all_workers)
            dense_w = [min_w + i * (max_w - min_w) / step_count for i in range(step_count + 1)]
            amdahl_vals = [1.0 / (amdahl_f + (1.0 - amdahl_f) / w) for w in dense_w]
            ax.plot(
                dense_w,
                amdahl_vals,
                ":",
                color="#e36209",
                label=f"Amdahl (f={amdahl_f*100:.2f}%)",
                linewidth=1.8,
            )

    ax.set_xlabel("Cores / Workers", fontsize=12)
    ax.set_ylabel("Speedup Assoluto (Tseq / Tp)", fontsize=12)
    speedup_title = config.get_title(
        f"Speedup Assoluto vs Cores ({target_phase})", config.speedup_title
    )
    ax.set_title(speedup_title, fontsize=14, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(config.path(config.speedup_name, "png"))
    plt.close(fig)

    # 1b. Scalabilità PNG (T_1 / T_p)
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for label, items in series.items():
        ws = [x["w"] for x in items if x["s_t1"] is not None]
        ss = [x["s_t1"] for x in items if x["s_t1"] is not None]
        if ws:
            ax.plot(ws, ss, marker="o", label=label, linewidth=2)

    if all_workers:
        ax.plot(
            all_workers,
            all_workers,
            "--",
            color="gray",
            label="Ideale (Scalab = p)",
            linewidth=1.5,
        )
        if amdahl_f is not None and 0.0 < amdahl_f < 1.0:
            step_count = 100
            min_w, max_w = min(all_workers), max(all_workers)
            dense_w = [min_w + i * (max_w - min_w) / step_count for i in range(step_count + 1)]
            amdahl_vals = [1.0 / (amdahl_f + (1.0 - amdahl_f) / w) for w in dense_w]
            ax.plot(
                dense_w,
                amdahl_vals,
                ":",
                color="#e36209",
                label=f"Amdahl (f={amdahl_f*100:.2f}%)",
                linewidth=1.8,
            )

    ax.set_xlabel("Cores / Workers", fontsize=12)
    ax.set_ylabel("Scalabilità (T1 / Tp)", fontsize=12)
    scalability_title = config.get_title(
        f"Scalabilità vs Cores ({target_phase})", config.scalability_title
    )
    ax.set_title(scalability_title, fontsize=14, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(config.path(config.scalability_name, "png"))
    plt.close(fig)

    # 2. Efficienza PNG
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for label, items in series.items():
        ws = [x["w"] for x in items if x["eff_t1"] is not None]
        es = [x["eff_t1"] for x in items if x["eff_t1"] is not None]
        if ws:
            ax.plot(ws, es, marker="s", label=label, linewidth=2)

    if all_workers:
        ax.axhline(
            100.0, linestyle="--", color="gray", label="Ideal (100%)", linewidth=1.5
        )

    ax.set_xlabel("Cores / Workers", fontsize=12)
    ax.set_ylabel("Efficienza Parallela (%)", fontsize=12)
    ax.set_ylim(0, 110)
    efficiency_title = config.get_title(
        f"Efficienza Parallela vs Cores ({target_phase})",
        config.efficiency_title,
    )
    ax.set_title(efficiency_title, fontsize=14, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(config.path(config.efficiency_name, "png"))
    plt.close(fig)

    # 3. Tempi PNG
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for label, items in series.items():
        ws = [x["w"] for x in items]
        ts = [x["time_s"] for x in items]
        ax.plot(ws, ts, marker="^", label=label, linewidth=2)

    ax.set_xlabel("Cores / Workers", fontsize=12)
    ax.set_ylabel("Completion time (secondi)", fontsize=12)
    time_title = config.get_title(
        f"Completion time vs Cores ({target_phase})", config.time_title
    )
    ax.set_title(time_title, fontsize=14, fontweight="bold")
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.grid(True, which="both", ls="--", alpha=0.6)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(config.path(config.time_name, "png"))
    plt.close(fig)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fmo_dir = os.path.abspath(os.path.join(script_dir, ".."))

    parser = argparse.ArgumentParser(
        description="Analisi e visualizzazione delle prestazioni di FMO (NSGA-II)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Esempi di utilizzo:
  # Esecuzione standard sulla cartella results/
  python3 plot_benchmark.py results/

  # Aggiungi un prefisso a tutti i nomi dei file generati (es. exp1_speedup.svg)
  python3 plot_benchmark.py results/ --prefix exp1_

  # Aggiungi un suffisso ai file generati (es. speedup_prostate.svg)
  python3 plot_benchmark.py results/ --suffix _prostate

  # Personalizza prefisso file e titolo dei grafici
  python3 plot_benchmark.py results/ --prefix run1_ --title-prefix "[Prostate N=256] "

  # Personalizza singoli nomi di file
  python3 plot_benchmark.py results/ --speedup-name speedup_cores --efficiency-name eff_cores

  # Personalizza i titoli dei grafici
  python3 plot_benchmark.py results/ --speedup-title "Speedup NSGA-II su Prostate"

  # Salva i grafici con nome personalizzato in un'altra cartella
  python3 plot_benchmark.py results/ -o ./my_plots/ --prefix fig_ --theme light
""",
    )
    parser.add_argument(
        "target",
        nargs="?",
        default=None,
        help="Cartella dei risultati o file CSV (default: cartella più recente in results/)",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default=None,
        help="Cartella dove salvare i grafici e i report (default: stessa cartella dei dati)",
    )
    parser.add_argument(
        "-t",
        "--theme",
        choices=["light", "dark"],
        default="light",
        help="Tema grafico SVG: 'light' (default) o 'dark'",
    )

    # Naming group
    naming = parser.add_argument_group("Personalizzazione Nomi File")
    naming.add_argument(
        "-p",
        "--prefix",
        default="",
        help="Prefisso per tutti i file generati (es. 'run1_', 'prostate_')",
    )
    naming.add_argument(
        "-s",
        "--suffix",
        default="",
        help="Suffisso per tutti i file generati (es. '_pop256', '_v1')",
    )
    naming.add_argument(
        "--speedup-name",
        default="speedup",
        help="Nome base file per speedup (default: 'speedup')",
    )
    naming.add_argument(
        "--scalability-name",
        default="scalability",
        help="Nome base file per scalabilità (default: 'scalability')",
    )
    naming.add_argument(
        "--efficiency-name",
        default="efficiency",
        help="Nome base file per efficienza (default: 'efficiency')",
    )
    naming.add_argument(
        "--time-name",
        default="execution_time",
        help="Nome base file per tempi di completamento (default: 'execution_time')",
    )
    naming.add_argument(
        "--stability-name",
        default="solution_stability",
        help="Nome base file per stabilità/RMSE (default: 'solution_stability')",
    )
    naming.add_argument(
        "--report-name",
        default="report",
        help="Nome base file per il report HTML (default: 'report')",
    )

    # Title group
    titles = parser.add_argument_group("Personalizzazione Titoli Grafici")
    titles.add_argument(
        "--title-prefix",
        default="",
        help="Prefisso da anteporre ai titoli di tutti i grafici (es. '[Prostate 256] ')",
    )
    titles.add_argument(
        "--title-suffix",
        default="",
        help="Suffisso da posporre ai titoli di tutti i grafici (es. ' - NSGA-II')",
    )
    titles.add_argument(
        "--speedup-title",
        default=None,
        help="Titolo personalizzato per il grafico dello speedup",
    )
    titles.add_argument(
        "--scalability-title",
        default=None,
        help="Titolo personalizzato per il grafico della scalabilità",
    )
    titles.add_argument(
        "--efficiency-title",
        default=None,
        help="Titolo personalizzato per il grafico dell'efficienza",
    )
    titles.add_argument(
        "--time-title",
        default=None,
        help="Titolo personalizzato per il grafico del tempo di completamento",
    )
    titles.add_argument(
        "--stability-title",
        default=None,
        help="Titolo personalizzato per il grafico di stabilità",
    )

    # Amdahl group
    amdahl_grp = parser.add_argument_group("Analisi Teorica di Amdahl")
    amdahl_grp.add_argument(
        "--amdahl-f",
        "--amhdahl-f",
        dest="amdahl_f",
        type=float,
        default=None,
        help="Frazione seriale non parallelizzabile (f) per il limite di Amdahl (es. 0.0106 oppure 1.06%%). Se omessa, viene calcolata automaticamente dal baseline sequenziale.",
    )
    amdahl_grp.add_argument(
        "--no-amdahl",
        action="store_true",
        help="Disabilita il tracciamento della curva teorica di Amdahl nei grafici di speedup.",
    )

    args = parser.parse_args()

    if args.amdahl_f is not None and args.amdahl_f > 1.0:
        # Se l'utente ha inserito una percentuale (es. 1.04 o 1.06 anziché 0.0104 o 0.0106), normalizza in frazione [0, 1]
        args.amdahl_f = args.amdahl_f / 100.0

    target = args.target
    if not target:
        target = find_latest_results_dir(fmo_dir)

    if not target or not os.path.exists(target):
        parser.print_help()
        print("\nErrore: Nessuna directory dei risultati trovata.")
        sys.exit(1)

    summary_rows, raw_rows, stab_rows, results_dir = load_data(target)
    output_dir = args.output_dir if args.output_dir else results_dir

    config = PlotConfig(
        output_dir=output_dir,
        prefix=args.prefix,
        suffix=args.suffix,
        speedup_name=args.speedup_name,
        scalability_name=args.scalability_name,
        efficiency_name=args.efficiency_name,
        time_name=args.time_name,
        stability_name=args.stability_name,
        report_name=args.report_name,
        title_prefix=args.title_prefix,
        title_suffix=args.title_suffix,
        speedup_title=args.speedup_title,
        scalability_title=args.scalability_title,
        efficiency_title=args.efficiency_title,
        time_title=args.time_title,
        stability_title=args.stability_title,
        theme=args.theme,
        amdahl_f=args.amdahl_f,
        show_amdahl=(not args.no_amdahl),
    )

    print(f"Caricamento dati da: {target}")
    print(f"Salvataggio grafici in: {config.output_dir} (Tema: {config.theme})")
    if config.prefix or config.suffix:
        print(f"Pattern nomi file: {config.prefix}<nome>{config.suffix}.svg")

    print_terminal_summary(summary_rows, stab_rows, config=config)
    generate_all_svg_charts(
        summary_rows, results_dir, stab_rows, theme=config.theme, config=config
    )
    generate_html_report(summary_rows, raw_rows, results_dir, stab_rows, config=config)
    generate_matplotlib_figures(summary_rows, results_dir, config=config)


if __name__ == "__main__":
    main()
