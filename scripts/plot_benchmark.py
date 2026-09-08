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

import sys
import os
import csv
import json
import math
from collections import defaultdict

# Controllo presenza matplotlib (opzionale)
HAS_MATPLOTLIB = False
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    pass


def find_latest_results_dir(base_dir):
    results_dir = os.path.join(base_dir, "results")
    if not os.path.exists(results_dir):
        return None
    runs = [os.path.join(results_dir, d) for d in os.listdir(results_dir) if d.startswith("run_")]
    if not runs:
        return None
    runs.sort(key=os.path.getmtime, reverse=True)
    return runs[0]


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
        with open(summary_csv, 'r', newline='') as f:
            reader = csv.DictReader(f)
            summary_rows = list(reader)

    raw_rows = []
    if os.path.exists(raw_csv):
        with open(raw_csv, 'r', newline='') as f:
            reader = csv.DictReader(f)
            raw_rows = list(reader)

    stab_summary_csv = os.path.join(results_dir, "stability_summary.csv")
    stab_raw_csv = os.path.join(results_dir, "stability_rmse.csv")

    stab_rows = []
    if os.path.exists(stab_summary_csv):
        try:
            with open(stab_summary_csv, 'r', newline='') as f:
                stab_rows = list(csv.DictReader(f))
        except Exception:
            pass
    elif os.path.exists(stab_raw_csv):
        try:
            with open(stab_raw_csv, 'r', newline='') as f:
                reader = csv.DictReader(f)
                groups = defaultdict(list)
                for r in reader:
                    w_val = int(r['workers'])
                    groups[(r['variant'], r['mode'], w_val)].append({
                        'rmse_ptv': float(r['rmse_ptv']),
                        'rmse_rectum': float(r['rmse_rectum']),
                        'rmse_bladder': float(r['rmse_bladder']),
                        'rmse_total_fitness': float(r['rmse_total_fitness']),
                        'rmse_genes': float(r['rmse_genes']),
                        'spread_par': float(r['spread_par']),
                    })
                for (v, m, w_val), items in sorted(groups.items()):
                    c = len(items)
                    avg_fit = sum(x['rmse_total_fitness'] for x in items) / c
                    avg_gen = sum(x['rmse_genes'] for x in items) / c
                    avg_ptv = sum(x['rmse_ptv'] for x in items) / c
                    avg_rec = sum(x['rmse_rectum'] for x in items) / c
                    avg_bla = sum(x['rmse_bladder'] for x in items) / c
                    avg_spread = sum(x['spread_par'] for x in items) / c
                    status = "BASELINE" if v == 'seq' else ("DETERMINISTICO (0 err)" if avg_fit == 0.0 and avg_gen == 0.0 else "STABILE (Dev. GA)")
                    stab_rows.append({
                        'variant': v,
                        'mode': m,
                        'workers': str(w_val),
                        'repetitions': str(c),
                        'avg_rmse_ptv': f"{avg_ptv:.6f}",
                        'avg_rmse_rectum': f"{avg_rec:.6f}",
                        'avg_rmse_bladder': f"{avg_bla:.6f}",
                        'avg_rmse_fitness': f"{avg_fit:.6f}",
                        'avg_rmse_genes': f"{avg_gen:.6f}",
                        'avg_spread': f"{avg_spread:.6f}",
                        'status': status
                    })
        except Exception:
            pass

    return summary_rows, raw_rows, stab_rows, results_dir


def print_terminal_summary(summary_rows, stab_rows=None):
    if not summary_rows and not stab_rows:
        print("Nessun dato di riepilogo disponibile.")
        return

    if summary_rows:
        phases = sorted(list({r['phase'] for r in summary_rows}))
        target_phase = 'nsga2_total' if 'nsga2_total' in phases else phases[0]

        print("\n" + "=" * 105)
        print(f"  RIEPILOGO PRESTAZIONI (Fase: {target_phase})")
        print("=" * 105)
        print(f"{'VARIANTE':<14} {'MODALITÀ':<14} {'WORKERS':<9} {'TEMPO MEDIO (s)':<18} {'SPEEDUP (seq)':<16} {'SPEEDUP (T1)':<16} {'EFFICIENZA':<14}")
        print("-" * 105)

        def sort_key(r):
            return (r['variant'], r['mode'], int(r['workers']))

        phase_rows = [r for r in summary_rows if r['phase'] == target_phase]
        phase_rows.sort(key=sort_key)

        for r in phase_rows:
            variant = r['variant']
            mode = r['mode']
            workers = r['workers']
            mean_s = float(r['avg_mean_ms']) / 1000.0

            s_seq = f"{float(r['speedup_seq']):.2f}x" if r.get('speedup_seq') else "-"
            s_t1 = f"{float(r['speedup_t1']):.2f}x" if r.get('speedup_t1') else "-"
            eff = f"{float(r['efficiency_t1']) * 100:.1f}%" if r.get('efficiency_t1') else "-"

            print(f"{variant:<14} {mode:<14} {workers:<9} {mean_s:<18.3f} {s_seq:<16} {s_t1:<16} {eff:<14}")

        print("=" * 105 + "\n")

    if stab_rows:
        print("=" * 110)
        print("  STABILITÀ DELLE SOLUZIONI (RMSE vs Baseline Sequenziale) & DETERMINISMO")
        print("=" * 110)
        print(f"{'VARIANTE':<10} {'MODE':<10} {'WORKERS':<9} {'RMSE_FITNESS':<14} {'RMSE_GENI':<14} {'SPREAD_PAR':<12} {'STATO DETERMINISMO':<30}")
        print("-" * 110)
        for r in stab_rows:
            fit_val = float(r.get('avg_rmse_fitness', 0.0))
            gen_val = float(r.get('avg_rmse_genes', 0.0))
            spread_val = float(r.get('avg_spread', 0.0))
            print(f"{r['variant']:<10} {r['mode']:<10} {r['workers']:<9} {fit_val:<14.6f} {gen_val:<14.6f} {spread_val:<12.4f} {r.get('status', '-'):<30}")
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
    theme="light"
):
    """
    Genera un grafico vettoriale SVG professionale e pulito, pronto per report e tesi.
    Supporta tema chiaro (sfondo bianco puro, ideale per report/stampe) e scuro,
    spaziatura uniforme log2 per le potenze di 2 e marcatori distinti per ogni serie.
    """
    w, h = 800, 480
    pad_left, pad_right, pad_top, pad_bottom = 75, 185, 55, 60
    plot_w = w - pad_left - pad_right
    plot_h = h - pad_top - pad_bottom

    # Palette ad alto contrasto adatta sia a schermo che a stampa
    palette = [
        {"color": "#0969da", "shape": "circle"},     # Blu
        {"color": "#1a7f37", "shape": "rect"},       # Verde
        {"color": "#d1242f", "shape": "triangle"},   # Rosso
        {"color": "#8250df", "shape": "diamond"},    # Viola
        {"color": "#e36209", "shape": "circle"},     # Arancione
        {"color": "#0891b2", "shape": "rect"},       # Ciano
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

    unique_x = sorted(list(set(all_x)))
    min_x, max_x = min(unique_x), max(unique_x)

    # Rilevamento scala log2 su X (evita sovrapposizioni su 1, 2, 4...)
    use_log_x = (min_x > 0 and (max_x / min_x >= 4) and all(x > 0 for x in unique_x))
    log2_min_x = math.log2(min_x) if use_log_x else 0
    log2_max_x = math.log2(max_x) if use_log_x else 0

    def to_px(val_x):
        if use_log_x:
            if log2_max_x == log2_min_x:
                return pad_left + plot_w / 2
            return pad_left + (math.log2(val_x) - log2_min_x) / (log2_max_x - log2_min_x) * plot_w
        else:
            if max_x == min_x:
                return pad_left + plot_w / 2
            return pad_left + (val_x - min_x) / (max_x - min_x) * plot_w

    use_log_y = False
    if is_speedup and speedup_mode == "log":
        use_log_y = True
        min_y = 1.0
        max_y = float(max_x)
        log2_min_y = 0.0
        log2_max_y = math.log2(max_y)
    elif is_speedup and speedup_mode == "linear":
        min_y = 0.0
        max_actual_y = max(all_y)
        max_y = math.ceil(max_actual_y * 1.15 / 5.0) * 5.0
        if max_y < 2.0: max_y = 2.0
    elif is_efficiency:
        min_y = 0.0
        max_y = 105.0
    elif is_time and time_mode == "log":
        use_log_y = True
        min_y = max(1.0, 10 ** math.floor(math.log10(min(all_y))))
        max_y = 10 ** math.ceil(math.log10(max(all_y)))
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
            return pad_top + plot_h - (math.log2(safe_y) - log2_min_y) / (log2_max_y - log2_min_y) * plot_h
        elif use_log_y and is_time:
            safe_y = max(val_y, min_y)
            return pad_top + plot_h - (math.log10(safe_y) - log10_min_y) / (log10_max_y - log10_min_y) * plot_h
        else:
            if max_y == min_y:
                return pad_top + plot_h / 2
            return pad_top + plot_h - (val_y - min_y) / (max_y - min_y) * plot_h

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">',
        f'<rect width="{w}" height="{h}" fill="{bg_canvas}" />',
        f'<text x="{pad_left}" y="36" fill="{text_title}" font-family="system-ui, -apple-system, sans-serif" font-size="16" font-weight="bold">{title}</text>',
        f'<rect x="{pad_left}" y="{pad_top}" width="{plot_w}" height="{plot_h}" fill="{bg_plot}" stroke="{border_plot}" stroke-width="1.2" />'
    ]

    # Asse Y
    if use_log_y and is_speedup:
        powers = range(int(log2_min_y), int(log2_max_y) + 1)
        for p in powers:
            val_y = 2 ** p
            py = to_py(val_y)
            svg.append(f'<line x1="{pad_left}" y1="{py:.1f}" x2="{pad_left + plot_w}" y2="{py:.1f}" stroke="{grid_major}" stroke-width="1" stroke-dasharray="3,3" />')
            svg.append(f'<text x="{pad_left - 10}" y="{py + 4:.1f}" fill="{text_tick}" font-family="system-ui, -apple-system, sans-serif" font-size="11" text-anchor="end">{int(val_y)}x</text>')
    elif use_log_y and is_time:
        for p in range(int(log10_min_y), int(log10_max_y) + 1):
            val_y = 10 ** p
            py = to_py(val_y)
            svg.append(f'<line x1="{pad_left}" y1="{py:.1f}" x2="{pad_left + plot_w}" y2="{py:.1f}" stroke="{grid_major}" stroke-width="1" stroke-dasharray="3,3" />')
            svg.append(f'<text x="{pad_left - 10}" y="{py + 4:.1f}" fill="{text_tick}" font-family="system-ui, -apple-system, sans-serif" font-size="11" text-anchor="end">{int(val_y)}s</text>')
    else:
        y_ticks = 5
        for i in range(y_ticks + 1):
            y_val = min_y + (max_y - min_y) * (i / y_ticks)
            py = to_py(y_val)
            svg.append(f'<line x1="{pad_left}" y1="{py:.1f}" x2="{pad_left + plot_w}" y2="{py:.1f}" stroke="{grid_major}" stroke-width="1" stroke-dasharray="3,3" />')
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
            svg.append(f'<text x="{pad_left - 10}" y="{py + 4:.1f}" fill="{text_tick}" font-family="system-ui, -apple-system, sans-serif" font-size="11" text-anchor="end">{y_text}</text>')

    # Asse X
    for ux in unique_x:
        px = to_px(ux)
        svg.append(f'<line x1="{px:.1f}" y1="{pad_top}" x2="{px:.1f}" y2="{pad_top + plot_h}" stroke="{grid_minor}" stroke-width="1" stroke-dasharray="2,2" />')
        svg.append(f'<line x1="{px:.1f}" y1="{pad_top + plot_h}" x2="{px:.1f}" y2="{pad_top + plot_h + 5}" stroke="{border_plot}" stroke-width="1.2" />')
        svg.append(f'<text x="{px:.1f}" y="{pad_top + plot_h + 18}" fill="{text_tick}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500" text-anchor="middle">{ux}</text>')

    # Linea ideale (Speedup o Efficienza)
    if is_speedup:
        if use_log_y:
            p1_x, p1_y = to_px(min_x), to_py(min_x)
            p2_x, p2_y = to_px(max_x), to_py(max_x)
            svg.append(f'<line x1="{p1_x:.1f}" y1="{p1_y:.1f}" x2="{p2_x:.1f}" y2="{p2_y:.1f}" stroke="{ideal_color}" stroke-width="1.6" stroke-dasharray="5,4" />')
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
                svg.append(f'<polyline fill="none" stroke="{ideal_color}" stroke-width="1.6" stroke-dasharray="5,4" points="{poly}" />')
    elif is_efficiency:
        py_100 = to_py(100.0)
        svg.append(f'<line x1="{pad_left}" y1="{py_100:.1f}" x2="{pad_left + plot_w}" y2="{py_100:.1f}" stroke="{ideal_color}" stroke-width="1.6" stroke-dasharray="5,4" />')

    # Tracciamento serie dati
    legend_entries = []
    for idx, (label, points) in enumerate(series_data.items()):
        style = palette[idx % len(palette)]
        color = style["color"]
        shape = style["shape"]
        valid_pts = [(px, py) for px, py in points if px is not None and py is not None]
        valid_pts.sort(key=lambda p: p[0])

        if len(valid_pts) > 1:
            poly_points = " ".join([f"{to_px(px):.1f},{to_py(py):.1f}" for px, py in valid_pts])
            svg.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.5" stroke-linecap="round" stroke-linejoin="round" points="{poly_points}" />')

        for px, py in valid_pts:
            cx, cy = to_px(px), to_py(py)
            if shape == "rect":
                svg.append(f'<rect x="{cx - 4:.1f}" y="{cy - 4:.1f}" width="8" height="8" fill="{color}" stroke="#ffffff" stroke-width="1.5" />')
            elif shape == "triangle":
                svg.append(f'<polygon points="{cx:.1f},{cy - 5.5:.1f} {cx + 5:.1f},{cy + 4.5:.1f} {cx - 5:.1f},{cy + 4.5:.1f}" fill="{color}" stroke="#ffffff" stroke-width="1.5" />')
            elif shape == "diamond":
                svg.append(f'<polygon points="{cx:.1f},{cy - 5.5:.1f} {cx + 5:.1f},{cy:.1f} {cx:.1f},{cy + 5.5:.1f} {cx - 5:.1f},{cy:.1f}" fill="{color}" stroke="#ffffff" stroke-width="1.5" />')
            else:
                svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4.5" fill="{color}" stroke="#ffffff" stroke-width="1.5" />')

        legend_entries.append((label, color, shape))

    # Box Legenda
    legend_h = len(legend_entries) * 24 + (24 if (is_speedup or is_efficiency) else 10) + 12
    legend_x = pad_left + plot_w + 14
    legend_y = pad_top
    svg.append(f'<rect x="{legend_x}" y="{legend_y}" width="{pad_right - 24}" height="{legend_h}" fill="{card_bg}" stroke="{card_border}" stroke-width="1" rx="6" />')

    cur_y = legend_y + 18
    for label, color, shape in legend_entries:
        icon_cx = legend_x + 16
        icon_cy = cur_y - 4
        svg.append(f'<line x1="{icon_cx - 8}" y1="{icon_cy}" x2="{icon_cx + 8}" y2="{icon_cy}" stroke="{color}" stroke-width="2" />')
        if shape == "rect":
            svg.append(f'<rect x="{icon_cx - 3.5}" y="{icon_cy - 3.5}" width="7" height="7" fill="{color}" stroke="#ffffff" stroke-width="1" />')
        elif shape == "triangle":
            svg.append(f'<polygon points="{icon_cx},{icon_cy - 4} {icon_cx + 4},{icon_cy + 3.5} {icon_cx - 4},{icon_cy + 3.5}" fill="{color}" stroke="#ffffff" stroke-width="1" />')
        elif shape == "diamond":
            svg.append(f'<polygon points="{icon_cx},{icon_cy - 4.5} {icon_cx + 4},{icon_cy} {icon_cx},{icon_cy + 4.5} {icon_cx - 4},{icon_cy}" fill="{color}" stroke="#ffffff" stroke-width="1" />')
        else:
            svg.append(f'<circle cx="{icon_cx}" cy="{icon_cy}" r="3.5" fill="{color}" stroke="#ffffff" stroke-width="1" />')

        svg.append(f'<text x="{legend_x + 32}" y="{cur_y}" fill="{legend_text}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500">{label}</text>')
        cur_y += 22

    if is_speedup:
        icon_cx = legend_x + 16
        icon_cy = cur_y - 4
        svg.append(f'<line x1="{icon_cx - 8}" y1="{icon_cy}" x2="{icon_cx + 8}" y2="{icon_cy}" stroke="{ideal_color}" stroke-dasharray="3,2" stroke-width="1.8" />')
        svg.append(f'<text x="{legend_x + 32}" y="{cur_y}" fill="{ideal_color}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500">Ideale (S = p)</text>')
    elif is_efficiency:
        icon_cx = legend_x + 16
        icon_cy = cur_y - 4
        svg.append(f'<line x1="{icon_cx - 8}" y1="{icon_cy}" x2="{icon_cx + 8}" y2="{icon_cy}" stroke="{ideal_color}" stroke-dasharray="3,2" stroke-width="1.8" />')
        svg.append(f'<text x="{legend_x + 32}" y="{cur_y}" fill="{ideal_color}" font-family="system-ui, -apple-system, sans-serif" font-size="11" font-weight="500">Ideale (100%)</text>')

    # Etichette assi
    svg.append(f'<text x="{pad_left + plot_w / 2}" y="{h - 15}" fill="{text_axis}" font-family="system-ui, -apple-system, sans-serif" font-size="12" font-weight="600" text-anchor="middle">{x_label}</text>')
    svg.append(f'<text x="22" y="{pad_top + plot_h / 2}" fill="{text_axis}" font-family="system-ui, -apple-system, sans-serif" font-size="12" font-weight="600" text-anchor="middle" transform="rotate(-90 22 {pad_top + plot_h / 2})">{y_label}</text>')

    svg.append('</svg>')

    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(svg))

    print(f"Grafico SVG salvato: {out_path}")


def generate_all_svg_charts(summary_rows, results_dir, stab_rows=None, theme="light"):
    """
    Genera i grafici SVG standalone per la fase principale e per la stabilità.
    """
    phases = sorted(list({r['phase'] for r in summary_rows}))
    target_phase = 'nsga2_total' if 'nsga2_total' in phases else phases[0]
    rows = [r for r in summary_rows if r['phase'] == target_phase]

    speedup_series = defaultdict(list)
    eff_series = defaultdict(list)
    time_series = defaultdict(list)

    for r in rows:
        label = f"{r['variant']} ({r['mode']})"
        w = int(r['workers'])
        time_s = float(r['avg_mean_ms']) / 1000.0
        s_t1 = float(r['speedup_t1']) if r.get('speedup_t1') else None
        eff_t1 = float(r['efficiency_t1']) * 100.0 if r.get('efficiency_t1') else None

        time_series[label].append((w, time_s))
        if s_t1 is not None:
            speedup_series[label].append((w, s_t1))
        if eff_t1 is not None:
            eff_series[label].append((w, eff_t1))

    # 1. Speedup principale (Log-Log per rappresentazione standard HPC)
    generate_svg_chart(
        f"Speedup vs Cores ({target_phase})",
        "Numero di Cores / Workers",
        "Speedup (T1 / Tp)",
        speedup_series,
        os.path.join(results_dir, "speedup.svg"),
        is_speedup=True,
        speedup_mode="log",
        theme=theme
    )

    # 2. Speedup scala lineare (alternativa per report)
    generate_svg_chart(
        f"Speedup vs Cores - Scala Lineare ({target_phase})",
        "Numero di Cores / Workers",
        "Speedup (T1 / Tp)",
        speedup_series,
        os.path.join(results_dir, "speedup_linear.svg"),
        is_speedup=True,
        speedup_mode="linear",
        theme=theme
    )

    # 3. Efficienza Parallela
    generate_svg_chart(
        f"Efficienza Parallela vs Cores ({target_phase})",
        "Numero di Cores / Workers",
        "Efficienza Parallela (%)",
        eff_series,
        os.path.join(results_dir, "efficiency.svg"),
        is_efficiency=True,
        theme=theme
    )

    # 4. Tempo di Esecuzione lineare
    generate_svg_chart(
        f"Tempo di Esecuzione vs Cores ({target_phase})",
        "Numero di Cores / Workers",
        "Tempo (secondi)",
        time_series,
        os.path.join(results_dir, "execution_time.svg"),
        is_time=True,
        time_mode="linear",
        theme=theme
    )

    # 5. Tempo di Esecuzione logaritmico
    generate_svg_chart(
        f"Tempo di Esecuzione vs Cores - Scala Log ({target_phase})",
        "Numero di Cores / Workers",
        "Tempo (secondi, scala log10)",
        time_series,
        os.path.join(results_dir, "execution_time_log.svg"),
        is_time=True,
        time_mode="log",
        theme=theme
    )

    # 6. Stabilità delle soluzioni
    if stab_rows:
        stab_series = defaultdict(list)
        for r in stab_rows:
            if r['variant'] == 'seq':
                continue
            label = f"{r['variant']} ({r['mode']})"
            w = int(r['workers'])
            fit_rmse = float(r.get('avg_rmse_fitness', 0.0))
            stab_series[label].append((w, fit_rmse))

        if stab_series:
            generate_svg_chart(
                "Stabilità delle Soluzioni: RMSE Fitness vs Cores",
                "Numero di Cores / Workers",
                "RMSE Fitness (vs Baseline)",
                stab_series,
                os.path.join(results_dir, "solution_stability.svg"),
                is_stability=True,
                theme=theme
            )


def generate_html_report(summary_rows, raw_rows, results_dir, stab_rows=None):
    """
    Genera un report HTML completo e interattivo con grafici Chart.js (online)
    e link diretti agli SVG vettoriali (offline).
    """
    html_path = os.path.join(results_dir, "report.html")

    phase_data = defaultdict(lambda: defaultdict(dict))
    for r in summary_rows:
        phase = r['phase']
        label = f"{r['variant']} ({r['mode']})"
        w = int(r['workers'])
        time_s = float(r['avg_mean_ms']) / 1000.0
        s_seq = float(r['speedup_seq']) if r.get('speedup_seq') else None
        s_t1 = float(r['speedup_t1']) if r.get('speedup_t1') else None
        eff_t1 = float(r['efficiency_t1']) * 100.0 if r.get('efficiency_t1') else None

        phase_data[phase][label][w] = {
            'time_s': time_s,
            'speedup_seq': s_seq,
            'speedup_t1': s_t1,
            'efficiency_t1': eff_t1
        }

    sys_info = ""
    sys_info_path = os.path.join(results_dir, "system_info.txt")
    if os.path.exists(sys_info_path):
        with open(sys_info_path, 'r') as f:
            sys_info = f.read()

    chart_payload = {}
    for phase, series_dict in phase_data.items():
        chart_payload[phase] = {}
        for label, w_dict in series_dict.items():
            sorted_w = sorted(w_dict.keys())
            chart_payload[phase][label] = {
                'workers': sorted_w,
                'times': [w_dict[w]['time_s'] for w in sorted_w],
                'speedup_seq': [w_dict[w]['speedup_seq'] for w in sorted_w],
                'speedup_t1': [w_dict[w]['speedup_t1'] for w in sorted_w],
                'efficiency_t1': [w_dict[w]['efficiency_t1'] for w in sorted_w]
            }

    stab_card_html = ""
    stab_link_html = ""
    if stab_rows:
        stab_link_html = """
            <a href="solution_stability.svg" target="_blank">🎯 solution_stability.svg</a>
            <a href="stability_summary.csv" download>📄 stability_summary.csv</a>
            <a href="stability_rmse.csv" download>📄 stability_rmse.csv</a>
        """
        rows_html = []
        for r in stab_rows:
            st = r.get('status', '')
            if 'DETERMINISTICO' in st:
                badge = f'<span class="badge" style="background: rgba(63, 185, 80, 0.2); color: #3fb950;">{st}</span>'
            elif 'BASELINE' in st:
                badge = f'<span class="badge" style="background: rgba(88, 166, 255, 0.2); color: #58a6ff;">{st}</span>'
            else:
                badge = f'<span class="badge" style="background: rgba(210, 153, 34, 0.2); color: #d29922;">{st}</span>'

            rows_html.append(f"""
                <tr>
                    <td><strong>{r.get('variant', '')}</strong></td>
                    <td>{r.get('mode', '')}</td>
                    <td>{r.get('workers', '1')}</td>
                    <td>{r.get('repetitions', '1')}</td>
                    <td><code>{r.get('avg_rmse_fitness', '-')}</code></td>
                    <td><code>{r.get('avg_rmse_genes', '-')}</code></td>
                    <td><code>{r.get('avg_rmse_ptv', '-')}</code></td>
                    <td><code>{r.get('avg_rmse_rectum', '-')}</code></td>
                    <td><code>{r.get('avg_rmse_bladder', '-')}</code></td>
                    <td><code>{r.get('avg_spread', '-')}</code></td>
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
                        {''.join(rows_html)}
                    </tbody>
                </table>
            </div>
        </div>
        """

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
            <a href="speedup.svg" target="_blank">📈 speedup.svg (log)</a>
            <a href="speedup_linear.svg" target="_blank">📈 speedup_linear.svg</a>
            <a href="efficiency.svg" target="_blank">⚡ efficiency.svg</a>
            <a href="execution_time.svg" target="_blank">⏱️ execution_time.svg</a>
            <a href="execution_time_log.svg" target="_blank">⏱️ execution_time_log.svg</a>
            {stab_link_html}
            <a href="summary.csv" download>📄 summary.csv</a>
            <a href="timings_raw.csv" download>📄 timings_raw.csv</a>
        </div>

        <label for="phaseSelect" style="font-weight: 600;">Seleziona Fase da Visualizzare: </label>
        <select id="phaseSelect" onchange="updateCharts()"></select>

        <div class="grid">
            <div class="card">
                <h2>📈 Speedup vs Cores (S = T1 / Tp)</h2>
                <div class="chart-container">
                    <canvas id="speedupChart"></canvas>
                </div>
            </div>
            <div class="card">
                <h2>⚡ Efficienza Parallela vs Cores (E = S / p)</h2>
                <div class="chart-container">
                    <canvas id="efficiencyChart"></canvas>
                </div>
            </div>
            <div class="card">
                <h2>⏱️ Tempo di Esecuzione vs Cores</h2>
                <div class="chart-container">
                    <canvas id="timeChart"></canvas>
                </div>
            </div>
            <div class="card">
                <h2>📊 Tabella di Riepilogo Prestazioni</h2>
                <div style="overflow-x: auto; max-height: 360px;">
                    <table id="summaryTable">
                        <thead>
                            <tr>
                                <th>Variante</th>
                                <th>Modo</th>
                                <th>Cores</th>
                                <th>Tempo (s)</th>
                                <th>Speedup</th>
                                <th>Efficienza</th>
                            </tr>
                        </thead>
                        <tbody id="tableBody"></tbody>
                    </table>
                </div>
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
        const phases = Object.keys(chartData);

        const selectElem = document.getElementById('phaseSelect');
        phases.forEach(ph => {{
            const opt = document.createElement('option');
            opt.value = ph;
            opt.textContent = ph;
            if (ph === 'nsga2_total') opt.selected = true;
            selectElem.appendChild(opt);
        }});

        let speedupChart, efficiencyChart, timeChart;
        const colors = [
            '#58a6ff', '#3fb950', '#f0883e', '#d29922', '#bc8cff', '#79c0ff'
        ];

        function createCharts(currentPhase) {{
            const phaseObj = chartData[currentPhase];
            if (!phaseObj) return;

            let allWorkers = new Set();
            Object.values(phaseObj).forEach(s => s.workers.forEach(w => allWorkers.add(w)));
            const sortedWorkers = Array.from(allWorkers).sort((a,b) => a - b);

            const speedupDatasets = [];
            let colorIdx = 0;
            for (const [label, data] of Object.entries(phaseObj)) {{
                const color = colors[colorIdx % colors.length];
                colorIdx++;
                speedupDatasets.push({{
                    label: label,
                    data: data.workers.map((w, i) => ({{ x: w, y: data.speedup_t1[i] }})),
                    borderColor: color,
                    backgroundColor: color,
                    tension: 0.1,
                    pointRadius: 4
                }});
            }}

            speedupDatasets.push({{
                label: 'Ideale (Lineare: S=p)',
                data: sortedWorkers.map(w => ({{ x: w, y: w }})),
                borderColor: '#8b949e',
                borderDash: [5, 5],
                pointRadius: 0,
                fill: false
            }});

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
                label: 'Ideale (100%)',
                data: sortedWorkers.map(w => ({{ x: w, y: 100 }})),
                borderColor: '#8b949e',
                borderDash: [5, 5],
                pointRadius: 0,
                fill: false
            }});

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

            if (speedupChart) speedupChart.destroy();
            speedupChart = new Chart(document.getElementById('speedupChart'), {{
                type: 'line',
                data: {{ datasets: speedupDatasets }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        x: {{ type: 'linear', title: {{ display: true, text: 'Numero di Cores / Workers', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
                        y: {{ type: 'linear', title: {{ display: true, text: 'Speedup', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }}
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
                        x: {{ type: 'linear', title: {{ display: true, text: 'Numero di Cores / Workers', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
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
                        x: {{ type: 'linear', title: {{ display: true, text: 'Numero di Cores / Workers', color: '#8b949e' }}, grid: {{ color: '#30363d' }} }},
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
                    const sVal = data.speedup_t1[i] !== null ? data.speedup_t1[i].toFixed(2) + 'x' : '-';
                    const eVal = data.efficiency_t1[i] !== null ? data.efficiency_t1[i].toFixed(1) + '%' : '-';

                    tr.innerHTML = `
                        <td><strong>${{vName}}</strong></td>
                        <td>${{vMode}}</td>
                        <td><span class="badge badge-primary">${{w}}</span></td>
                        <td>${{tVal}} s</td>
                        <td>${{sVal}}</td>
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
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Report interattivo HTML salvato in: {html_path}")


def generate_matplotlib_figures(summary_rows, results_dir):
    if not HAS_MATPLOTLIB:
        print("Note: matplotlib non installato. Figure PNG saltate (usati i file vettoriali SVG e report.html).")
        return

    phases = sorted(list({r['phase'] for r in summary_rows}))
    target_phase = 'nsga2_total' if 'nsga2_total' in phases else phases[0]
    rows = [r for r in summary_rows if r['phase'] == target_phase]
    if not rows:
        return

    series = defaultdict(list)
    for r in rows:
        key = f"{r['variant']} ({r['mode']})"
        series[key].append({
            'w': int(r['workers']),
            'time_s': float(r['avg_mean_ms']) / 1000.0,
            's_t1': float(r['speedup_t1']) if r.get('speedup_t1') else None,
            'eff_t1': float(r['efficiency_t1']) * 100.0 if r.get('efficiency_t1') else None,
        })

    for k in series:
        series[k].sort(key=lambda x: x['w'])

    all_workers = sorted(list({x['w'] for s in series.values() for x in s}))

    # 1. Speedup PNG
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for label, items in series.items():
        ws = [x['w'] for x in items if x['s_t1'] is not None]
        ss = [x['s_t1'] for x in items if x['s_t1'] is not None]
        if ws:
            ax.plot(ws, ss, marker='o', label=label, linewidth=2)

    if all_workers:
        ax.plot(all_workers, all_workers, '--', color='gray', label='Ideale (S = p)', linewidth=1.5)

    ax.set_xlabel('Numero di Cores / Workers', fontsize=12)
    ax.set_ylabel('Speedup (T1 / Tp)', fontsize=12)
    ax.set_title(f'Speedup vs Cores ({target_phase})', fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "speedup.png"))
    plt.close(fig)

    # 2. Efficienza PNG
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for label, items in series.items():
        ws = [x['w'] for x in items if x['eff_t1'] is not None]
        es = [x['eff_t1'] for x in items if x['eff_t1'] is not None]
        if ws:
            ax.plot(ws, es, marker='s', label=label, linewidth=2)

    if all_workers:
        ax.axhline(100.0, linestyle='--', color='gray', label='Ideale (100%)', linewidth=1.5)

    ax.set_xlabel('Numero di Cores / Workers', fontsize=12)
    ax.set_ylabel('Efficienza Parallela (%)', fontsize=12)
    ax.set_ylim(0, 110)
    ax.set_title(f'Efficienza Parallela vs Cores ({target_phase})', fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.6)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "efficiency.png"))
    plt.close(fig)

    # 3. Tempi PNG
    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
    for label, items in series.items():
        ws = [x['w'] for x in items]
        ts = [x['time_s'] for x in items]
        ax.plot(ws, ts, marker='^', label=label, linewidth=2)

    ax.set_xlabel('Numero di Cores / Workers', fontsize=12)
    ax.set_ylabel('Tempo di Esecuzione (secondi)', fontsize=12)
    ax.set_title(f'Tempo di Esecuzione vs Cores ({target_phase})', fontsize=14, fontweight='bold')
    ax.set_xscale('log', base=2)
    ax.set_yscale('log')
    ax.grid(True, which="both", ls="--", alpha=0.6)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(os.path.join(results_dir, "execution_time.png"))
    plt.close(fig)


def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    fmo_dir = os.path.abspath(os.path.join(script_dir, ".."))

    target = None
    theme = "light"

    args = sys.argv[1:]
    i = 0
    while i < len(args):
        arg = args[i]
        if arg in ("--theme", "-t") and i + 1 < len(args):
            theme = args[i + 1].lower()
            i += 2
        elif arg.startswith("--theme="):
            theme = arg.split("=")[1].lower()
            i += 1
        elif not arg.startswith("-"):
            target = arg
            i += 1
        else:
            i += 1

    if not target:
        target = find_latest_results_dir(fmo_dir)

    if not target or not os.path.exists(target):
        print("Uso: python3 plot_benchmark.py [PATH_CARTELLA_RISULTATI_O_CSV] [--theme light|dark]")
        print("Nessuna directory dei risultati trovata.")
        sys.exit(1)

    print(f"Caricamento dati da: {target} (Tema grafici SVG: {theme})")
    summary_rows, raw_rows, stab_rows, results_dir = load_data(target)

    print_terminal_summary(summary_rows, stab_rows)
    generate_all_svg_charts(summary_rows, results_dir, stab_rows, theme=theme)
    generate_html_report(summary_rows, raw_rows, results_dir, stab_rows)
    generate_matplotlib_figures(summary_rows, results_dir)


if __name__ == "__main__":
    main()
