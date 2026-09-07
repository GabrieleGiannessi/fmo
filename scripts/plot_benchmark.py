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


def generate_svg_chart(title, x_label, y_label, series_data, out_path, is_speedup=False, is_efficiency=False):
    """
    Genera un grafico vettoriale SVG autosufficiente e moderno senza dipendenze esterne.
    """
    w, h = 760, 460
    pad_left, pad_right, pad_top, pad_bottom = 80, 190, 60, 60
    plot_w = w - pad_left - pad_right
    plot_h = h - pad_top - pad_bottom

    # Raccogli tutti i punti X e Y
    all_x = []
    all_y = []
    for label, points in series_data.items():
        for px, py in points:
            if px is not None and py is not None:
                all_x.append(px)
                all_y.append(py)

    if not all_x or not all_y:
        return

    min_x, max_x = min(all_x), max(all_x)
    min_y = 0.0
    max_y = max(all_y)

    if is_speedup:
        max_y = max(max_y, max_x)
    elif is_efficiency:
        max_y = max(max_y, 105.0)

    if max_x == min_x:
        max_x = min_x + 1
    if max_y == min_y:
        max_y = min_y + 1

    # Funzioni coordinate pixel
    def to_px(val_x):
        return pad_left + (val_x - min_x) / (max_x - min_x) * plot_w

    def to_py(val_y):
        return pad_top + plot_h - (val_y - min_y) / (max_y - min_y) * plot_h

    # Colori serie
    palette = ['#58a6ff', '#3fb950', '#f0883e', '#d29922', '#bc8cff', '#79c0ff']

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}">',
        f'<rect width="{w}" height="{h}" fill="#0d1117" rx="8" />',
        # Titolo
        f'<text x="{pad_left}" y="36" fill="#f0f6fc" font-family="sans-serif" font-size="18" font-weight="bold">{title}</text>',
        # Riquadro area plot
        f'<rect x="{pad_left}" y="{pad_top}" width="{plot_w}" height="{plot_h}" fill="#161b22" stroke="#30363d" stroke-width="1" />'
    ]

    # Linee griglia orizzontale Y
    y_ticks = 5
    for i in range(y_ticks + 1):
        y_val = min_y + (max_y - min_y) * (i / y_ticks)
        py = to_py(y_val)
        svg.append(f'<line x1="{pad_left}" y1="{py}" x2="{pad_left + plot_w}" y2="{py}" stroke="#21262d" stroke-width="1" />')
        y_text = f"{y_val:.1f}" if is_speedup or is_efficiency else (f"{y_val:.0f}" if y_val >= 10 else f"{y_val:.2f}")
        svg.append(f'<text x="{pad_left - 10}" y="{py + 4}" fill="#8b949e" font-family="sans-serif" font-size="12" text-anchor="end">{y_text}</text>')

    # Linea ideale
    if is_speedup:
        svg.append(f'<line x1="{to_px(min_x)}" y1="{to_py(min_x)}" x2="{to_px(max_x)}" y2="{to_py(max_x)}" stroke="#8b949e" stroke-width="1.5" stroke-dasharray="5,5" />')
    elif is_efficiency:
        svg.append(f'<line x1="{pad_left}" y1="{to_py(100.0)}" x2="{pad_left + plot_w}" y2="{to_py(100.0)}" stroke="#8b949e" stroke-width="1.5" stroke-dasharray="5,5" />')

    # Tracciamento serie
    legend_y = pad_top + 15
    for idx, (label, points) in enumerate(series_data.items()):
        color = palette[idx % len(palette)]
        valid_pts = [(px, py) for px, py in points if px is not None and py is not None]
        valid_pts.sort(key=lambda p: p[0])

        if len(valid_pts) > 1:
            poly_points = " ".join([f"{to_px(px):.1f},{to_py(py):.1f}" for px, py in valid_pts])
            svg.append(f'<polyline fill="none" stroke="{color}" stroke-width="2.5" points="{poly_points}" />')

        # Marcatori
        for px, py in valid_pts:
            cx, cy = to_px(px), to_py(py)
            svg.append(f'<circle cx="{cx:.1f}" cy="{cy:.1f}" r="4" fill="{color}" stroke="#0d1117" stroke-width="1.5" />')

        # Legenda
        svg.append(f'<rect x="{pad_left + plot_w + 15}" y="{legend_y - 10}" width="12" height="12" fill="{color}" rx="2" />')
        svg.append(f'<text x="{pad_left + plot_w + 35}" y="{legend_y}" fill="#c9d1d9" font-family="sans-serif" font-size="12">{label}</text>')
        legend_y += 24

    if is_speedup:
        svg.append(f'<line x1="{pad_left + plot_w + 15}" y1="{legend_y - 4}" x2="{pad_left + plot_w + 27}" y2="{legend_y - 4}" stroke="#8b949e" stroke-dasharray="3,3" stroke-width="2" />')
        svg.append(f'<text x="{pad_left + plot_w + 35}" y="{legend_y}" fill="#8b949e" font-family="sans-serif" font-size="12">Ideale (S = p)</text>')
    elif is_efficiency:
        svg.append(f'<line x1="{pad_left + plot_w + 15}" y1="{legend_y - 4}" x2="{pad_left + plot_w + 27}" y2="{legend_y - 4}" stroke="#8b949e" stroke-dasharray="3,3" stroke-width="2" />')
        svg.append(f'<text x="{pad_left + plot_w + 35}" y="{legend_y}" fill="#8b949e" font-family="sans-serif" font-size="12">Ideale (100%)</text>')

    # Etichette assi
    svg.append(f'<text x="{pad_left + plot_w / 2}" y="{h - 15}" fill="#8b949e" font-family="sans-serif" font-size="13" text-anchor="middle">{x_label}</text>')
    svg.append(f'<text x="22" y="{pad_top + plot_h / 2}" fill="#8b949e" font-family="sans-serif" font-size="13" text-anchor="middle" transform="rotate(-90 22 {pad_top + plot_h / 2})">{y_label}</text>')

    # Tick X
    unique_x = sorted(list(set(all_x)))
    for ux in unique_x:
        px = to_px(ux)
        svg.append(f'<line x1="{px}" y1="{pad_top + plot_h}" x2="{px}" y2="{pad_top + plot_h + 5}" stroke="#8b949e" stroke-width="1" />')
        svg.append(f'<text x="{px}" y="{pad_top + plot_h + 20}" fill="#8b949e" font-family="sans-serif" font-size="11" text-anchor="middle">{ux}</text>')

    svg.append('</svg>')

    with open(out_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(svg))

    print(f"Grafico SVG salvato: {out_path}")


def generate_all_svg_charts(summary_rows, results_dir, stab_rows=None):
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

    generate_svg_chart(
        f"Speedup vs Cores ({target_phase})",
        "Numero di Cores / Workers",
        "Speedup (T1 / Tp)",
        speedup_series,
        os.path.join(results_dir, "speedup.svg"),
        is_speedup=True
    )

    generate_svg_chart(
        f"Efficienza Parallela vs Cores ({target_phase})",
        "Numero di Cores / Workers",
        "Efficienza Parallela (%)",
        eff_series,
        os.path.join(results_dir, "efficiency.svg"),
        is_efficiency=True
    )

    generate_svg_chart(
        f"Tempo di Esecuzione vs Cores ({target_phase})",
        "Numero di Cores / Workers",
        "Tempo (secondi)",
        time_series,
        os.path.join(results_dir, "execution_time.svg"),
        is_speedup=False
    )

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
                is_speedup=False
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
            <strong>Grafici e Dati:</strong>
            <a href="speedup.svg" target="_blank">📈 speedup.svg</a>
            <a href="efficiency.svg" target="_blank">⚡ efficiency.svg</a>
            <a href="execution_time.svg" target="_blank">⏱️ execution_time.svg</a>
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

    target = sys.argv[1] if len(sys.argv) > 1 else find_latest_results_dir(fmo_dir)

    if not target or not os.path.exists(target):
        print("Uso: python3 plot_benchmark.py <PATH_CARTELLA_RISULTATI_O_CSV>")
        print("Nessuna directory dei risultati trovata.")
        sys.exit(1)

    print(f"Caricamento dati da: {target}")
    summary_rows, raw_rows, stab_rows, results_dir = load_data(target)

    print_terminal_summary(summary_rows, stab_rows)
    generate_all_svg_charts(summary_rows, results_dir, stab_rows)
    generate_html_report(summary_rows, raw_rows, results_dir, stab_rows)
    generate_matplotlib_figures(summary_rows, results_dir)


if __name__ == "__main__":
    main()
