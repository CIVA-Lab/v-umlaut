"""Reproduce Figure 4 of the paper: synthetic experiments without RANSAC.

Runs the v_umlaut binary once per data point (scene x noise model x sigma) in parallel,
saves the error quartiles to <out>.csv and draws the figure to <out>.pdf and <out>.png.
"""
import argparse
import csv
import itertools
import os
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parent
SCENES = {'general': 'Standard', 'wami': 'WAMI'}
NOISES = {'independent': '(a) Pixel noise: 5 independent correspondences, 2 views',
          'dependent': '(b) Pixel noise: 2 dependent correspondences along 2 lines, 1 view'}
SIGMAS = [0.0, 0.5, 1.0, 1.5, 2.0]
METHODS = {'PoseLib-5pt': '5pE', 'PoseLib-7pt': '7pF', 'V-Umlaut': '5pF-V-Umlaut'}
STYLES = {'5pE': ('#1f77b4', 's', '--'), '7pF': ('#2ca02c', '^', '-.'), '5pF-V-Umlaut': ('#d62728', 'o', '-')}


def run(binary, trials, scene, noise, sigma):
    noise_args = [f'--{noise}-point-noise', str(sigma)] if sigma > 0 else []   # sigma 0: the binary's noise-free default
    out = subprocess.run([binary, scene, str(trials), *noise_args], capture_output=True, text=True, check=True).stdout
    rows, method = [], None
    for line in out.splitlines():
        line = line.strip()
        if line.startswith('['):
            method = METHODS[line[1:-1]]
        elif 'error (deg)' in line:
            stats = dict(re.findall(r'(\w+)=(\S+)', line))
            rows.append({'scene': scene, 'noise': noise, 'sigma': sigma, 'method': method,
                         'error': line.split()[0].lower(), 'q25': float(stats['Q25']),
                         'median': float(stats['median']), 'q75': float(stats['Q75'])})
    return rows


def plot(rows, out):
    plt.rcParams.update({'font.size': 21, 'axes.titlesize': 21, 'axes.titleweight': 'bold', 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, 4, figsize=(24.5, 10.9))
    fig.subplots_adjust(left=0.038, right=0.998, bottom=0.13, top=0.86, wspace=0.18, hspace=0.24)
    for col, (noise, scene) in enumerate(itertools.product(NOISES, SCENES)):
        for row, error in enumerate(['rotation', 'translation']):
            ax = axes[row, col]
            # Noise on the dependent points does not affect 5pE and 7pF, so (b) shows V-Umlaut only.
            for method in STYLES if noise == 'independent' else ['5pF-V-Umlaut']:
                color, marker, linestyle = STYLES[method]
                pts = sorted((r['sigma'], r['q25'], r['median'], r['q75']) for r in rows
                             if (r['scene'], r['noise'], r['method'], r['error']) == (scene, noise, method, error))
                sigma, q25, median, q75 = zip(*pts)
                ax.fill_between(sigma, q25, q75, color=color, alpha=0.12, linewidth=0)
                ax.plot(sigma, median, color=color, marker=marker, linestyle=linestyle, linewidth=2)
            ax.set(xlim=(-0.1, 2.1), ylim=(-3.5, 75), xticks=SIGMAS, yticks=[0, 20, 40, 60],
                   xlabel='Gaussian noise σ (pixels)', ylabel=f'Median {error.title()} Error ($^\\circ$)')
            ax.grid(linestyle='--', linewidth=0.4, alpha=0.3)
            if row == 0:
                ax.set_title(SCENES[scene])
    for i, caption in enumerate(NOISES.values()):
        left, right = axes[1, 2 * i].get_position(), axes[1, 2 * i + 1].get_position()
        fig.text((left.x0 + right.x1) / 2, 0.01, caption, ha='center', family='STIXGeneral', fontsize=30)
    handles = [Line2D([], [], color=c, marker=m, linestyle=ls, linewidth=2.6, markersize=18)
               for c, m, ls in STYLES.values()]
    fig.legend(handles, STYLES, loc='upper center', ncol=3, frameon=False, fontsize=28, handlelength=2.8)
    for ext in ['pdf', 'png']:
        fig.savefig(f'{out}.{ext}', bbox_inches='tight', pad_inches=0.05, dpi=150)


def main():
    exe = 'Release/v_umlaut.exe' if os.name == 'nt' else 'v_umlaut'
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--binary', default=str(ROOT / 'build' / exe), help='v_umlaut binary (default: %(default)s)')
    parser.add_argument('--trials', type=int, default=10000, help='trials per data point (default: %(default)s)')
    parser.add_argument('--out', default=str(ROOT / 'build' / 'figure4'), help='output prefix (default: %(default)s)')
    args = parser.parse_args()

    jobs = list(itertools.product(SCENES, NOISES, SIGMAS))
    with ThreadPoolExecutor(os.cpu_count()) as pool:
        rows = sum(pool.map(lambda job: run(args.binary, args.trials, *job), jobs), [])
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(f'{args.out}.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    plot(rows, args.out)
    print(f'Saved {args.out}.csv, {args.out}.pdf and {args.out}.png')


if __name__ == '__main__':
    main()
