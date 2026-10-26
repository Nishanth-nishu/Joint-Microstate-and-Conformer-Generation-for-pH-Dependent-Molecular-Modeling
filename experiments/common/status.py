"""GPU-h ledger and STATUS.md for Stage A (plan_P3 §3.6; round5/DECISION.md graft 2: "GPU-h ledger with a 25 GPU-h
stop, STATUS.md"). Sources: the ledger written by ``fx_ledger`` (one line per finished task row) and ``sacct`` for the
jobs named ``fx_*`` (their elapsed GPU time counts too, so failed and running jobs are not free), plus ``FX_FAIL`` /
``FX_OOM_RETRY`` lines in the logs.

``python experiments/common/status.py --res $RES --logs $LOGS [--planned GPU_H]`` writes ``$RES/round5/STATUS.md`` and
exits 4 (FX_BUDGET_STOP) if used + planned exceeds the cap.
"""
from __future__ import annotations

import subprocess
import sys
from argparse import ArgumentParser
from pathlib import Path

import pandas as pd

CAP_GPU_H = 25.0   # round5/DECISION.md item 10: "If still over 25 GPU-h, stop and ask the user"
STAGE_A_START = '2026-10-16'


def ledger(res: Path) -> pd.DataFrame:
    """Finished task rows: job, tag, start, end, GPU-h."""
    p = res / 'round5' / 'gpu_ledger.tsv'
    if not p.exists():
        return pd.DataFrame(columns=['job', 'tag', 't0', 't1', 'gpu_h'])
    return pd.read_csv(p, sep='\t', names=['job', 'tag', 't0', 't1', 'gpu_h'])


def slurm_gpu_h() -> pd.DataFrame:
    """GPU-h of every fx_* job since the Stage A start from sacct (elapsed × allocated GPUs)."""
    cmd = ['sacct', '-X', '-n', '-P', '-S', STAGE_A_START, '--name=fx_inf,fx_smoke,fx_m1,fx_cpu',
           '--format=JobID,JobName,State,ElapsedRaw,AllocTRES']
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return pd.DataFrame(columns=['job', 'name', 'state', 'gpu_h'])
    rows = []
    for line in out.splitlines():
        job, name, state, sec, tres = (line.split('|') + [''] * 5)[:5]
        ngpu = sum(int(t.split('=')[1]) for t in tres.split(',') if t.startswith('gres/gpu=') and '=' in t)
        rows.append(dict(job=job, name=name, state=state, gpu_h=int(sec or 0) * ngpu / 3600))
    return pd.DataFrame(rows)


def alerts(logs: Path) -> list[str]:
    """FX_FAIL / FX_OOM_RETRY / FX_PREFLIGHT failure lines of the fx_* job logs."""
    out = []
    for p in sorted(logs.glob('fx_*.log')):
        for line in p.read_text(errors='replace').splitlines():
            if line.startswith(('FX_FAIL', 'FX_OOM_RETRY')) or (line.startswith('FX_PREFLIGHT') and 'ok' not in line):
                out.append(f'{p.name}: {line}')
    return out


def main() -> None:
    ap = ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument('--res', required=True)
    ap.add_argument('--logs', required=True)
    ap.add_argument('--planned', type=float, default=0.0)
    a = ap.parse_args()
    res = Path(a.res)
    L, S = ledger(res), slurm_gpu_h()
    used = float(S.gpu_h.sum()) if len(S) else float(L.gpu_h.sum())
    lines = ['# Stage A STATUS', '', f'GPU-h used {used:.2f} of cap {CAP_GPU_H} (sacct, fx_* jobs); '
             f'ledger rows {len(L)} = {float(L.gpu_h.sum()):.2f} GPU-h; planned {a.planned:.2f}', '',
             '## Jobs', S.to_string(index=False) if len(S) else '(none)', '',
             '## Finished task rows', L[['tag', 'gpu_h']].to_string(index=False) if len(L) else '(none)', '',
             '## Alerts'] + (alerts(Path(a.logs)) or ['(none)'])
    stop = used + a.planned > CAP_GPU_H
    if stop:
        lines += ['', f'FX_BUDGET_STOP used {used:.2f} + planned {a.planned:.2f} > {CAP_GPU_H}: ask the user']
    (res / 'round5').mkdir(parents=True, exist_ok=True)
    (res / 'round5' / 'STATUS.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines[:3]))
    sys.exit(4 if stop else 0)


if __name__ == '__main__':
    main()
