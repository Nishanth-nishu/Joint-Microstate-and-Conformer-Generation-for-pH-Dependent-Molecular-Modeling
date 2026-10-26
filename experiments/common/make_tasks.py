"""Writes the Stage A task tables read by infer.sbatch (stageA_tasks.tsv, smoke_tasks.tsv). Estimates: one QM9 test
2K pass = Σ2K 27,462 conformers × 0.034 GPU-s ≈ 0.26 GPU-h [E-SAMPLE-053, -054]; pools scale with their conformer
count (seed-0 max(128, 4K) = 151,088 → 1.43; 4K = 54,924 → 0.52, EXP-5 §12.1/§12.6); PG carries a kernel overhead
(×1.8, measured by the GPU smoke: 0.050 vs 0.028 s per conformer); Exp 3 affected subsets ≈ 10 % (102 molecules)
and 21 % (211) of a 2K pass
(EXP-3 §9.1); DRUGS arms under the 4.5-5 GPU-h cap (EXP-3 §9.4), ≈ 0.75 each, measured first by the timing row.
Waves: w1 validation + TD test pools; (CPU: analysis.py validation decides selector and PG setting); w2 PG test pools,
Exp 3 QM9 arms; w3 DRUGS (after the download). Usage: python experiments/common/make_tasks.py
"""
from pathlib import Path

HERE = Path(__file__).resolve().parent
# KEYS of a CPU post-hoc row (F_flip) names its base's subset, so the evaluator is not cut to the first rows
COLS = ['TAG', 'DATASET', 'ARM', 'SEED', 'N', 'KEYS', 'BASE', 'POST', 'LIMIT', 'GPUH_EST', 'WAVE']
PG_SETTINGS = ['inv_recall', 'inv_precision', 'noninv_recall', 'noninv_precision']
EXP3_ARMS = [('K1', 'keys_double.txt', 0.03), ('K1_wrong', 'keys_double.txt', 0.03),
             ('K1_canonical', 'keys_double_or_triple.txt', 0.06), ('K2_inf', 'keys_double.txt', 0.03),
             ('F_filter', 'keys_double.txt', 0.035)]


def stage_a(limit=0, scale=1.0, drugs=True):
    rows = [('exp5/val/td_s0', 'qm9val', 'O', 0, '4K', '-', '-', 'eval,features,validity', limit, 0.52, 'w1')]
    rows += [(f'exp5/val/pg_{s}_s0', 'qm9val', f'PG:{s}', 0, '2K', '-', '-', 'eval', limit, 0.47, 'w1')
             for s in PG_SETTINGS]
    rows += [(f'exp5/test/td_s{s}', 'qm9test', 'O', s, 'pool' if s == 0 else '4K', '-', '-',
              'eval,features,validity', limit, 1.43 if s == 0 else 0.52, 'w1') for s in (0, 1, 2)]
    for s in (0, 1, 2):
        rows += [(f'exp5/test/pg2k_s{s}', 'qm9test', 'PG:CHOSEN', s, '2K', '-', '-', 'eval', limit, 0.47, 'w2'),
                 (f'exp5/test/pg4k_s{s}', 'qm9test', 'PG:CHOSEN', s, '4K', '-', '-', 'eval,features,validity',
                  limit, 0.94, 'w2'),
                 (f'exp3/qm9/O_s{s}', 'qm9test', 'O_prefix', s, '2K', '-', f'exp5/test/td_s{s}', 'eval,stereo', limit,
                  0.0, 'w2'),
                 (f'exp3/qm9/F_flip_s{s}', 'qm9test', 'F_flip', s, '2K', 'keys_double.txt', f'exp3/qm9/O_s{s}',
                  'eval,stereo',
                  limit, 0.0, 'w2')]
        rows += [(f'exp3/qm9/{a}_s{s}', 'qm9test', a, s, '2K', k, f'exp3/qm9/O_s{s}', 'eval,stereo', limit, g, 'w2')
                 for a, k, g in EXP3_ARMS]
    if drugs:
        rows += [('exp3/drugs/O_s0', 'drugstest', 'O', 0, '2K', 'keys_drugs_double_or_triple.txt', '-', 'eval,stereo',
                  limit, 0.8, 'w3')]
        rows += [(f'exp3/drugs/{a}_s0', 'drugstest', a, 0, '2K', k.replace('keys_', 'keys_drugs_'), 'exp3/drugs/O_s0',
                  'eval,stereo', limit, 0.75, 'w3') for a, k, _ in EXP3_ARMS if a != 'F_filter'] + \
                [('exp3/drugs/F_flip_s0', 'drugstest', 'F_flip', 0, '2K', 'keys_drugs_double_or_triple.txt',
                  'exp3/drugs/O_s0', 'eval,stereo', limit,
                  0.0, 'w3'),
                 ('exp3/drugs/F_filter_s0', 'drugstest', 'F_filter', 0, '2K', 'keys_drugs_double.txt',
                  'exp3/drugs/O_s0', 'eval,stereo', limit, 0.3, 'w3')]
    return [r[:9] + (round(r[9] * scale, 3),) + r[10:] for r in rows]


def smoke():
    """20 molecules per new arm (QM9 test; Exp 3 arms on the double-bond keys), one wave, ≤ 0.5 GPU-h in total."""
    rows = [r for r in stage_a(limit=20, scale=20 / 1000, drugs=False) if r[1] == 'qm9test' and r[3] == 0]
    rows = [tuple(x.replace('exp5/', 'smoke/exp5/').replace('exp3/', 'smoke/exp3/') if isinstance(x, str) else x
                  for x in r[:-1]) + ('smoke',) for r in rows]
    rows = [r[:2] + ('PG:noninv_recall',) + r[3:] if r[2] == 'PG:CHOSEN' else r for r in rows]
    # Exp 3 O on the same 20 double-bond molecules as the arms (the pool covers the first 20 csv rows instead)
    rows = [r[:2] + ('O', r[3], r[4], 'keys_double.txt', '-') + r[7:] if r[2] == 'O_prefix' else r for r in rows]
    # EXP-3 §9.1 condition 3: an arm run on unaffected molecules must equal O (here: the pool prefix) there
    rows.append(('smoke/exp3/qm9/O_unaffected_s0', 'qm9test', 'O', 0, '2K', 'keys_unaffected.txt', '-', 'eval', 5,
                 0.002, 'smoke'))
    rows.append(('smoke/exp3/qm9/K1_unaffected_s0', 'qm9test', 'K1', 0, '2K', 'keys_unaffected.txt',
                 'smoke/exp3/qm9/O_unaffected_s0', 'identity', 5, 0.002, 'smoke'))
    rows.append(('smoke/exp5/val/pg_inv_recall_s0', 'qm9val', 'PG:inv_recall', 0, '2K', '-', '-', 'eval', 20, 0.006,
                 'smoke'))
    return rows


def write(path, rows):
    with open(path, 'w', newline='\n') as f:
        f.write('\t'.join(COLS) + '\n')
        for r in rows:
            f.write('\t'.join(str(x) for x in r) + '\n')


if __name__ == '__main__':
    write(HERE / 'stageA_tasks.tsv', stage_a())
    write(HERE / 'smoke_tasks.tsv', smoke())
    print('wrote', HERE / 'stageA_tasks.tsv', HERE / 'smoke_tasks.tsv')
