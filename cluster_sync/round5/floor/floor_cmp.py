import pandas as pd, numpy as np
R = '/scratch/nishanth.r/tordiff/results/'
for tag, old in (('floor', 'analysis/local_structure_test.csv'), ('floor_mmff', 'analysis/local_structure_test_mmff.csv')):
    N = pd.read_csv(R + f'round5/{tag}/local_structure_test_fixed.csv')
    O = pd.read_csv(R + old)
    err = N[N.get('error').notna()] if 'error' in N else N.iloc[:0]
    print(f'== {tag}: rows {len(N)}, errors {len(err)}:', err.error.str.split(':').str[0].value_counts().to_dict())
    N = N[N.get('error').isna()] if 'error' in N else N
    m = N.merge(O[['corrected_smiles', 'gt_idx', 'floor_gtother', 'n_gt']], on=['corrected_smiles', 'gt_idx'],
                suffixes=('', '_old'))
    big = m[m.n_gt > 3]
    d = (big.floor_gtother - big.floor_gtother_old).abs()
    print(f'G9 n_gt>3 rows matched {len(big)}: max |new-old| {d.max():.2e}, share equal (<1e-6) {(d < 1e-6).mean():.4f}')
    sm = m[m.n_gt <= 3]
    print(f'n_gt<=3 rows {len(sm)}: old floor_gtother mean {sm.floor_gtother_old.mean():.5f} (bug), fixed mean '
          f'{sm.floor_gtother.mean():.4f}, per-molecule {sm.groupby("corrected_smiles").floor_gtother.mean().mean():.4f}')
    for col in ('floor_gtother', 'floor_gtother_half', 'floor_gtother_half_all', 'floor_best_sym'):
        v = N[col].dropna()
        print(f'  {col}: rows {len(v)} mean {v.mean():.4f} per-molecule {N.dropna(subset=[col]).groupby("corrected_smiles")[col].mean().mean():.4f}')
