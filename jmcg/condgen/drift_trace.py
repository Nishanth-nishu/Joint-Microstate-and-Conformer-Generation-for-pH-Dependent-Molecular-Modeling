"""Cross-process drift diagnosis: trace one TDGenerator.generate call stage by stage, so two processes can be diffed.

Records (all as numpy): the ETKDG seed coordinates, the conformer positions after perturb_seeds, and for every model
call the input positions, the model output and a hash of the numpy / torch / python RNG states before the call; then
the final coordinates. `compare` reports the first stage where two traces differ and by how much.

Run: python -m jmcg.condgen.drift_trace trace --model_dir <dir> --state 'C/C(O)=C/C(C)=O' --seed 0 --out a.pkl [--repeat 2]
     python -m jmcg.condgen.drift_trace compare a.pkl b.pkl
--repeat 2 traces the same call twice in one process (trace k -> <out>.k).
"""
import argparse
import hashlib
import os
import pickle
import random

import numpy as np
import torch


def _h(x):
    return hashlib.sha1(np.ascontiguousarray(x).tobytes()).hexdigest()[:12]


def _rng_hash():
    return dict(np=_h(np.random.get_state()[1]), torch=_h(torch.random.get_rng_state().numpy()),
                py=hashlib.sha1(repr(random.getstate()).encode()).hexdigest()[:12])


def _np(x):
    if torch.is_tensor(x):
        return x.detach().cpu().double().numpy().copy()
    return None


def trace_once(gen, mol, n_confs, seed, hook_modules=False):
    from jmcg.condgen import td_backend as tb
    td = tb.td_sampling
    rec = dict(model_calls=[])
    orig_embed, orig_perturb, orig_fwd = td.embed_seeds, td.perturb_seeds, gen.model.forward

    def embed_seeds(*a, **k):
        confs, pdb = orig_embed(*a, **k)
        rec['seed_pos'] = [_np(c.pos) for c in confs]
        return confs, pdb

    def perturb_seeds(confs, *a, **k):
        out = orig_perturb(confs, *a, **k)
        rec['perturbed_pos'] = [_np(c.pos) for c in out]
        return out

    hooks, layer_out = [], []
    if hook_modules:  # per-module outputs, first model call only, in execution order
        def mk(name):
            def hook(mod, inp, out):
                if len(rec['model_calls']) == 0:
                    ts = [out] if torch.is_tensor(out) else [t for t in (out if isinstance(out, (tuple, list)) else [])
                                                             if torch.is_tensor(t)]
                    layer_out.append((name, type(mod).__name__, [_np(t) for t in ts],
                                      [_np(t) for t in inp if torch.is_tensor(t)]))
            return hook
        hooks = [m.register_forward_hook(mk(n)) for n, m in gen.model.named_modules()
                 if n and not isinstance(m, torch.jit.ScriptModule)]  # e3nn compiles its own submodules
    rec['layer_out'] = layer_out

    def forward(data, *a, **k):
        entry = dict(rng=_rng_hash(), pos_in=_np(data.pos))
        for name in ('node_sigma', 'edge_mask'):
            if hasattr(data, name):
                entry[name] = _np(getattr(data, name))
        out = orig_fwd(data, *a, **k)
        if len(rec['model_calls']) == 0:
            entry['data_after'] = {k_: _np(v) for k_, v in data.to_dict().items() if torch.is_tensor(v)}
        res = out if torch.is_tensor(out) else getattr(out, 'edge_pred', None)
        entry['out'] = _np(res)
        rec['model_calls'].append(entry)
        return out

    td.embed_seeds, td.perturb_seeds, gen.model.forward = embed_seeds, perturb_seeds, forward
    try:
        mols, _ = gen.generate(mol, n_confs, seed)
    finally:
        td.embed_seeds, td.perturb_seeds, gen.model.forward = orig_embed, orig_perturb, orig_fwd
        for h in hooks:
            h.remove()
    rec['final_pos'] = [m.GetConformer().GetPositions() for m in mols]
    return rec


def cmd_trace(a):
    from jmcg.condgen import states as st
    from jmcg.condgen.td_backend import TDGenerator
    gen = TDGenerator(a.model_dir, device='cpu')
    mol, _ = st.build_state(a.state)
    env = {k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'MKL_CBWR', 'PYTHONHASHSEED',
                                           'ATEN_CPU_CAPABILITY', 'DNNL_MAX_CPU_ISA')}
    meta = dict(env=env, torch_threads=torch.get_num_threads(), torch_config=torch.__config__.show(),
                deterministic=torch.are_deterministic_algorithms_enabled())
    for k in range(a.repeat):
        rec = trace_once(gen, mol, a.n_confs, a.seed, a.hook_modules)
        rec['meta'] = meta
        path = a.out if a.repeat == 1 else f'{a.out}.{k}'
        pickle.dump(rec, open(path, 'wb'))
        print(f'trace {path}: {len(rec["model_calls"])} model calls, env {env}, torch threads {meta["torch_threads"]}')


def _maxdiff(x, y):
    if x is None or y is None:
        return None if x is y else float('inf')
    x, y = np.asarray(x), np.asarray(y)
    return float('inf') if x.shape != y.shape else float(np.abs(x - y).max())


def cmd_compare(a):
    A, B = pickle.load(open(a.a, 'rb')), pickle.load(open(a.b, 'rb'))
    first = None
    rows = [('seed_pos', max(_maxdiff(x, y) for x, y in zip(A['seed_pos'], B['seed_pos']))),
            ('perturbed_pos', max(_maxdiff(x, y) for x, y in zip(A['perturbed_pos'], B['perturbed_pos'])))]
    for i, (p, q) in enumerate(zip(A['model_calls'], B['model_calls'])):
        rows.append((f'call{i:02d}.rng_same', p['rng'] == q['rng']))
        rows.append((f'call{i:02d}.pos_in', _maxdiff(p['pos_in'], q['pos_in'])))
        rows.append((f'call{i:02d}.out', _maxdiff(p['out'], q['out'])))
    rows.append(('final_pos', max(_maxdiff(x, y) for x, y in zip(A['final_pos'], B['final_pos']))))
    for name, v in rows:
        diverged = (v is False) or (isinstance(v, float) and v > 0)
        if diverged and first is None:
            first = name
        if a.verbose or name.startswith(('seed', 'pert', 'final')) or (diverged and name == first):
            print(f'{name:22s} {v}')
    outs = [v for n, v in rows if n.endswith('.out') and isinstance(v, float)]
    print(f'n_model_calls {len(A["model_calls"])}/{len(B["model_calls"])}  max|out diff| per call: '
          + ' '.join(f'{v:.1e}' for v in outs))
    print('FIRST_DIVERGENCE', first or 'none (bit-identical)')
    if A.get('data_after') or A['model_calls'] and 'data_after' in A['model_calls'][0]:
        da, db = A['model_calls'][0]['data_after'], B['model_calls'][0]['data_after']
        for k in sorted(set(da) | set(db)):
            print(f'data_after_call0.{k:20s} shape {getattr(da.get(k), "shape", None)} '
                  f'vs {getattr(db.get(k), "shape", None)} maxdiff {_maxdiff(da.get(k), db.get(k))}')
    if A.get('layer_out'):
        n_shown = 0
        for (na, ty, oa, ia), (nb, _, ob, ib) in zip(A['layer_out'], B['layer_out']):
            d_out = [_maxdiff(x, y) for x, y in zip(oa, ob)]
            d_in = [_maxdiff(x, y) for x, y in zip(ia, ib)]
            if na != nb or any(v for v in d_out + d_in) or a.verbose:
                print(f'LAYER {na:45s} {ty:24s} in {["%.1e" % v for v in d_in]} out {["%.1e" % v for v in d_out]}'
                      + ('' if na == nb else f' (order differs: {nb})'))
                n_shown += 1
                if n_shown >= 25 and not a.verbose:
                    break
        print(f'layers traced {len(A["layer_out"])}/{len(B["layer_out"])}')


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest='cmd', required=True)
    s = sp.add_parser('trace')
    s.add_argument('--model_dir', required=True); s.add_argument('--state', required=True)
    s.add_argument('--seed', type=int, default=0); s.add_argument('--n_confs', type=int, default=8)
    s.add_argument('--repeat', type=int, default=1); s.add_argument('--out', required=True)
    s.add_argument('--hook_modules', action='store_true')
    s = sp.add_parser('compare')
    s.add_argument('a'); s.add_argument('b'); s.add_argument('--verbose', action='store_true')
    a = ap.parse_args()
    dict(trace=cmd_trace, compare=cmd_compare)[a.cmd](a)


if __name__ == '__main__':
    main()
