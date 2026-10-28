# Smoke 4: exact commands (job 45272, gnode118, 2026-10-28 16:55:59–17:00:12 UTC / 22:25:59–22:30:12 IST)

Submitted from the staging checkout:

    cd /home2/nishanth.r/jmcg_stage && sbatch slurm_jmcg/cg_smoke4.sbatch

The job (`slurm_jmcg/cg_smoke4.sbatch`) runs on gnode118, CPU only, 4 cores, with `OMP_NUM_THREADS=4`:

    rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/        # P=/scratch/nishanth.r/jmcg (git HEAD df5210a)
    source /scratch/nishanth.r/tordiff/venv/bin/activate   # torch 1.13.1, rdkit 2022.09.5, pyg 2.0.4, e3nn 0.5.1
    export PYTHONPATH=$P TD_REPO=/scratch/nishanth.r/tordiff/torsional-diffusion CUDA_VISIBLE_DEVICES=""
    python -m jmcg.condgen.test_states
    python -m jmcg.condgen.inventory --workdir $TD/workdir --qm9_raw $TD/data/QM9/qm9/ --qm9_split $TD/data/QM9/split.npy --out $OUT/inventory.json
    QM9_DIR=$TD/workdir/pretrained/qm9_default             # sha256 6a9a6ee9b0e51bda... checked, job stops on mismatch
    python -m jmcg.condgen.test_backend --model_dir $QM9_DIR --states_csv jmcg/condgen/smoke_states.csv
    python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/qm9       --states_csv jmcg/condgen/smoke_states.csv --n_confs 8 --seed 0 --lock_arm
    python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/qm9_proc2 --states_csv jmcg/condgen/smoke_states.csv --n_confs 8 --seed 0 --lock_arm
    python -m jmcg.condgen.revalidate --run_dir $OUT/qm9 --out $OUT/qm9_reread
    python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/rep_a     --states_csv $OUT/repro_states.csv --n_confs 8 --seed 0
    python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/rep_b     --states_csv $OUT/repro_states.csv --n_confs 8 --seed 0
    python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/rep_seed1 --states_csv $OUT/repro_states.csv --n_confs 8 --seed 1

- `OUT` is `/scratch/nishanth.r/jmcg/results/cond_gen/smoke4`.
- Small files (json, csv, sdf, log) are copied to this directory.
- Inputs: `jmcg/condgen/smoke_states.csv`, sha256 db070079…, 24 rows.
- Seed 0, 8 conformers per method, 20 SDE steps, no force-field refinement.
- TD's torus import seed is 20221024, and the table hash is 8cf17268cc0b0dee (recorded in `qm9/run_info.json`).
- Full log: `/home2/nishanth.r/jmcg_stage/logs/jmcg_cgsmoke4_45272.log`.
