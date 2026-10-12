#!/usr/bin/env python3
"""Score ROM bands against a brute-force Monte Carlo dumped by
test_rom_coverage_structures's OneMemberDump mode (H17).

Usage: score_h17_mc.py DUMP_DIR ROM_TAG [ROM_TAG ...]
DUMP_DIR holds mc0.txt..mcN.txt (heads) and <ROM_TAG>.txt (heads, q05, q95, flag).
Scoring mirrors tests/regression/mc_quantiles.hpp: midpoint (Hazen) quantiles,
member coverage = fraction of members inside [q05, q95], samples with t > 60 s
and an MC 5-95 width >= 1e-3 ft (unit-agnostic: heads are internal feet).
"""
import glob, os, sys
import numpy as np

FLAGS = {0: "steady", 1: "no_flow", 2: "transient"}

def load(path):
    with open(path) as f:
        times = np.array([float(x) for x in f.readline().split()])
        data = {}
        for line in f:
            parts = line.split()
            data[(parts[0], parts[1])] = np.array([float(x) for x in parts[2:]])
    return times, data

def main():
    d, tags = sys.argv[1], sys.argv[2:]
    mcs = sorted(glob.glob(os.path.join(d, "mc*.txt")), key=lambda p: int(os.path.basename(p)[2:-4]))
    mc = [load(p) for p in mcs]
    nodes = sorted({k[0] for k in mc[0][1]})
    n_t = min(len(m[0]) for m in mc)
    print(f"{len(mcs)} members, {len(nodes)} nodes, {n_t} report times")
    H = np.stack([np.stack([m[1][(nm, 'h')][:n_t] for nm in nodes]) for m in mc])  # M x N x T
    Hs = np.sort(H, axis=0)
    q05 = np.quantile(Hs, 0.05, axis=0, method="hazen")
    q95 = np.quantile(Hs, 0.95, axis=0, method="hazen")
    mcw = q95 - q05
    times = mc[0][0][:n_t]
    for tag in tags:
        tr, rd = load(os.path.join(d, tag + ".txt"))
        n = min(n_t, len(tr))
        R05 = np.stack([rd[(nm, 'q05')][:n] for nm in nodes])
        R95 = np.stack([rd[(nm, 'q95')][:n] for nm in nodes])
        FL = np.stack([rd[(nm, 'flag')][:n] for nm in nodes]).astype(int)
        mask = (times[None, :n] > 60.0) & (mcw[:, :n] >= 1e-3)
        inside = ((H[:, :, :n] >= R05[None]) & (H[:, :, :n] <= R95[None])).mean(axis=0)
        ratio = (R95 - R05) / np.where(mcw[:, :n] > 0, mcw[:, :n], 1.0)
        def line(label, m):
            if m.sum() == 0:
                print(f"  {label:<26} samples=     0"); return
            print(f"  {label:<26} samples={m.sum():6d} member-coverage={inside[m].mean():.3f} "
                  f"width-med={np.median(ratio[m]):.3f}")
        print(f"[{tag}]")
        line("all", mask)
        line("material (MC >= 0.05 ft)", mask & (mcw[:, :n] >= 0.05))
        for f, nm in FLAGS.items():
            line("flag " + nm, mask & (FL == f))
        # per-node worst offenders by width (material samples)
        mm = mask & (mcw[:, :n] >= 0.05)
        rows = []
        for i, nm in enumerate(nodes):
            if mm[i].sum() >= 30:
                rows.append((np.median(ratio[i][mm[i]]), inside[i][mm[i]].mean(), nm, mm[i].sum()))
        rows.sort()
        if rows:
            print("  narrowest nodes:", ", ".join(f"{r[2]} {r[0]:.2f}/{r[1]:.2f}" for r in rows[:5]))
            print("  widest nodes:   ", ", ".join(f"{r[2]} {r[0]:.2f}/{r[1]:.2f}" for r in rows[-5:]))

if __name__ == "__main__":
    main()
