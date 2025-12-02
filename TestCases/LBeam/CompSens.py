#!/usr/bin/env python3
import os
import csv
import shutil
import argparse
from pathlib import Path

CALTOP = "calTop.exe"
FILTER_FILES = ["dcol.bin", "drow.bin", "dval.bin", "dnnz.bin", "dsum.bin"]

# ------------------------ Utilities ------------------------

def run_caltop(job, p, r, f, pexp, sigmin, sigrelax):
    cmd = (
        f"{CALTOP} -i {job} -p {p} -r {r} -f {f} "
        f"--pexp {pexp} --sigmin {sigmin} --sigrelax {sigrelax}"
    )
    rc = os.system(cmd)
    if rc != 0:
        raise RuntimeError(f"Command failed (rc={rc}): {cmd}")

def copy_filter_bins(target_dir: Path):
    """
    Copy the required filter binaries into target_dir, if not already present.
    Searches current working directory first, then its parent.
    """
    target_dir.mkdir(parents=True, exist_ok=True)
    search_dirs = [Path.cwd(), Path.cwd().parent]

    for fname in FILTER_FILES:
        dest = target_dir / fname
        if dest.exists():
            continue  # already there; skip to avoid SameFileError

        src = None
        for sd in search_dirs:
            candidate = sd / fname
            if candidate.exists():
                src = candidate
                break
        if src is None:
            raise FileNotFoundError(f"Missing required filter file: {fname}")
        shutil.copy(src, dest)

def read_density(path="density.dat"):
    with open(path, "r") as fh:
        return [float(line.strip()) for line in fh if line.strip()]

def write_density(dens, path="density.dat"):
    with open(path, "w") as fh:
        for v in dens:
            fh.write(f"{float(v)}\n")

def read_adjoint_stress(path="stress_sens.csv"):
    """Reads first column as adjoint sensitivity per element (skips header)."""
    vals = []
    with open(path, "r", newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd, None)  # assume a header row exists
        for row in rd:
            if not row:
                continue
            vals.append(float(row[0]))
    return vals

def read_pnorm_from_objectives(path="objectives.csv"):
    """Reads PNORM column from objectives.csv (expects header row)."""
    with open(path, "r", newline="") as fh:
        rows = list(csv.reader(fh))
    if len(rows) < 2:
        raise ValueError(f"Unexpected format in {path}")
    header = [h.strip() for h in rows[0]]
    vals   = [v.strip() for v in rows[1]]
    try:
        j = header.index("PNORM")
    except ValueError:
        raise ValueError("PNORM column not found in objectives.csv header")
    return float(vals[j])

# ------------------------ Core Routine ------------------------

def verify_stress_sens_cfd(job, p, r, f, pexp, sigmin, sigrelax,
                           h=1e-6, outdir="Sensitivity_stress_CFD",
                           max_elems=None, start_idx=0):
    """
    Central finite-difference verification of P-norm stress sensitivities.
    - Uses constant absolute perturbation h (default 1e-6).
    - Respects bounds: rho in [0,1]; falls back to asymmetric or one-sided if needed.
    - Writes result_stress_cfd.csv with [index, rho, adjoint, CFD, rel_error].
      rel_error = |CFD - adjoint| / max(1e-12, |adjoint|)
    """

    # 0) Run once in CWD to (re)generate density and required files
    run_caltop(job, p, r, f, pexp, sigmin, sigrelax)
    cwd = Path.cwd()
    inp_path = cwd / f"{job}.inp"
    if not inp_path.exists():
        raise FileNotFoundError(f"Missing input file: {inp_path}")
    density_cwd = read_density("density.dat")

    # 1) Create sandbox and copy required files once
    odir = Path(outdir)
    odir.mkdir(exist_ok=True)
    shutil.copy(inp_path, odir / f"{job}.inp")
    shutil.copy(cwd / "density.dat", odir / "density.dat")
    copy_filter_bins(odir)

    # 2) Work inside sandbox
    os.chdir(odir)

    # Baseline run (gets a consistent state and outputs)
    run_caltop(job, p, r, f, pexp, sigmin, sigrelax)
    adjoint = read_adjoint_stress("stress_sens.csv")
    density = read_density("density.dat")
    pnorm0  = read_pnorm_from_objectives("objectives.csv")  # baseline PNORM (used if one-sided)

    n = min(len(density), len(adjoint))
    if max_elems is not None:
        n = min(n, start_idx + max_elems)
    idxs = range(start_idx, n)

    CFD = []
    ABS_REL_ERR = []
    RATIO = []
    eps = 1e-12

    for i in idxs:
        rho = density[i]
        # Respect bounds for +/- h
        h_plus  = min(h, 1.0 - rho)
        h_minus = min(h, rho)

        if h_plus > 0.0 and h_minus > 0.0 and abs(h_plus - h_minus) < 1e-15:
            # Symmetric central difference
            dens_plus = density[:]
            dens_plus[i] = rho + h_plus
            write_density(dens_plus, "density.dat")
            run_caltop(job, p, r, f, pexp, sigmin, sigrelax)
            p_plus = read_pnorm_from_objectives("objectives.csv")

            dens_minus = density[:]
            dens_minus[i] = rho - h_minus
            write_density(dens_minus, "density.dat")
            run_caltop(job, p, r, f, pexp, sigmin, sigrelax)
            p_minus = read_pnorm_from_objectives("objectives.csv")
            cfd = (p_plus - p_minus) / (2* h)

        else:
            # rho exactly at bound and h=0 both ways
            cfd = 0.0

        # Restore baseline density for next element
        write_density(density, "density.dat")

        CFD.append(cfd)
        denom = max(eps, abs(adjoint[i]))
        rel_err = abs((cfd - adjoint[i]) / denom) * 100
        if cfd!=0:
            ratio = adjoint[i]/cfd
        else:
            ratio =999999
        ABS_REL_ERR.append(rel_err)
        RATIO.append(ratio)
        print(f"[{i}] rho={rho:.6f}  CFD={cfd:.6e}  adj={adjoint[i]:.6e}  rel_err={rel_err:.3e} ratio = {ratio:.6f}")

    # 3) Write comparison CSV (with relative error instead of simple difference)
    with open("result_stress_cfd.csv", "w", newline="") as fh:
        wr = csv.writer(fh)
        wr.writerow(["index", "rho", "adjoint", "CFD", "rel_error", "ratio"])
        for k, i in enumerate(idxs):
            wr.writerow([i, density[i], adjoint[i], CFD[k], ABS_REL_ERR[k], RATIO[k]])

    print("Wrote result_stress_cfd.csv")

# ------------------------ CLI ------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Central FD check for P-norm stress sensitivities (uses constant h=1e-6 by default)."
    )
    ap.add_argument("-i", "--input", default="LBeam", help="Job/input name without extension (default: SCB)")
    ap.add_argument("-p", "--penal", type=float, default=3.0)
    ap.add_argument("-r", "--radius", type=float, default=0.05)
    ap.add_argument("-f", "--filternnz", type=int, default=500)
    ap.add_argument("--pexp", type=int, default=3)
    ap.add_argument("--sigmin", type=float, default=100)
    ap.add_argument("--sigrelax", type=float, default=0.001)
    ap.add_argument("--h", type=float, default=1e-6, help="Constant absolute perturbation (default 1e-6)")
    ap.add_argument("--outdir", default="Sensitivity_stress_CFD")
    ap.add_argument("--max-elems", type=int, default=5000, help="Limit number of elements")
    ap.add_argument("--start-idx", type=int, default=10999, help="Start element index")
    args = ap.parse_args()

    verify_stress_sens_cfd(
        job=args.input, p=args.penal, r=args.radius, f=args.filternnz,
        pexp=args.pexp, sigmin=args.sigmin, sigrelax=args.sigrelax,
        h=args.h, outdir=args.outdir,
        max_elems=args.max_elems, start_idx=args.start_idx
    )

if __name__ == "__main__":
    main()
