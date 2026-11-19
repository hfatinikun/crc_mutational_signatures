#!/usr/bin/env python3
"""Python fitters: SigProfilerAssignment and MuSiCal."""

from __future__ import annotations

import shutil
import tempfile
from pathlib import Path
from typing import Tuple

import pandas as pd

from helpers_py import align_contexts


def fit_SigProfilerAssignment(catalogue: pd.DataFrame, cosmic: pd.DataFrame, genome_build: str = "GRCh38") -> Tuple[pd.DataFrame, pd.DataFrame]:
    from SigProfilerAssignment import Analyzer as Analyze

    cat, cos = align_contexts(catalogue, cosmic)
    tmp = Path(tempfile.mkdtemp(prefix="spa_"))
    try:
        matrix_path = tmp / "matrix.txt"
        # SPA expects MutationType as first column name sometimes; keep index name
        cat_out = cat.copy()
        cat_out.index.name = cat_out.index.name or "MutationType"
        cat_out.to_csv(matrix_path, sep="\t")

        cosmic_path = tmp / "cosmic.txt"
        cos_out = cos.copy()
        cos_out.index.name = cos_out.index.name or "Type"
        cos_out.to_csv(cosmic_path, sep="\t")

        out_dir = tmp / "assignment"
        Analyze.cosmic_fit(
            samples=str(matrix_path),
            output=str(out_dir),
            input_type="matrix",
            signature_database=str(cosmic_path),
            genome_build=genome_build,
            make_plots=False,
        )

        # Locate activities file
        activities = None
        for p in out_dir.rglob("Assignment_Solution_Activities.txt"):
            activities = p
            break
        if activities is None:
            for p in out_dir.rglob("*Activities*.txt"):
                activities = p
                break
        if activities is None:
            raise FileNotFoundError(f"SPA activities not found under {out_dir}")

        act = pd.read_csv(activities, sep="\t", index_col=0)
        # SPA activities: samples × signatures
        if act.shape[0] == cat.shape[1] and act.shape[1] == cos.shape[1]:
            contrib = act.T
        elif act.shape[1] == cat.shape[1] and act.shape[0] == cos.shape[1]:
            contrib = act
        else:
            # try transpose heuristics by name overlap
            if set(act.columns).issubset(set(cos.columns)) or set(act.columns) & set(cos.columns):
                contrib = act.T
            else:
                contrib = act

        # Align signatures present in COSMIC
        common_sigs = [s for s in cos.columns if s in contrib.index]
        if not common_sigs:
            # rename / intersect loosely
            common_sigs = list(contrib.index.intersection(cos.columns))
        if not common_sigs:
            raise ValueError("No overlapping signatures between SPA activities and COSMIC")
        contrib = contrib.loc[common_sigs]
        # Align samples
        common_smp = [s for s in cat.columns if s in contrib.columns]
        if not common_smp:
            # SPA may strip sample names oddly
            if contrib.shape[1] == cat.shape[1]:
                contrib.columns = cat.columns
                common_smp = list(cat.columns)
            else:
                raise ValueError("SPA sample names do not match catalogue")
        contrib = contrib.loc[:, common_smp]
        cos_use = cos.loc[:, common_sigs]
        recon = cos_use @ contrib
        recon.index = cos_use.index
        recon.columns = contrib.columns
        return contrib, recon
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def fit_MuSiCal(catalogue: pd.DataFrame, cosmic: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    import musical

    cat, cos = align_contexts(catalogue, cosmic)
    H, _model = musical.refit.refit(
        cat,
        cos,
        method="likelihood_bidirectional",
        thresh=0.001,
    )
    # H typically signatures × samples or samples × signatures
    if isinstance(H, pd.DataFrame):
        contrib = H
    else:
        contrib = pd.DataFrame(H, index=cos.columns, columns=cat.columns)

    if contrib.shape[0] == cat.shape[1] and contrib.shape[1] == cos.shape[1]:
        contrib = contrib.T
    if list(contrib.index) != list(cos.columns):
        # try set names
        if contrib.shape[0] == cos.shape[1]:
            contrib.index = cos.columns
        if contrib.shape[1] == cat.shape[1]:
            contrib.columns = cat.columns

    common_sigs = [s for s in cos.columns if s in contrib.index]
    contrib = contrib.loc[common_sigs]
    common_smp = [s for s in cat.columns if s in contrib.columns]
    if not common_smp and contrib.shape[1] == cat.shape[1]:
        contrib.columns = cat.columns
        common_smp = list(cat.columns)
    contrib = contrib.loc[:, common_smp]
    cos_use = cos.loc[:, common_sigs]
    recon = cos_use @ contrib
    recon.index = cos_use.index
    recon.columns = contrib.columns
    return contrib, recon


FITTERS_PY = {
    "SigProfilerAssignment": fit_SigProfilerAssignment,
    "MuSiCal": fit_MuSiCal,
}
