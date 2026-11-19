#!/usr/bin/env python3
"""Shared helpers for Python mutational-signature fitting."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import numpy as np
import pandas as pd


def find_mywork_root(start: Optional[Path] = None) -> Path:
    start = (start or Path.cwd()).resolve()
    for cand in [start, *start.parents]:
        if (cand / "config" / "pipeline.json").exists():
            return cand
        my = cand / "mywork"
        if (my / "config" / "pipeline.json").exists():
            return my
    raise FileNotFoundError("Cannot locate mywork/ (config/pipeline.json). Run from mywork/.")


def load_pipeline_config(root: Optional[Path] = None) -> dict:
    root = root or find_mywork_root()
    with open(root / "config" / "pipeline.json") as f:
        return json.load(f)


def load_catalogue(dataset: str, cfg: dict, root: Path) -> pd.DataFrame:
    rel = cfg["datasets"][dataset]["file"]
    path = root / cfg["paths"]["input_dir"] / rel
    return pd.read_csv(path, sep="\t", index_col=0)


def load_cosmic(cosmic_key: str, cfg: dict, root: Path) -> pd.DataFrame:
    rel = cfg["cosmic"][cosmic_key]
    path = root / cfg["paths"]["cosmic_dir"] / rel
    return pd.read_csv(path, sep="\t", index_col=0)


def align_contexts(catalogue: pd.DataFrame, cosmic: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    common = [c for c in cosmic.index if c in catalogue.index]
    if not common:
        raise ValueError("No overlapping mutation contexts")
    return catalogue.loc[common], cosmic.loc[common]


def normalize_contributions(contrib: pd.DataFrame) -> pd.DataFrame:
    colsum = contrib.sum(axis=0).replace(0, 1.0)
    return contrib.div(colsum, axis=1)


def cosine_per_sample(observed: pd.DataFrame, reconstructed: pd.DataFrame) -> pd.Series:
    common_ctx = observed.index.intersection(reconstructed.index)
    common_smp = observed.columns.intersection(reconstructed.columns)
    X = observed.loc[common_ctx, common_smp].to_numpy(dtype=float)
    Y = reconstructed.loc[common_ctx, common_smp].to_numpy(dtype=float)
    out = {}
    for j, smp in enumerate(common_smp):
        x = X[:, j]
        y = Y[:, j]
        nx = np.linalg.norm(x)
        ny = np.linalg.norm(y)
        out[smp] = np.nan if nx == 0 or ny == 0 else float(np.dot(x, y) / (nx * ny))
    return pd.Series(out, name="Cosine_Similarity")


def write_matrix(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t")


def save_fit_outputs(
    contribution: pd.DataFrame,
    reconstructed: pd.DataFrame,
    catalogue: pd.DataFrame,
    tool: str,
    cosmic_key: str,
    dataset: str,
    data_type: str,
    root: Path,
) -> pd.DataFrame:
    # contribution: signatures × samples; reconstructed/catalogue: contexts × samples
    contrib_dir = root / "results" / "contributions" / tool / cosmic_key
    recon_dir = root / "results" / "reconstructed" / tool / cosmic_key
    write_matrix(contribution, contrib_dir / f"{dataset}_contribution.tsv")
    write_matrix(normalize_contributions(contribution), contrib_dir / f"{dataset}_contribution_normalized.tsv")
    write_matrix(reconstructed, recon_dir / f"{dataset}_reconstructed.tsv")

    cos = cosine_per_sample(catalogue, reconstructed)
    mut_n = catalogue.sum(axis=0)
    return pd.DataFrame(
        {
            "Sample": cos.index,
            "Cosine_Similarity": cos.values,
            "Dataset": dataset,
            "Data_Type": data_type,
            "Tool": tool,
            "COSMIC_Version": cosmic_key,
            "Mutation_Count": mut_n.loc[cos.index].values,
        }
    )


def append_cosine_tsv(df: pd.DataFrame, tool: str, root: Path) -> Path:
    out_dir = root / "results" / "cosine"
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"CosineResults_{tool}_ALL.tsv"
    if path.exists():
        old = pd.read_csv(path, sep="\t")
        key_new = set(zip(df["Dataset"], df["Data_Type"], df["COSMIC_Version"]))
        mask = ~old.apply(
            lambda r: (r["Dataset"], r["Data_Type"], r["COSMIC_Version"]) in key_new,
            axis=1,
        )
        old = old.loc[mask]
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(path, sep="\t", index=False)
    return path
