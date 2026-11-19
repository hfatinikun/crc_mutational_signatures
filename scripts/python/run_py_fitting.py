#!/usr/bin/env python3
"""Run Python fitting grid (SPA + MuSiCal).

Usage (from mywork/):
  python scripts/python/run_py_fitting.py
  python scripts/python/run_py_fitting.py --tools MuSiCal --cosmic v2,v3.2 --datasets Preclinical
"""

from __future__ import annotations

import argparse
import sys
import traceback
from datetime import datetime
from pathlib import Path

# Allow importing sibling modules when invoked as a script
sys.path.insert(0, str(Path(__file__).resolve().parent))

from fitters_py import FITTERS_PY  # noqa: E402
from helpers_py import (  # noqa: E402
    append_cosine_tsv,
    find_mywork_root,
    load_catalogue,
    load_cosmic,
    load_pipeline_config,
    save_fit_outputs,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tools", default=None, help="Comma-separated tool IDs")
    parser.add_argument("--cosmic", default=None, help="Comma-separated COSMIC keys")
    parser.add_argument("--datasets", default=None, help="Comma-separated dataset names")
    parser.add_argument("--phase", type=int, default=None, help="Only datasets tagged with this phase in config")
    args = parser.parse_args()

    root = find_mywork_root()
    cfg = load_pipeline_config(root)

    tools = [
        name
        for name, meta in cfg["tools"].items()
        if meta.get("enabled") and meta.get("language") == "Python" and name in FITTERS_PY
    ]
    if args.tools:
        tools = [t for t in args.tools.split(",") if t in tools]

    cosmic_keys = list(cfg["cosmic"].keys())
    if args.cosmic:
        cosmic_keys = [c for c in args.cosmic.split(",") if c in cosmic_keys]

    datasets = list(cfg["datasets"].keys())
    if args.phase is not None:
        datasets = [
            d for d in datasets
            if cfg["datasets"][d].get("phase") == args.phase
        ]
    if args.datasets:
        datasets = [d for d in args.datasets.split(",") if d in datasets]
    if not datasets:
        raise SystemExit("No datasets selected (check --phase / --datasets).")

    log_dir = root / "results" / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"py_fitting_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

    def log(msg: str) -> None:
        print(msg, flush=True)
        with open(log_path, "a") as f:
            f.write(msg + "\n")

    log(f"Tools: {tools}")
    log(f"COSMIC: {cosmic_keys}")
    log(f"Datasets: {datasets}")

    catalogues = {d: load_catalogue(d, cfg, root) for d in datasets}

    for tool in tools:
        fit_fun = FITTERS_PY[tool]
        log(f"\n==== Tool: {tool} ====")
        for ck in cosmic_keys:
            log(f"-- COSMIC {ck}")
            try:
                cosmic = load_cosmic(ck, cfg, root)
            except Exception as e:
                log(f"  ERROR loading COSMIC: {e}")
                continue
            for ds in datasets:
                log(f"  Dataset: {ds}")
                data_type = cfg["datasets"][ds]["data_type"]
                cat = catalogues[ds]
                try:
                    contrib, recon = fit_fun(cat, cosmic)
                    cos_df = save_fit_outputs(
                        contrib, recon, cat.loc[recon.index, recon.columns], tool, ck, ds, data_type, root
                    )
                    append_cosine_tsv(cos_df, tool, root)
                    med = float(cos_df["Cosine_Similarity"].median())
                    assigned = int((cos_df["Cosine_Similarity"] >= 0.9).sum())
                    log(f"    OK: median cosine={med:.4g} assigned(>=0.9)={assigned}/{len(cos_df)}")
                except Exception as e:
                    log(f"    FIT/SAVE ERROR: {e}")
                    log(traceback.format_exc())

    log(f"\nDone. Log: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
