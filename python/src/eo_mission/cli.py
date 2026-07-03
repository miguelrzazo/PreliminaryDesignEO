"""Command line interface for the EO mission portfolio pipeline."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _default_config() -> Path:
    candidates = (
        Path.cwd() / "python" / "config" / "portfolio.yaml",
        Path(__file__).parents[2] / "config" / "portfolio.yaml",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[-1]


def _run(args: argparse.Namespace) -> None:
    from .pipeline import run_pipeline

    result = run_pipeline(args.config, args.output_dir, initialize_orekit=args.orekit)
    print(f"Wrote artifacts to {result.output_dir}")
    print(f"Feasible solutions: {len(result.ranked)}")


def _plot(args: argparse.Namespace) -> None:
    from .artifacts import load_tensor
    from .plots import plot_summary_heatmaps

    out = Path(args.output_dir)
    snr = load_tensor(out / "tensors" / "snr.npz")
    mtf = load_tensor(out / "tensors" / "mtf.npz")
    coverage = load_tensor(out / "tensors" / "coverage.npz")
    ranked_path = out / "tables" / "ranked_solutions.csv"
    ranked = pd.read_csv(ranked_path) if ranked_path.exists() else pd.DataFrame()
    plot_summary_heatmaps(snr, mtf, coverage, ranked, out)
    print(f"Wrote heatmaps to {out / 'heatmaps'}")


def _export(args: argparse.Namespace) -> None:
    out = Path(args.output_dir)
    ranked = pd.read_csv(out / "tables" / "ranked_solutions.csv")
    if args.gsd is not None:
        ranked = ranked[ranked["gsd_m"] == args.gsd]
    ranked = ranked.head(args.limit)
    target = Path(args.path)
    target.parent.mkdir(parents=True, exist_ok=True)
    ranked.to_csv(target, index=False)
    print(f"Wrote {len(ranked)} ranked solutions to {target}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="EO mission portfolio pipeline")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="run tensors, ranking and plots")
    run.add_argument("--config", type=Path, default=_default_config())
    run.add_argument("--output-dir", type=Path, default=Path("output/python-orekit"))
    run.add_argument("--orekit", action="store_true", help="initialize Orekit JVM")
    run.set_defaults(func=_run)

    plot = sub.add_parser("plot", help="regenerate heatmaps from saved tensors")
    plot.add_argument("--output-dir", type=Path, default=Path("output/python-orekit"))
    plot.set_defaults(func=_plot)

    export = sub.add_parser("export", help="export top ranked feasible solutions")
    export.add_argument("--output-dir", type=Path, default=Path("output/python-orekit"))
    export.add_argument("--path", type=Path, default=Path("output/python-orekit/top_solutions.csv"))
    export.add_argument("--gsd", type=float)
    export.add_argument("--limit", type=int, default=50)
    export.set_defaults(func=_export)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
