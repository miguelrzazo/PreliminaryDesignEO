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


def _mission(args: argparse.Namespace) -> None:
    from .mission import run_selected_mission

    artifact = run_selected_mission(
        args.config, args.output_dir, data_path=args.orekit_data,
        html=not args.no_html, backend=args.backend,
    )
    print(f"Mission analysis: {artifact}")


def _design(args: argparse.Namespace) -> None:
    import json
    from .design_optimizer import run_design
    artifact = run_design(args.output_dir, master_path=args.master,
        cloud_fraction=args.cloud_fraction, persistence=args.persistence,
        mc_samples=args.mc_samples, seed=args.seed, duration_days=args.duration_days,
        data_path=args.orekit_data, criterion=args.criterion)
    result = json.loads(artifact.read_text())
    print(f"Design search: {artifact}")
    print(f"Status: {result['status']}; hardware candidates: {result['hardware_candidates']}")
    print("Selected result reused: no")


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

    design = sub.add_parser("design", help="recompute and optimize the full MATLAB trade space with Orekit/Markov coverage")
    design.add_argument("--master", type=Path, help="authoritative MATLAB master.m inputs")
    design.add_argument("--output-dir", type=Path, default=Path("python/outputs/design"))
    design.add_argument("--orekit-data", type=Path)
    design.add_argument("--cloud-fraction", type=float, help="default: the same cloud fraction as MATLAB master.m")
    design.add_argument("--persistence", type=float, default=0.5)
    design.add_argument("--mc-samples", type=int, default=200)
    design.add_argument("--seed", type=int, default=12345)
    design.add_argument("--duration-days", type=float, default=30)
    design.add_argument("--criterion", choices=("matlab-mean", "p95"), default="matlab-mean")
    design.set_defaults(func=_design)

    mission = sub.add_parser("mission", help="propagate the selected thesis mission and compute timed accesses")
    mission.add_argument("--config", type=Path, default=_default_config().with_name("selected_mission.yaml"))
    mission.add_argument("--output-dir", type=Path, default=Path("python/outputs"))
    mission.add_argument("--orekit-data", type=Path, help="Orekit data archive/directory (or OREKIT_DATA_PATH)")
    mission.add_argument("--backend", choices=("orekit", "scipy-reference"), default="orekit")
    mission.add_argument("--no-html", action="store_true", help="export numerical artifacts without the explorer")
    mission.set_defaults(func=_mission)
    return parser


def main(argv: list[str] | None = None) -> None:
    parser = build_parser()
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
