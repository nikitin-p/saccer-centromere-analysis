#!/usr/bin/env python3
import argparse
import subprocess
import sys
from pathlib import Path


STEPS = [
    "01_preprocess.py",
    "02_filter.py",
    "03_align.py",
    "04_analyze.py",
    "05_visualize.py",
]


def main():
    parser = argparse.ArgumentParser(description="Run full analysis pipeline")
    parser.add_argument("--config", required=True, help="Path to config.yaml")
    parser.add_argument(
        "--start-from",
        type=int,
        default=1,
        help="Start from step N (1-5)",
    )
    args = parser.parse_args()

    scripts_dir = Path(__file__).parent

    for i, script in enumerate(STEPS, start=1):
        if i < args.start_from:
            print(f"Skipping step {i}: {script}")
            continue

        print(f"\n{'='*60}")
        print(f"Step {i}: {script}")
        print(f"{'='*60}\n")

        result = subprocess.run(
            [sys.executable, str(scripts_dir / script), "--config", args.config],
            check=False,
        )

        if result.returncode != 0:
            print(f"\nError: Step {i} failed with exit code {result.returncode}")
            print(f"Fix the issue and re-run with --start-from {i}")
            sys.exit(result.returncode)

    print(f"\n{'='*60}")
    print("Pipeline complete!")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
