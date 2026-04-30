"""
SNPic Python CLI & API wrapper.

Provides a Pythonic interface to the SNPic R pipeline.
Requires R >= 4.5 with the required packages installed.
"""

import subprocess
import os
import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, Literal


@dataclass
class SNPicConfig:
    """Configuration for a SNPic run."""

    input_folder: str
    """Path to folder containing SNP list files."""
    master_map: str
    """Path to master_disease_mapping.csv."""
    out_prefix: str
    """Output prefix (path + prefix) for results."""
    snp_gene_map: Optional[str] = None
    """Path to SNP-gene map file (required for mode='gene')."""
    k_min: int = 5
    """Minimum number of topics."""
    k_max: int = 20
    """Maximum number of topics."""
    n_bootstrap: int = 50
    """Number of bootstrap iterations."""
    cores: Optional[int] = None
    """Number of CPU cores. Auto-detect if None."""
    seed: int = 123
    """Random seed."""
    keep_all_traits: bool = False
    """Skip stability filter."""
    k_only: Optional[int] = None
    """Fast mode: single K value, skip bootstrap."""
    mode: Literal["gene", "ss"] = "gene"
    """Matrix mode: 'gene' (Gene-as-Word) or 'ss' (Sumstat-as-Word)."""
    model: Literal["lda", "gaussian", "both"] = "lda"
    """Probabilistic model."""

    def to_cmd_args(self) -> list[str]:
        args = [
            "--input_folder", self.input_folder,
            "--master_map", self.master_map,
            "--out_prefix", self.out_prefix,
            "--k_min", str(self.k_min),
            "--k_max", str(self.k_max),
            "--n_bootstrap", str(self.n_bootstrap),
            "--seed", str(self.seed),
            "--mode", self.mode,
            "--model", self.model,
        ]
        if self.snp_gene_map:
            args += ["--snp_gene_map", self.snp_gene_map]
        if self.cores:
            args += ["--cores", str(self.cores)]
        if self.keep_all_traits:
            args.append("--keep_all_traits")
        if self.k_only is not None:
            args += ["--k_only", str(self.k_only)]
        return args


def find_snpic_rscript() -> Path:
    """Locate the SNPic R runner script."""
    # Look relative to this file first
    here = Path(__file__).resolve().parent
    candidates = [
        here.parent / "code" / "core" / "run_snpic.R",
        Path.cwd() / "code" / "core" / "run_snpic.R",
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(
        "run_snpic.R not found. Ensure SNPic is installed correctly."
    )


def run_snpic(
    config: SNPicConfig,
    rscript: str = "Rscript",
    capture_output: bool = False,
    check: bool = True,
) -> subprocess.CompletedProcess:
    """
    Execute a SNPic analysis.

    Args:
        config: SNPicConfig with run parameters.
        rscript: Path/name of Rscript executable.
        capture_output: If True, capture stdout/stderr.
        check: If True, raise CalledProcessError on non-zero exit.

    Returns:
        subprocess.CompletedProcess with result.

    Example:
        >>> config = SNPicConfig(
        ...     input_folder="./data/sig_snp_list",
        ...     master_map="./data/master_map/master_disease_mapping.csv",
        ...     out_prefix="./results/my_analysis",
        ...     snp_gene_map="./data/snp_gene_map.txt",
        ...     mode="gene",
        ...     model="lda",
        ...     k_only=12,
        ... )
        >>> run_snpic(config)
    """
    runner = find_snpic_rscript()
    cmd = [rscript, str(runner)] + config.to_cmd_args()

    if capture_output:
        return subprocess.run(cmd, capture_output=True, text=True, check=check)
    else:
        return subprocess.run(cmd, check=check)


# ─── CLI entry point ──────────────────────────────────────────────────────────

def cli():
    """Command-line interface for SNPic."""
    import argparse

    parser = argparse.ArgumentParser(
        description="SNPic: SNP Topic Modeling for Interpretable Clustering"
    )
    parser.add_argument("--input_folder", required=True)
    parser.add_argument("--master_map", required=True)
    parser.add_argument("--out_prefix", required=True)
    parser.add_argument("--snp_gene_map")
    parser.add_argument("--k_min", type=int, default=5)
    parser.add_argument("--k_max", type=int, default=20)
    parser.add_argument("--n_bootstrap", type=int, default=50)
    parser.add_argument("--cores", type=int)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--keep_all_traits", action="store_true")
    parser.add_argument("--k_only", type=int)
    parser.add_argument("--mode", choices=["gene", "ss"], default="gene")
    parser.add_argument("--model", choices=["lda", "gaussian", "both"], default="lda")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--verbose", "-v", action="store_true")

    args = parser.parse_args()

    config = SNPicConfig(
        input_folder=args.input_folder,
        master_map=args.master_map,
        out_prefix=args.out_prefix,
        snp_gene_map=args.snp_gene_map,
        k_min=args.k_min,
        k_max=args.k_max,
        n_bootstrap=args.n_bootstrap,
        cores=args.cores,
        seed=args.seed,
        keep_all_traits=args.keep_all_traits,
        k_only=args.k_only,
        mode=args.mode,
        model=args.model,
    )

    if args.verbose:
        print(f"🧬 SNPic v{__import__('snpic').__version__}")
        print(f"   Runner: {find_snpic_rscript()}")
        print(f"   Config: {config}")

    result = run_snpic(config, rscript=args.rscript)
    if args.verbose:
        print(f"✅ SNPic completed (exit code: {result.returncode})")


if __name__ == "__main__":
    cli()
