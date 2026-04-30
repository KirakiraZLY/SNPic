"""
SNPic: SNP Topic Modeling for Interpretable Clustering of Complex Phenotypes.

A Python wrapper for the SNPic R pipeline. Provides programmatic access
to run SNPic analyses and parse results.
"""

__version__ = "0.1.0"
__author__ = "Leyi Zhang"

from .cli import run_snpic, SNPicConfig
