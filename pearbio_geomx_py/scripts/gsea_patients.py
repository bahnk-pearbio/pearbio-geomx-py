"""
Runs GSEA on patients comparisons.
"""

# coding: utf-8

import click

from pearbio_geomx_py.gsea import run_gsea_for_patients


# pylint: disable=no-value-for-parameter
# pylint: disable=too-many-arguments
@click.command()
@click.option(
    "-p",
    "--permutations",
    "permutation_num",
    default=10000,
    type=int,
    help="Number of permutations",
)
@click.option(
    "-c",
    "--correction",
    "correction",
    default="bonferroni",
    type=str,
    help="Correction method",
)
@click.option(
    "-t",
    "--threads",
    "threads",
    default=8,
    type=int,
    help="Number of threads",
)
@click.option("-g", "--gmt", "gmts", help="GMT file paths", multiple=True)
@click.option("-k", "--name", "names", help="GMT names", multiple=True)
@click.argument("rds_path")
@click.argument("csv_path")
def main(permutation_num, correction, threads, gmts, names, rds_path, csv_path):
    """
    Takes a `GeoMxExperiment` object with patient comparisons
    log2 fold changes as a `RDS` file and a set of signatures
    as `GMT` files. Then, runs GSEA on every comparisons, performs
    multiple hypothesis correction and saves results as `CSV` file.
    """
    # add names to signatures
    signatures = dict(zip(names, gmts))

    # run GSEA
    args = {
        "path": rds_path,
        "signatures": signatures,
        "permutation_num": permutation_num,
        "correction_method": correction,
        "threads": threads,
    }
    dframe = run_gsea_for_patients(**args)

    # save
    dframe.to_csv(csv_path, index=False)


if __name__ == "__main__":
    main()
