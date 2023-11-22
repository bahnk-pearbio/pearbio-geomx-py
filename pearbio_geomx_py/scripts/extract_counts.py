"""
Takes GeoMx `XLSX` result files and converts results to `AnnData` and
`SummarizedExperiment` objects.
"""

# coding: utf-8

import click

from pearbio_geomx_py.results import GeoMxResultFiles

# pylint: disable=no-value-for-parameter
@click.command()
@click.option(
    "--exp-factors",
    default="Type,Infiltration",
    help="experimental factors to extract separated by commas",
)
@click.argument("rna_excel_file")
@click.argument("prot_excel_file")
@click.argument("base_path")
def main(exp_factors, rna_excel_file, prot_excel_file, base_path):
    """
    Opens `XLSX` GeoMx results files and extracts segment info and counts.
    Then, saves everything as `AnnData` and `SummarizedExperiment` objects.
    """
    geomx = GeoMxResultFiles(
        rna_path = rna_excel_file,
        prot_path = prot_excel_file,
        experimental_factors = exp_factors.split(",")
    )

    geomx.create_rna_anndata()
    geomx.create_prot_anndata()

    geomx.save(base_path)

if __name__ == "__main__":
    main()
