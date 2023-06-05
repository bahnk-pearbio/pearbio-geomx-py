"""
Takes GeoMx `XLSX` results file and converts results to `AnnData` and
`SummarizedExperiment` objects.
"""

# coding: utf-8

import click

from pearbio_geomx_py.segments import GeoMxResultsFile


# pylint: disable=no-value-for-parameter
@click.command()
@click.option(
    "--exp-factors",
    default="Type,Infiltration",
    help="experimental factors to extract separated by commas",
)
@click.argument("excel_file")
@click.argument("base_path")
def main(exp_factors, excel_file, base_path):
    """
    Opens `XLSX` GeoMx results file and extracts segment info and counts.
    Then, saves everything as `AnnData` and `SummarizedExperiment` objects.
    """
    factors = exp_factors.split(",")
    results = GeoMxResultsFile(excel_file, factors)
    adata = results.create_anndata(excel_file, factors)
    adata.write_h5ad(f"{base_path}.h5ad")
    results.save_anndata_as_rds(adata, f"{base_path}.rds")


if __name__ == "__main__":
    main()
