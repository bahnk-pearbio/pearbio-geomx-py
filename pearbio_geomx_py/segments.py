"""
For GeoMx segments sample sheet validity.

GeoMx calls segments the ROIs the scientist drew for generating RNA-seq
libraries.
GeoMx process FASTQ files and provides counts counts data and information
about segments as `XLSX` file.
This file contains methods to parse segments metadata, extract counts and
and saves all these as usable python and R objects for downstream analysis.
"""

# coding: utf-8

from typing import List, Dict
from pathlib import Path

import pandas as pd
import anndata as ad

from rpy2 import robjects
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

# pylint: disable=no-name-in-module
from pydantic import BaseModel


# pylint: disable=too-few-public-methods
class Segment(BaseModel):
    """Data class for GeoMx sample sheet row."""

    Sample_ID: str
    AOISurfaceArea: float
    AOINucleiCount: int
    RawReads: int
    AlignedReads: int
    experimental_factors: dict

    def dict(self) -> Dict:
        """Returns `dict` containing only required info."""
        segment = super().dict()
        _ = segment.pop("experimental_factors")

        for factor, level in self.experimental_factors.items():
            segment[factor] = level

        return segment


class GeoMxResultsFile:
    """GeoMx `XLSX` results file."""

    path: Path = None
    segments: pd.DataFrame = None
    experimental_factors: List[str] = None

    def __init__(self, path: str, experimental_factors: List[str]) -> None:
        """
        Constructor taking `XLSX` results file path.

        :param path: Path of the sample sheet `XLSX` file.
        :type path: str
        :param experimental_factors: Experimental factors (columns) to extract.
        :type experimental_factors: List[str]
        :raise FileNoFoundError: If sample sheet file doesn't exist.
        :rtype: None
        """
        path = Path(path)

        if not path.exists():
            raise FileNotFoundError(f"GeoMx results file {path} doesn't exist.")

        self.path = path
        self.experimental_factors = experimental_factors
        self.segments = self.get_segments_data_frame(
            self.path, self.experimental_factors
        )

    def _get_excel_sheet(self, path: Path, sheet_name: str) -> pd.DataFrame:
        """
        Extract spreadsheet from GeoMx `XLSX` results file.

        Method assumes the GeoMx `XLSX` results file exists.

        :param path: Path of GeoMx `XLSX` results file.
        :type path: Path
        :param sheet_name: Name of spreadsheet to extract.
        :type sheet_name: str
        :raise IOError: If `XLSX` results file can't be read.
        :rtype: pandas.DataFrame
        """
        # pylint: disable=fixme
        # TODO: figure out why this raises an OSError instead
        try:
            dframe = pd.read_excel(path, sheet_name=sheet_name)
        except Exception as exc:
            raise IOError(f"GeoMx results file {path} can't be read.") from exc

        return dframe

    def get_segments_data_frame(
        self, path: Path, experimental_factors: List[str]
    ) -> pd.DataFrame:
        """
        Extract segments information from `XLSX` results file.

        Method creates data frame with required experimental factors (columns)
        for downstream processing.

        :param path: Path of GeoMx `XLSX` results file.
        :type path: Path
        :param experimental_factors: Experimental factors (columns) to extract.
        :type experimental_factors: List[str]
        :raise KeyError: If an experimental factor can't be found.
        :rtype: pd.DataFrame
        """
        # sheet name should stay the same
        dframe = self._get_excel_sheet(path=path, sheet_name="SegmentProperties")

        # list of dict with info for each segment
        segments = []

        for _, row in dframe.iterrows():
            factors = {}

            # experimental factors are columns of the spreadsheet,
            # user specifies the column she/he wants to extract
            for factor in experimental_factors:
                try:
                    level = row[factor]
                except Exception as exc:
                    raise KeyError(
                        f"Experimental factor {factor} is not a column"
                    ) from exc

                factors[factor] = level

            segment = Segment(**row.to_dict(), experimental_factors=factors)
            segments.append(segment.dict())

        return pd.DataFrame.from_records(segments)

    def create_anndata(self, path: Path, experimental_factors: List[str]) -> ad.AnnData:
        """
        Create `anndata` object containing experiment.

        :param path: Path of GeoMx `XLSX` results file.
        :type path: Path
        :param experimental_factors: Experimental factors (columns) to extract.
        :type experimental_factors: List[str]
        :rtype: ad.AnnData
        """
        # segments
        self.segments = self.get_segments_data_frame(
            path=path, experimental_factors=experimental_factors
        )

        # genes
        genes = (
            self._get_excel_sheet(path=path, sheet_name="TargetProperties")
            .rename(columns={"TargetName": "Gene"})
            .loc[:, ["Gene", "GeneID"]]
            .sort_values(by="Gene")
            .reset_index(drop=True)
        )

        # raw_counts
        raw_counts = (
            self._get_excel_sheet(path=path, sheet_name="TargetCountMatrix")
            .rename(columns={"TargetName": "Gene"})
            .set_index("Gene")
            .loc[genes.Gene, self.segments.Sample_ID]
        )

        # normalized counts
        norm_counts = (
            self._get_excel_sheet(path=path, sheet_name="Q3 TargetCountMatrix")
            .set_index("Gene")
            .loc[genes.Gene, self.segments.Sample_ID]
        )

        args = {
            "X": raw_counts.T,
            "obs": self.segments.set_index("Sample_ID", drop=False),
            "var": genes.set_index("Gene", drop=False),
            "layers": {"norm_X": norm_counts.T},
        }

        return ad.AnnData(**args)

    def save_anndata_as_rds(self, adata: ad.AnnData, path: str) -> None:
        """
        Save `AnnData` object as `R` `SummarizedExperiment object in `RDS`
        file.

        :param adata: `AnnData` object to save.
        :type adata: anndata.AnnData
        :param path: Path of output `RDS` file.
        :type path: str
        :rtype: None
        """
        # R packages
        pandas2ri.activate()
        r_base = importr("base")
        r_sumexp = importr("SummarizedExperiment")

        # counts and normalized counts
        r_counts = numpy2ri.py2rpy(adata.X.T)
        r_norm_counts = numpy2ri.py2rpy(adata.layers["norm_X"].T)
        r_assays = robjects.vectors.ListVector(
            {"counts": r_counts, "norm_counts": r_norm_counts}
        )

        # convert column metadata to R data frame
        col_metadata = pandas2ri.py2rpy(adata.obs)

        # convert row metadata to R data frame
        row_metadata = pandas2ri.py2rpy(adata.var)

        # create a SummarizedExperiment object in R
        r_sumexp_obj = r_sumexp.SummarizedExperiment(
            assays=r_assays, colData=col_metadata, rowData=row_metadata
        )

        # save SummarizedExperiment as RDS file
        r_base.saveRDS(r_sumexp_obj, path)
