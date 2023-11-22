"""
For GeoMx segments sample sheet validity.

GeoMx calls segments the ROIs the scientist drew for generating RNA-seq
libraries.
GeoMx process FASTQ files and provides counts counts data and information
about segments as `XLSX` file.
GeoMx also gives protein expression data for each segment.
This file contains methods to parse segments metadata, extract counts and
and saves all these as usable python and R objects for downstream analysis.
"""

# coding: utf-8

import re
import pandas as pd
import anndata as ad

from typing import List, Dict
from pathlib import Path

from rpy2 import robjects
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.packages import importr

from pearbio_geomx_py.segments import get_segments_data_frame

class GeoMxResultFiles:
    """GeoMx `XLSX` results file."""

    rna_path: Path = None
    prot_path: Path = None

    segments: pd.DataFrame = None
    experimental_factors: List[str] = None

    rna: ad.AnnData =  None
    prot: ad.AnnData = None

    prot_sheets: Dict[str, str] = {
        'qc':           'QC'            ,
        'area':         'Area'          ,
        'housekeeping': 'Housekeeping'  ,
        'neg_norm':     'Neg Norm'      ,
        'nuclei':       'Nuclei'        ,
    }

    rna_slots: Dict[str, str] = {
        'raw':  'counts'     ,
        'norm': 'norm_counts'    ,
    }

    def __init__(
        self,
        rna_path: str,
        prot_path: str,
        experimental_factors: List[str]
        ) -> None:
        """
        Constructor taking `XLSX` results file path.

        :param rna_path: Path of the RNA sample sheet `XLSX` file.
        :type path: str
        :param prot_path: Path of the protein panel sample sheet `XLSX` file.
        :type path: str
        :param experimental_factors: Experimental factors (columns) to extract.
        :type experimental_factors: List[str]
        :raise FileNoFoundError: If sample sheet files don't exist.
        :rtype: None
        """
        rna_path = Path(rna_path)
        prot_path = Path(prot_path)

        if not rna_path.exists():
            raise FileNotFoundError(f'GeoMx results file {path} doesn\'t exist.')

        if not prot_path.exists():
            raise FileNotFoundError(f'GeoMx results file {path} doesn\'t exist.')

        self.rna_path = rna_path
        self.prot_path = prot_path
        self.experimental_factors = experimental_factors

        self.segments = get_segments_data_frame(
            rna_path = self.rna_path,
            prot_path = self.prot_path,
            experimental_factors = self.experimental_factors
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
            dframe = pd.read_excel(io = path, sheet_name = sheet_name)
        except Exception as exc:
            raise IOError(f'GeoMx results file {path} can\'t be read.') from exc

        return dframe

    def create_rna_anndata(self) -> None:
        """Update `anndata` object containing RNA expression data."""
        self.rna = self._create_rna_anndata(
            segments = self.segments,
            slot_names = self.rna_slots,
            path = self.rna_path
        )

    def _create_rna_anndata(
        self,
        segments: pd.DataFrame,
        slot_names: List[str],
        path: Path
        ) -> ad.AnnData:
        """
        Create `anndata` object containing RNA expression data.

        :param path: Path of GeoMx `XLSX` results file.
        :type path: Path
        :rtype: ad.AnnData
        """
        # genes
        genes = self\
            ._get_excel_sheet(path = path, sheet_name = 'TargetProperties')\
            .rename(columns = {'TargetName': 'Gene'})\
            .loc[:, ['Gene', 'GeneID']]\
            .sort_values(by = 'Gene')\
            .reset_index(drop = True)

        # raw_counts
        raw_counts = self\
            ._get_excel_sheet(path = path, sheet_name = 'TargetCountMatrix')\
            .rename(columns = {'TargetName': 'Gene'})\
            .set_index('Gene')\
            .loc[genes.Gene, segments.sample_id]

        # normalized counts
        norm_counts = self\
            ._get_excel_sheet(path = path, sheet_name = 'Q3 TargetCountMatrix')\
            .set_index('Gene')\
            .loc[genes.Gene, segments.sample_id]
        
        args = {
            'X': raw_counts.T,
            'obs': segments.set_index('sample_id', drop = False),
            'var': genes.set_index('Gene', drop = False),
            'layers': {
                slot_names['raw']: raw_counts.T,
                slot_names['norm']: norm_counts.T,
            },
        }

        return ad.AnnData(**args)

    def create_prot_anndata(self) -> None:
        """Update `anndata` object containing RNA expression data."""
        self.prot = self._create_prot_anndata(
            segments = self.segments,
            sheet_names = self.prot_sheets,
            path = self.prot_path
        )

    def _create_prot_anndata(
        self,
        segments: pd.DataFrame,
        sheet_names: Dict[str, str],
        path: Path
        ) -> ad.AnnData:
        """
        Create `anndata` object containing protein expression data.

        :param path: Path of GeoMx `XLSX` results file.
        :type path: Path
        :rtype: ad.AnnData
        """
        prot = self\
            ._get_excel_sheet(path = path, sheet_name = 'Probe Groups')\
            .rename(columns = {'Target name (display name)': 'Protein'})\
            .rename(columns = lambda x: re.sub('#', '', x))\
            .sort_values(by = 'Protein')\
            .set_index('Protein')\
            .reset_index()

        sheets = {}

        for sheet_slot, sheet_name, in self.prot_sheets.items():

            dframe = self\
                ._get_excel_sheet(path = path, sheet_name = sheet_name)\
                .set_index('ROI_ID')\
                .filter(items = prot.Protein.to_list())\
                .loc[self.segments.roi_label,prot.Protein]

            # indexes
            dframe.columns.name = 'Protein'
            dframe.index = segments.sample_id
            dframe.index.name = 'sample_id'

            sheets[sheet_slot] = dframe

        args = {
            'X': list(sheets.values())[0],
            'obs': segments.set_index('sample_id', drop = False),
            'var': prot.set_index('Protein', drop = False),
            'layers': sheets,
        }

        return ad.AnnData(**args)

    def save(self, basepath: str) -> None:
        """Save results as HDF5 Python and `R` RDS files."""
        # R RDS files
        self._save_as_rds(self.rna, f'{basepath}.rna.rds')
        self._save_as_rds(self.prot, f'{basepath}.prot.rds')

        # HDF5 python files
        self.rna.write_h5ad(f'{basepath}.rna.h5ad')
        self.prot.write_h5ad(f'{basepath}.prot.h5ad')

    def _save_as_rds(self, adata: ad.AnnData, path: str) -> None:
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
        r_base = importr('base')
        r_sumexp = importr('SummarizedExperiment')

        r_assays = {}
        
        for layer, values in adata.layers.items():
            r_assays[layer] = numpy2ri.py2rpy(values.T)

        r_assays = robjects.vectors.ListVector(r_assays)

        # convert column metadata to R data frame
        col_metadata = pandas2ri.py2rpy(adata.obs)

        # convert row metadata to R data frame
        row_metadata = pandas2ri.py2rpy(adata.var)

        # create a SummarizedExperiment object in R
        r_sumexp_obj = r_sumexp.SummarizedExperiment(
            assays = r_assays,
            colData = col_metadata,
            rowData = row_metadata
        )

        # save SummarizedExperiment as RDS file
        r_base.saveRDS(r_sumexp_obj, path)

if __name__ == '__main__':

    geomx = GeoMxResultFiles(
        rna_path = '/home/nourdine/tmp/GX0000_rna.xlsx',
        prot_path = '/home/nourdine/tmp/GX0000_prot.xlsx',
        experimental_factors = ['Type', 'Infiltration', 'Patient']
    )

    geomx.create_rna_anndata()
    geomx.create_prot_anndata()
    geomx.save('/tmp/test')
