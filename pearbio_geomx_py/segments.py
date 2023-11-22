"""
GeoMx segments info validation class.
"""

# coding: utf-8

import pandas as pd

from typing import Dict, List
from pathlib import Path

# pylint: disable=no-name-in-module
from pydantic import BaseModel

# pylint: disable=too-few-public-methods
class Segment(BaseModel):
    """
    Data class for GeoMx segement.

    It store information for both RNA and protein.
    """

    sample_id: str
    roi_label: int

    rna_aoi_surface_area: float
    rna_aoi_nuclei_count: int

    rna_raw_reads: int
    rna_aligned_reads: int
    rna_deduplicated_reads: int
    rna_trimmed_reads: int
    rna_stitched_reads: int

    rna_roi_coordinate_x: int
    rna_roi_coordinate_y: int

    prot_roi_coordinate_x: int
    prot_roi_coordinate_y: int

    prot_aoi_surface_area: int
    prot_aoi_nuclei_count: int

    experimental_factors: Dict[str, str]

    def __init__(self, **data):
        """
        Customize data import.

        This is the most direct way to control attribute names.
        """
        return super().__init__(**{

            'sample_id':                data['Sample_ID']               ,
            'roi_label':                data['ROILabel']                ,

            'rna_aoi_surface_area':     data['AOISurfaceArea']          ,
            'rna_aoi_nuclei_count':     data['AOINucleiCount']          ,

            'rna_raw_reads':            data['RawReads']                ,
            'rna_aligned_reads':        data['AlignedReads']            ,
            'rna_deduplicated_reads':   data['DeduplicatedReads']       ,
            'rna_trimmed_reads':        data['TrimmedReads']            ,
            'rna_stitched_reads':       data['StitchedReads']           ,

            'rna_roi_coordinate_x':     data['ROICoordinateX']          ,
            'rna_roi_coordinate_y':     data['ROICoordinateY']          ,

            'prot_roi_coordinate_x':    data['ROI X Coordinate']        ,
            'prot_roi_coordinate_y':    data['ROI Y Coordinate']        ,

            'prot_aoi_surface_area':    data['AOI surface area']        ,
            'prot_aoi_nuclei_count':    data['AOI nuclei count']        ,

            'experimental_factors':     data['experimental_factors']    ,
        })

    def dict(self) -> Dict:
        """Returns `dict` containing only required info."""
        segment = super().dict()
        _ = segment.pop('experimental_factors')

        for factor, level in self.experimental_factors.items():
            segment[factor] = level

        return segment

def get_segments_data_frame(
    rna_path: Path,
    prot_path: Path,
    experimental_factors: List[str]
) -> pd.DataFrame:
    """
    Extract segments information from `XLSX` results file.

    Method creates data frame with required experimental factors (columns)
    for downstream processing.

    :param rna_path: Path of the RNA sample sheet `XLSX` file.
    :type path: str
    :param prot_path: Path of the protein panel sample sheet `XLSX` file.
    :type path: str
    :param experimental_factors: Experimental factors (columns) to extract.
    :type experimental_factors: List[str]
    :raise KeyError: If an experimental factor can't be found.
    :raise ValueError: If missing protein info for some segments.
    :rtype: pd.DataFrame
    """
    rna = pd.read_excel(io = rna_path, sheet_name = 'SegmentProperties')
    prot = pd.read_excel(io = prot_path, sheet_name = 'QC')

    if rna.shape[0] > prot.shape[0]:
        raise ValueError('Missing protein info')

    data = pd.merge(
        left = rna,
        right = prot,
        how = 'inner',
        left_on = 'ROILabel',
        right_on = 'ROI_ID'
    )

    segments = []

    for _, row in data.iterrows():

        factors = {}

        # experimental factors are columns of the spreadsheet,
        # user specifies the column she/he wants to extract
        for factor in experimental_factors:

            try:
                level = row[factor]
            except Exception as exc:
                raise KeyError(f'{factor} is not a column') from exc

            factors[factor] = level

        segment = Segment(**row.to_dict(), experimental_factors = factors)
        segments.append(segment.dict())

    return pd.DataFrame.from_records(segments)

if __name__ == '__main__':

    dframe = get_segments_data_frame(
        rna_path = '/home/nourdine/tmp/GX0000_rna.xlsx',
        prot_path = '/home/nourdine/tmp/GX0000_prot.xlsx',
        experimental_factors = ['Type', 'Infiltration', 'Patient']
    )

    print(dframe)
