"""
Module to handle Gene Set Enrichment Analysis.
"""

# coding: utf-8

import gseapy as gp

from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

import pandas as pd

from statsmodels.sandbox.stats.multicomp import multipletests


# pylint: disable=too-many-arguments
# pylint: disable=too-many-locals
def run_gsea_for_patients(
    path: str,
    signatures: dict,
    slot_name: str = "patients_comparisons",
    lfc_slot_name: str = "lfc",
    lfc_col_name: str = "lfc",
    patient1_slot_name: str = "patient1",
    patient2_slot_name: str = "patient2",
    correction_method: str = "bonferroni",
    threads: int = 8,
    permutation_num: int = 10000,
) -> pd.DataFrame:
    """
    Run GSEA for patient comparisons and for multiple gene sets.

    :param path: RDS path containing `GeoMxExperiment` object.
    :type path: str
    :param signatures: Dictionary containing paths and names of gene sets.
    :type signatures: dict
    :param slot_name: Patients comparisons slot name.
    :type slot_name: str
    :param lfc_slot_name: Log2 fold change slot name containing data frame.
    :type lfc_slot_name: str
    :param lfc_col_name: Log2 fold change column name.
    :type lfc_col_name: str
    :param patient1_col_name: Patient 1 name slot name.
    :type patient1_col_name: str
    :param patient2_col_name: Patient 2 name slot name.
    :type patient2_col_name: str
    :param correction_method: Multiple hypothesis correction method.
    :type correction_method: str
    :param threads: Number of threads.
    :type threads: int
    :param permutation_num: Number of permutations.
    :type permutation_num: int
    """
    # R requirements
    r_base = importr("base")
    experiment = r_base.readRDS(path)
    pandas2ri.activate()

    # GSEA results
    dframes = []

    for label, comparisons in experiment.slots[slot_name].items():
        for _, comparison in comparisons.items():
            # patients
            patient1 = comparison.rx2(patient1_slot_name)[0]
            patient2 = comparison.rx2(patient2_slot_name)[0]

            # gene rank
            rnk = pandas2ri.rpy2py_dataframe(comparison.rx2(lfc_slot_name))
            rnk = rnk.sort_values(lfc_col_name, ascending=False)
            rnk = rnk.set_index("gene")

            for signature, signature_path in signatures.items():
                gene_sets = gp.read_gmt(path=signature_path)

                args = {
                    "threads": threads,
                    "min_size": 5,
                    "max_size": 1000,
                    "permutation_num": permutation_num,
                    "seed": 123,
                    "verbose": True,
                    "rnk": rnk,
                    "gene_sets": gene_sets,
                }

                pre_res = gp.prerank(**args)

                dframe = pre_res.res2d
                dframe["Label"] = label
                dframe["Patient1"] = patient1
                dframe["Patient2"] = patient2
                dframe["signature"] = signature

                dframes.append(dframe)

    gsea = pd.concat(dframes)

    # thresholding p-values
    # pylint: disable=unnecessary-lambda-assignment
    pval = lambda x: 1 / permutation_num if x == 0 else x
    gsea["PValue"] = gsea["NOM p-val"].apply(pval)

    # multiple hypothesis testing correction
    args = {"pvals": gsea["PValue"], "method": correction_method}
    gsea["AdjustedPValue"] = multipletests(**args)[1]

    return gsea
