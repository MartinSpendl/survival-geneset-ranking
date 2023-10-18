import os
import sys

import numpy as np
import pandas as pd
import scipy.stats as ss
from tqdm import tqdm

sys.path.append(os.path.dirname(os.path.abspath("__file__")))
from src.utils.constants import TCGA_DATASETS
from src.ssGSEA.ssGSEA_related import calculate_ssGSEA


if __name__ == "__main__":
    # from: http://current.geneontology.org/products/pages/downloads.html

    GENESETS = {  # name: filename
        "GO-BiologicalProcess": "c5.go.bp.v2023.1.Hs.symbols.gmt",
        "Hallmarks": "h.all.v2023.1.Hs.symbols.gmt",
    }

    MIN_N_GENES = 10
    OVERRIDE = True
    OUTPUT_DATA_DIR = "data"
    NUMBER_OF_CORES = 8

    for geneset_name, geneset_filename in GENESETS.items():

        # if not, create geneset-named directory
        if geneset_name not in os.listdir(OUTPUT_DATA_DIR):
            os.mkdir(f"{OUTPUT_DATA_DIR}/{geneset_name}")

        for tcga in tqdm(TCGA_DATASETS):

            # if already calculated, skip
            if (
                f"TCGA-{tcga}-ssGSEA.tsv"
                in os.listdir(f"{OUTPUT_DATA_DIR}/{geneset_name}")
                and not OVERRIDE
            ):
                continue

            else:
                calculate_ssGSEA(
                    tcga,
                    geneset_name=geneset_name,
                    geneset_filename=geneset_filename,
                    min_n_genes=MIN_N_GENES,
                    path_to_input_data="tcga-data/data",
                    path_to_output_data=OUTPUT_DATA_DIR,
                    number_of_cores=NUMBER_OF_CORES,
                )
