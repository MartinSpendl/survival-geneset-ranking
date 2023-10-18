import numpy as np
import pandas as pd
from single_sample_gsea import ss_gsea as ssGSEA


def get_GSEA_genesets(expressed_genes: list, filename: str, dir="."):
    genesets = {}

    with open(f"{dir}/data/genesets/{filename}") as f:
        lines = f.readlines()

        for line in lines:
            values = line.split("\t")
            genesets[values[0]] = values[2:]

    return {
        name: set([g for g in values if g in expressed_genes])
        for name, values in genesets.items()
    }


def calculate_ssGSEA(
    tcga: str,
    geneset_name: str,
    geneset_filename: str,
    min_n_genes: int = 10,
    path_to_input_data: str = "tcga-data",
    path_to_output_data: str = "data",
    number_of_cores: int = 1,
):
    expressions = pd.read_csv(
        f"{path_to_input_data}/TCGA-{tcga}.csv", sep=",", index_col=[0]
    )

    genesets = get_GSEA_genesets(expressions.columns, geneset_filename)
    genesets = {k: v for k, v in genesets.items() if len(v) >= min_n_genes}

    expressions = expressions.T
    expressions.index.name = "gene"

    ssGSEA_scores = ssGSEA(expressions, genesets, num_cores=number_of_cores)

    ssGSEA_scores.reset_index().to_csv(
        f"{path_to_output_data}/{geneset_name}/TCGA-{tcga}-ssGSEA.tsv",
        sep="\t",
        index=False,
    )
