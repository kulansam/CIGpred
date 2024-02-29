# Welcome to CIGpred: A tool for identifying cell identity genes (CIGs) and their master transcription factors at both bulk and single-cell levels, utilizing genetic sequence and RNA expression patterns
In every cell type, a unique set of cell identity genes (CIGs) plays a pivotal role in defining its specific characteristics. Alongside other gene categories like housekeeping genes and heat shock genes, cells express their CIGs, crucial for guiding cellular differentiation and the formation of tissues and organs in organisms. Thus, we have developed two logistic regression-based machine-learning methodologies for the identification of cell identity genes (CIGs) and their master transcription factors:

1. CIG-Pred: This tool accurately discerns cell identity genes by integrating genetic sequence and RNA expression profiles in a given cell type.
2. CIG-reg-Pred: Subsequently, CIG-reg-Pred determines the master transcription factors governing the expression of identified cell identity genes in a given cell type.

# What constitutes the input data for CIG-Pred and CIG-reg-Pred algorithms?

The input for CIG-Pred consists of either raw read counts from bulk-seq or unique molecular identifier (UMI) counts from single-cell sequencing data. In addition, CIG-Pred automatically utilizes precomputed genetic sequence features during its prediction process.

For CIG-reg-Pred, the required inputs include the predicted cell identity scores obtained from the CIG-Pred algorithm, along with gene regulatory network (GRN) information. It's worth noting that the GRN information is already integrated into the model.
