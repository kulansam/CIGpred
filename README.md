# Welcome to CIGpred: A tool for identifying cell identity genes (CIGs) and their master transcription factors at both bulk and single-cell levels, utilizing genetic sequence and RNA expression patterns
In every cell type, a unique set of cell identity genes (CIGs) plays a pivotal role in defining its specific characteristics. Alongside other gene categories like housekeeping genes and heat shock genes, cells express their CIGs, crucial for guiding cellular differentiation and the formation of tissues and organs in organisms. Thus, we have developed two logistic regression-based machine-learning methodologies for the identification of cell identity genes (CIGs) and their master transcription factors:

1. CIG-Pred: This tool accurately discerns cell identity genes by integrating genetic sequence and RNA expression profiles in a given cell type.
2. CIG-reg-Pred: Subsequently, CIG-reg-Pred determines the master transcription factors governing the expression of identified cell identity genes in a given cell type.

# What constitutes the input data for CIG-Pred and CIG-reg-Pred algorithms?

The input for CIG-Pred consists of either raw read counts from bulk-seq or gene per cell unique molecular identifier (UMI) count matrix from single-cell sequencing data. In addition, CIG-Pred automatically utilizes precomputed genetic sequence features during its prediction process.

For CIG-reg-Pred, the required inputs include the predicted cell identity scores obtained from the CIG-Pred algorithm, along with gene regulatory network (GRN) information. It's worth noting that the GRN information is already integrated into the model.
# Installation
- Download CIG-Pred package from GitHub:
  ```sh
  git clone https://github.com/kulansam/CIGpred.git
  ```
  ```sh
  cd CIGpred
  ```
- Install requirements
  ```sh
  pip install -r requirements.txt
  ```
- Install CIGpred
  ```sh
  pip install .
  ```
# Tutorial 
The following command could used to identify the cell identity genes (CIGs) in the human genome using expression determined by bulk-RNA sequencing profiles:
  ```sh
    python CIG_pred.py -organism hs -assaytype bulk -inputtype rawcount -file expression_data.txt (tab separated)
  ```

# Cite us
Please cite us at <a href='#' target='_blank'>bioRxiv</a> if you find CIG-Pred  is useful to your project.</p>

# Term of Usage
By accessing CIG-Pred data, you agree to the following terms:

1. You agree not to share the CIG-Pred data, whether modified or unmodified, with individuals outside your research group. This includes preventing access by unauthorized individuals and refraining from directly providing the data to others.

2. You agree not to develop another website or methods using the CIG-Pred data without prior permission. Contact us for any such intentions.

3. You agree to appropriately cite the CIG-Pred paper and its original contributions if utilized in your work.

4. You certify that you are authorized to accept these terms on behalf of your institution.

# Contact
KulandaiSamy.Arulsamy@childrens.harvard.edu
or
Kaifu.Chen@childrens.harvard.edu
<hr>
Copy Right @ Kaifu Chen Lab @ Boston Childrens Hospital / Harvard Medical School

