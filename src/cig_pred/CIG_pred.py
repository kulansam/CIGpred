#!/usr/bin/env python

import sys
import argparse
import pandas as pd
import os
import numpy as np
import anndata
import scanpy as sc
from scipy import io
import warnings
warnings.filterwarnings("ignore")
from cig_pred_human_functions import *
from cig_pred_mouse_functions import *



#INPUT: python CIG_pred.py -organism hs -assaytype bulk -inputtype tpm -file=./test_exp_yang.txt
path='/Users/dongyuzhao/Library/CloudStorage/GoogleDrive-kulandaisamy.arulsamy@enders.tch.harvard.edu/My Drive/CIG_machine_learning/ML_model_master_TF_pred/our_model_apply/Final_DNA_seq_features_jan2023/final_model_setup/cig_predtest/data/'
os.chdir(path)
def printHelp():
    print ("Welcome!")
    print ("Functions:\n")
    print ("CIG-Pred: Predict cell identity genes or gene of interests based on different signals provided (RNAseq)\n")
    print ("CIG-reg-Pred: Predict master transcription factors based Gene regulatory network and CIG-Pred prediction result\n")
    print ("Kulandaisamy, et al. kulandaisamy@gmail.com, Dr.Kaifu Chen lab, @ Boston Children's Hospital, Harvard Medical School")
    print()

#print(sys.argv)
#print(len(sys.argv))

if (len(sys.argv) ==9):
    parser = argparse.ArgumentParser(usage="...",description='', epilog="Chen lab, Houston Methodist")
    parser.add_argument('-organism', dest='organism', type=str, help='Mention the organism name of a sample/data')
    parser.add_argument('-assaytype', dest='assaytype', type=str, help='Type of the data either bulk/pseudobulk or single-cell RNA expression data')
    parser.add_argument('-inputtype', dest='inputtype', type=str, help='Mention whether the expression values are in TPM or in raw count')
    parser.add_argument('-file', dest='file', type=str, help='Mention a filename that contains the expression values of genes')


    args = parser.parse_args()
    if args is not None:
        organism_name_pred = args.organism
        assaytype_pred = args.assaytype
        inputtype_pred = args.inputtype
        exp_filename_pred = args.file
        
    ##add sequence features and model training table
    if organism_name_pred =='hs':
        all_seq_features_human = pd.read_csv("Requried_seq_features_cigpred_human.txt",sep="\t")
        gene_name_length_human=pd.read_csv("./gene_name_map_with_length_human.txt",sep="\t")
        training_table_human=pd.read_csv("training_table_human.txt",sep="\t")   
        training_table_cigreg_human=pd.read_csv("training_table_master_tf_cigs_human.txt",sep="\t")
        all_db_grn_human=pd.read_csv("./GRN_network_compiled_human.txt",sep="\t")
        tf_human=pd.read_csv("./tf_human_list.txt",sep="\t")
    elif organism_name_pred =='mm':
        all_seq_features_mouse = pd.read_csv("Requried_seq_features_cigpred_mouse.txt",sep="\t")
        gene_name_length_mouse=pd.read_csv("gene_name_map_with_length_mouse.txt",sep="\t")
        #gene_name_length_mouse['Genename'] = gene_name_length_mouse['Genename'].str.upper()
        training_table_mouse=pd.read_csv("training_table_mouse.txt",sep="\t")
        training_table_cigreg_mouse=pd.read_csv("training_table_master_tf_cigs_mouse.txt",sep="\t")
        all_db_grn_mouse=pd.read_csv("all_db_GRN_combined_mouse.txt",sep="\t")
        tf_mouse=pd.read_csv("./tf_mouse_list_geneid.txt",sep="\t")

    
    if organism_name_pred == 'hs' and assaytype_pred == 'bulk' and exp_filename_pred is not None:
        print ("HUMAN",organism_name_pred,assaytype_pred,inputtype_pred,exp_filename_pred)
        exp_matrix=pd.read_csv(exp_filename_pred,sep="\t")
        celltype_names=exp_matrix.columns.to_list()
        print(exp_matrix.shape)
        cig_pred_output_table=pesudobulk_cigpred_human(inputtype_pred,exp_matrix,all_seq_features_human,gene_name_length_human,training_table_human)
        cig_pred_output_table.to_csv(str(exp_filename_pred)+'cig_pred_result.out', index = False,sep="\t")
        
        cig_reg_pred_output_table=pesudobulk_cig_reg_pred_human(cig_pred_output_table,all_db_grn_human,tf_human,training_table_cigreg_human,celltype_names)
        print(cig_reg_pred_output_table.shape)
        cig_reg_pred_output_table.to_csv(str(exp_filename_pred)+'cig_REG_pred_result.out', index = False,sep="\t")

    elif organism_name_pred == 'hs' and assaytype_pred == 'single' and inputtype_pred == 'umicount' and exp_filename_pred is not None:
        print ("HUMAN",organism_name_pred,assaytype_pred,inputtype_pred,exp_filename_pred)
        barcodes = pd.read_csv(str(exp_filename_pred)+"barcodes.tsv",sep="\t",header=None)
        barcodes.columns=['barcode']
        features = pd.read_csv(str(exp_filename_pred)+"features.tsv",sep="\t",header=None)
        features.columns=['Geneid','gene_symbols','note']
        counts_mat = io.mmread(str(exp_filename_pred)+"matrix.mtx")
        counts_mat= counts_mat.toarray()
        counts_mat= np.matrix(counts_mat.transpose())
        # create anndata object
        adata = anndata.AnnData(counts_mat,obs=barcodes['barcode'].tolist(),var=features)
        adata.obs.index = barcodes['barcode'].tolist()
        adata.var.index = features['Geneid'].tolist()
        adata.var_names_make_unique()##REMOVE DUPLICATES
        adata_single_cell_cig=cig_pred_singlecell_human(adata,features,all_seq_features_human,training_table_human,exp_filename_pred)
        adata_single_cell_cig.T.write_h5ad(str(exp_filename_pred)+"_cig_matrix_out.h5ad")

    elif organism_name_pred == 'mm' and assaytype_pred == 'bulk' and exp_filename_pred is not None:
        print ("MOUSE",organism_name_pred,assaytype_pred,inputtype_pred,exp_filename_pred)
        exp_matrix=pd.read_csv(exp_filename_pred,sep="\t")
        celltype_names=exp_matrix.columns.to_list()
        print(exp_matrix.shape)
        cig_pred_output_table=pesudobulk_cigpred_mouse(inputtype_pred,exp_matrix,all_seq_features_mouse,gene_name_length_mouse,training_table_mouse)
        cig_pred_output_table.to_csv(str(exp_filename_pred)+'cig_pred_result.out', index = False,sep="\t")
        
        cig_reg_pred_output_table=pesudobulk_cig_reg_pred_mouse(cig_pred_output_table,all_db_grn_mouse,tf_mouse,training_table_cigreg_mouse,celltype_names)
        cig_reg_pred_output_table.to_csv(str(exp_filename_pred)+'cig_REG_pred_result.out', index = False,sep="\t")

    elif organism_name_pred == 'mm' and assaytype_pred == 'single' and inputtype_pred == 'umicount' and exp_filename_pred is not None:
        print ("MOUSE",organism_name_pred,assaytype_pred,inputtype_pred,exp_filename_pred)
        barcodes = pd.read_csv(str(exp_filename_pred)+"barcodes.tsv",sep="\t",header=None)
        barcodes.columns=['barcode']
        features = pd.read_csv(str(exp_filename_pred)+"features.tsv",sep="\t",header=None)
        features.columns=['Geneid','gene_symbols','note']
        counts_mat = io.mmread(str(exp_filename_pred)+"matrix.mtx")
        counts_mat= counts_mat.toarray()
        counts_mat= np.matrix(counts_mat.transpose())
        # create anndata object
        adata = anndata.AnnData(counts_mat,obs=barcodes['barcode'].tolist(),var=features)
        adata.obs.index = barcodes['barcode'].tolist()
        adata.var.index = features['Geneid'].tolist()
        adata.var_names_make_unique()##REMOVE DUPLICATES
        adata_single_cell_cig=cig_pred_singlecell_mouse(adata,features,all_seq_features_mouse,training_table_mouse,exp_filename_pred)
        adata_single_cell_cig.T.write_h5ad(str(exp_filename_pred)+"_cig_matrix_out.h5ad")

else:  
    print ("\nFor help, please try: CIG-Pred -h\n")
    printHelp()