#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 23:57:22 2023

@author: abhishek
"""



import pandas as pd
from Bio import SeqIO
import pysam
from pyensembl import EnsemblRelease
from genelocator import get_genelocator


def return_gene_number_v38(data):
    
    df_in = data.copy()
    
    for i in range(df_in.shape[0]):

        if df_in.loc[i, "Chromosome"]=="1":
            df_in.loc[i, "Chr"]="NC_000001.11"
        elif df_in.loc[i, "Chromosome"]=="2":
            df_in.loc[i, "Chr"]="NC_000002.12"
        elif df_in.loc[i, "Chromosome"]=="3":
            df_in.loc[i, "Chr"]="NC_000003.12"
        elif df_in.loc[i, "Chromosome"]=="4":
            df_in.loc[i, "Chr"]="NC_000004.12"
        elif df_in.loc[i, "Chromosome"]=="5":
            df_in.loc[i, "Chr"]="NC_000005.10"
        elif df_in.loc[i, "Chromosome"]=="6":
            df_in.loc[i, "Chr"]="NC_000006.12"
        elif df_in.loc[i, "Chromosome"]=="7":
            df_in.loc[i, "Chr"]="NC_000007.14"
        elif df_in.loc[i, "Chromosome"]=="8":
            df_in.loc[i, "Chr"]="NC_000008.11"
        elif df_in.loc[i, "Chromosome"]=="9":
            df_in.loc[i, "Chr"]="NC_000009.12"
        elif df_in.loc[i, "Chromosome"]=="10":
            df_in.loc[i, "Chr"]="NC_000010.11"
        elif df_in.loc[i, "Chromosome"]=="11":
            df_in.loc[i, "Chr"]="NC_000011.10"
        elif df_in.loc[i, "Chromosome"]=="12":
            df_in.loc[i, "Chr"]="NC_000012.12"
        
        elif df_in.loc[i, "Chromosome"]=="13":
            df_in.loc[i, "Chr"]="NC_000013.11"
        elif df_in.loc[i, "Chromosome"]=="14":
            df_in.loc[i, "Chr"]="NC_000014.9"
        elif df_in.loc[i, "Chromosome"]=="15":
            df_in.loc[i, "Chr"]="NC_000015.10"
        elif df_in.loc[i, "Chromosome"]=="16":
            df_in.loc[i, "Chr"]="NC_000016.10"
        elif df_in.loc[i, "Chromosome"]=="17":
            df_in.loc[i, "Chr"]="NC_000017.11"
        elif df_in.loc[i, "Chromosome"]=="18":
            df_in.loc[i, "Chr"]="NC_000018.10"
        elif df_in.loc[i, "Chromosome"]=="19":
            df_in.loc[i, "Chr"]="NC_000019.10"
        elif df_in.loc[i, "Chromosome"]=="20":
            df_in.loc[i, "Chr"]="NC_000020.11"
        elif df_in.loc[i, "Chromosome"]=="21":
            df_in.loc[i, "Chr"]="NC_000021.9"
        elif df_in.loc[i, "Chromosome"]=="22":
            df_in.loc[i, "Chr"]="NC_000022.11"
            
        elif df_in.loc[i, "Chromosome"]=="X":
            df_in.loc[i, "Chr"]="NC_000023.11"
            
        elif df_in.loc[i, "Chromosome"]=="Y":
            df_in.loc[i, "Chr"]="NC_000024.10"
        else:
            df_in.loc[i, "Chr"]=df_in.loc[i, "Chromosome"]
        # if df_in.loc[i, "CHR"]==12:
        #     df_in.loc[i, "chromosome"]="NC_000012.12"
        
    return df_in


def return_gene_number(data):
    
    df_in = data.copy()
    
    for i in range(df_in.shape[0]):

        if df_in.loc[i, "chr"]=="1":
            df_in.loc[i, "chromosome"]="NC_000001.10"
        if df_in.loc[i, "chr"]=="2":
            df_in.loc[i, "chromosome"]="NC_000002.11"
        if df_in.loc[i, "chr"]=="3":
            df_in.loc[i, "chromosome"]="NC_000003.11"
        if df_in.loc[i, "chr"]=="4":
            df_in.loc[i, "chromosome"]="NC_000004.11"
        if df_in.loc[i, "chr"]=="5":
            df_in.loc[i, "chromosome"]="NC_000005.9"
        if df_in.loc[i, "chr"]=="6":
            df_in.loc[i, "chromosome"]="NC_000006.11"
        if df_in.loc[i, "chr"]=="7":
            df_in.loc[i, "chromosome"]="NC_000007.13"
        if df_in.loc[i, "chr"]=="8":
            df_in.loc[i, "chromosome"]="NC_000008.10"
        if df_in.loc[i, "chr"]=="9":
            df_in.loc[i, "chromosome"]="NC_000009.11"
        if df_in.loc[i, "chr"]=="10":
            df_in.loc[i, "chromosome"]="NC_000010.10"
        if df_in.loc[i, "chr"]=="11":
            df_in.loc[i, "chromosome"]="NC_000011.9"
        if df_in.loc[i, "chr"]=="12":
            df_in.loc[i, "chromosome"]="NC_000012.11"
        
        if df_in.loc[i, "chr"]=="13":
            df_in.loc[i, "chromosome"]="NC_000013.10"
        if df_in.loc[i, "chr"]=="14":
            df_in.loc[i, "chromosome"]="NC_000014.8"
        if df_in.loc[i, "chr"]=="15":
            df_in.loc[i, "chromosome"]="NC_000015.9"
        if df_in.loc[i, "chr"]=="16":
            df_in.loc[i, "chromosome"]="NC_000016.9"
        if df_in.loc[i, "chr"]=="17":
            df_in.loc[i, "chromosome"]="NC_000017.10"
        if df_in.loc[i, "chr"]=="18":
            df_in.loc[i, "chromosome"]="NC_000018.9"
        if df_in.loc[i, "chr"]=="19":
            df_in.loc[i, "chromosome"]="NC_000019.9"
        if df_in.loc[i, "chr"]=="20":
            df_in.loc[i, "chromosome"]="NC_000020.10"
        if df_in.loc[i, "chr"]=="21":
            df_in.loc[i, "chromosome"]="NC_000021.8"
        if df_in.loc[i, "chr"]=="22":
            df_in.loc[i, "chromosome"]="NC_000022.10"
            
        if df_in.loc[i, "chr"]=="X":
            df_in.loc[i, "chromosome"]="NC_000023.10"
            
        if df_in.loc[i, "chr"]=="Y":
            df_in.loc[i, "chromosome"]="NC_000024.9"
        # if df_in.loc[i, "CHR"]==12:
        #     df_in.loc[i, "chromosome"]="NC_000012.12"
        
    return df_in
        
        
def get_gene_positions(genbank_file):
    gene_positions = []

    with open(genbank_file) as handle:
        for record in SeqIO.parse(handle, 'genbank'):
            # print(record.features)
            for feature in record.features:
                if feature.type == 'gene':
                    gene_positions.append({
                        'chromosome': record.id,
                        'gene': feature.qualifiers.get('gene', ['Unknown'])[0],
                        'start': feature.location.start,
                        'end': feature.location.end,
                        'strand': feature.location.strand  # Add strand information
                    })
        # break

    return pd.DataFrame(gene_positions)


def genbank_to_dataframe(genbank_file):
    records = list(SeqIO.parse(genbank_file, 'genbank'))
    features = []

    for record in records:
        for feature in record.features:
            if feature.type == 'CDS' and 'gene' in feature.qualifiers:
                gene_id = feature.qualifiers['gene'][0]
                feature_dict = {
                    'seq_id': record.id,
                    'feature_type': feature.type,
                    'location': feature.location,
                    'qualifiers': feature.qualifiers,
                    'gene_id': gene_id
                }
                
                features.append(feature_dict)
                break

    return pd.DataFrame(features)
