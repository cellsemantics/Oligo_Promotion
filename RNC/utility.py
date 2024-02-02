#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:49:20 2023

@author: abhishekh
"""
import pandas as pd
import re
import numpy as np

def return_groupwise_amino_acid(filepath: str, sheet_name: str):
    
    
    """
    The function takes a file path and sheet name to extract score data, and outputs a DataFrame where 
    rows represent mutant codons, columns indicate mutation positions, and cell values indicate the 
    mean score for each mutation at the respective position and amino acid, after removing the stop codon.
    
    Args:
        filepath (str): Path to the metasheet.
        sheet_name (str): The name of the sheet in the metasheet where the score of interest is located.  
    
    Returns:
        pandas.DataFrame: A DataFrame where each row index indicates the mutant amino acid, each column 
        index indicates the position of mutation, and the value in each cell indicates the mean score 
        value for the mutation at a position with the mutant amino acid present at the DataFrame row index. 
    """
    
    df = pd.read_excel(filepath, sheet_name = sheet_name) # Read the score sheet in dataframe

    df = df[df.loc[:, "aa"]!="Stop"] # exclude the stop codon 
    original_order = df['aa'].unique()
    df['aa'] = pd.Categorical(df['aa'], categories=original_order, ordered=True)
    grouped_amino_acid= df.iloc[:, 2:].groupby('aa').mean().reset_index() 
    grouped_amino_acid.index = grouped_amino_acid["aa"]
    grouped_amino_acid.drop(["aa"], inplace=True, axis=1)
    grouped_amino_acid.columns = range(len(grouped_amino_acid.columns))
    return grouped_amino_acid

def pre_process_score(data:pd.DataFrame()):
    
    """
    Processes fitness scores data.
    
    Parameters:
        data (pandas DataFrame): Contains fitness scores with columns 'Amino acid', 'Codon', and 'aa'.
    
    Returns:
        df_tmp (pandas DataFrame): Contains all point mutations.
        complete_amino_acid (str): Represents the reference amino acid sequence after cleaning.
    """
    
    data = data.copy()
    data = data.iloc[:, 2:]  # drop column 'Amino acid', 'Codon'
    #### remove any stop codon from data
    data = data[data.loc[:, "aa"]!="Stop"]
    #### group the fitness w.r.t aa (amino acid) by mean value
    grouped_amino_acid = data.groupby('aa').mean().reset_index()
    #### get the list of amino acid in the ref sequence
    ref_seq_amino_acid = grouped_amino_acid.columns.to_list()
    ref_seq_amino_acid.remove("aa")


    lst = list()
    df_tmp = pd.DataFrame(columns=["mutant"])

    def remove_numbers_special_chars(s):
        return re.sub('[^A-Za-z]', '', s)  # Keep only alphabetic characters


    # Remove numbers and special characters from each string
    cleaned_list = [remove_numbers_special_chars(s) for s in ref_seq_amino_acid]
    
    complete_amino_acid = ''.join(cleaned_list) # this is the actual reference amino acid

    mutated_codon_list = list(np.unique(grouped_amino_acid["aa"])) # list of mutated codon


    for i in range(len(cleaned_list)):
        for j in range(len(mutated_codon_list)):

            lst.append(cleaned_list[i]+ str(i) + mutated_codon_list[j])
            
    df_tmp["mutant"] = lst
#     df_tmp["score"] = np.array(grouped_amino_acid.iloc[:, 1:]).T.ravel()
    return df_tmp, complete_amino_acid # #df_tmp : mutant dataframe and complete_amino_acid: actual referce amino acid sequence

def return_ztransformed_dataframe(data:pd.DataFrame()):

    """
    Z-Transforms a pandas dataframe for enhanced comparability and comprehension.
    
    Parameters:
        data (pandas dataframe): The dataframe to be transformed.

    Returns:
        z_transform (pandas dataframe): The transformed dataframe after applying the Z-transformation.
    """

    data = data.copy()
    z_transform = (data - np.nanmean(np.array(data)))/np.nanstd(np.array(data))
    return z_transform

def custom_figure_axis(ax, fontsize=10, show_ticks = True):

    """
    Customize the appearance of matplotlib axis for a figure.

    Parameters:
        ax (matplotlib.axes._subplots.AxesSubplot): The axis to be customized.
        fontsize (int, optional): Font size for axis labels and ticks. Default is 10.
        show_ticks (bool, optional): Whether to display ticks and labels. Default is True.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: The customized axis.
    """

    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['top'].set_color('white')
    ax.spines['right'].set_color('white')
    ax.grid(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(0.2)
    ax.spines['left'].set_linewidth(0.2)
    ax.tick_params(axis='x', labelsize=fontsize, rotation=90)
    ax.tick_params(axis='y', labelsize=fontsize)
    ax.tick_params(axis='both', which='both', width=0.5)
    ax.xaxis.label.set_fontsize(fontsize)
    ax.yaxis.label.set_fontsize(fontsize)
    
    if show_ticks==False:
        ax.set_xlabel(None)
        ax.set_ylabel(None)
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        
    return ax

def man_whiteney(group1, group2):   

    """
    Perform a one-sided Mann-Whitney U test and return the p-value.

    This function compares two independent groups and tests if 'group1' tends to have larger values than 'group2'.

    Parameters:
        group1 (array-like): The data of the first group.
        group2 (array-like): The data of the second group.

    Returns:
        float: The p-value from the one-sided Mann-Whitney U test.

    Notes:
        The null hypothesis is that the distribution of 'group1' is not greater than 'group2'.
    """


    from scipy.stats import mannwhitneyu
    statistic, p_value = mannwhitneyu(group1, group2, alternative='greater')
    return p_value



def box_plot_esm_vs_all_three_score(data: pd.DataFrame(), ax):

    """
    Generate 3 custom box plots based on median cutoffs for fitness, functional score, and transformation score.

    Parameters:
        data (pd.DataFrame): Input dataframe containing various scores.
        ax (matplotlib.axes._subplots.AxesSubplot): The axis for plotting.

    Returns:
        matplotlib.axes._subplots.AxesSubplot: The customized axis after plotting.
    """

    import seaborn as sns
    import scipy.stats as stats

    data = data.copy()

    mean_funtional_cutoff = np.round(data.loc[:, "mean_funtional"].median(), 3)
    weighted_mean_fitness_cutoff =    np.round(data.loc[:, "weighted_mean_fitness"].median(), 3)
    transformation_score_cutoff =        np.round(data.loc[:, "transformation score"].median(), 3)

    lst_mean_funtional = list()
    lst_weighted_mean_fitness = list()
    lst_transformation_score = list()

    lst_weighted_mean_fitness.append(data[data.loc[:, "weighted_mean_fitness"]<=weighted_mean_fitness_cutoff]["esm"])
    lst_weighted_mean_fitness.append(data[data.loc[:, "weighted_mean_fitness"]>weighted_mean_fitness_cutoff]["esm"])

    lst_mean_funtional.append(data[data.loc[:, "mean_funtional"]<=mean_funtional_cutoff]["esm"])
    lst_mean_funtional.append(data[data.loc[:, "mean_funtional"]>mean_funtional_cutoff]["esm"])

    lst_transformation_score.append(data[data.loc[:,  'transformation score']<=transformation_score_cutoff]["esm"])
    lst_transformation_score.append(data[data.loc[:,  'transformation score']>transformation_score_cutoff]["esm"])

    sns.boxplot(lst_weighted_mean_fitness, ax=ax[0], color = "lime", boxprops=dict(edgecolor="black", linewidth=0),flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, medianprops={"color": "black","linewidth":1})
    sns.boxplot(lst_mean_funtional, ax=ax[1], color = "lime", boxprops=dict(edgecolor="black", linewidth=0),flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, medianprops={"color": "black","linewidth":1})
    sns.boxplot(lst_transformation_score, ax=ax[2], color = "lime", boxprops=dict(edgecolor="black", linewidth=0),flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, medianprops={"color": "black","linewidth":1})
    ax[0].set_ylabel("esm score", fontsize=5)
    ax[0].set_xlabel("weighted mean fitness")
    ax[1].set_xlabel("mean functional score")
    ax[2].set_xlabel("Transformation score")

    ax[0].set_xticklabels(["<=" + str(weighted_mean_fitness_cutoff), ">" + str(weighted_mean_fitness_cutoff)], fontsize=5)
    ax[1].set_xticklabels(["<=" + str(mean_funtional_cutoff), ">" + str(mean_funtional_cutoff)], fontsize=5)
    ax[2].set_xticklabels(["<=" + str(transformation_score_cutoff), ">" + str(transformation_score_cutoff)], fontsize=5)

    offset = 0.01

    formatted_p = "{:.2e}".format(man_whiteney(lst_weighted_mean_fitness[1], lst_weighted_mean_fitness[0]))
    ax[0].text((max(ax[0].get_xlim()) - offset), (max(ax[0].get_ylim()) - offset), "One sided p:" + formatted_p, fontsize=6, color='red', ha='right', va='top')
    # print(formatted_p)
    formatted_p = "{:.2e}".format(man_whiteney(lst_mean_funtional[1], lst_mean_funtional[0]))
    ax[1].text((max(ax[1].get_xlim()) - offset), (max(ax[1].get_ylim()) - offset), "One sided p:" + formatted_p, fontsize=6, color='red', ha='right', va='top')
    # print(formatted_p)

    formatted_p = "{:.2e}".format(man_whiteney(lst_transformation_score[1], lst_transformation_score[0]))
    ax[2].text((max(ax[2].get_xlim()) - offset), (max(ax[2].get_ylim()) - offset), "One sided p:" + formatted_p, fontsize=6, color='red', ha='right', va='top')

    for i in range(3):
        ax[i] =  custom_figure_axis(ax[i], fontsize=5, show_ticks = True)

    return ax

