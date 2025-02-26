#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:49:20 2023

@author: abhishekh
"""
import pandas as pd
import numpy as np

###### return one sided  mannwhitneyu test p value

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
    statistic, p_value = mannwhitneyu(group1, group2, alternative='greater', nan_policy = "omit")
    return p_value


##### return box plot in median sorted order from low to high

def box_plot(data, x, y, ax, hue = None):

    """
    Create a box plot with seaborn.

    Parameters:
    - data (DataFrame): The input DataFrame containing the data.
    - x (str): The variable to be plotted on the x-axis.
    - y (str): The variable to be plotted on the y-axis.
    - ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axes on which to draw the plot.
    - hue (str, optional): The variable that defines the subsets of the data, which will be drawn on separate facets in the plot.

    Returns:
    - ax (matplotlib.axes._subplots.AxesSubplot): The Matplotlib axes on which the box plot is drawn in median sorted order.
    """
    import seaborn as sns
    median_per_category = data.groupby(x)[y].median().sort_values()
    sns.boxplot(x=x, y=y, data=data, hue=hue, boxprops=dict(edgecolor="black", linewidth=0),ax=ax,flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, dodge=True, order=median_per_category.index, medianprops={"color": "black","linewidth":1})
    return ax

##### return axis in custom form
def custom_figure_axis(ax, fontsize=10, show_ticks = True, rotation=90):
    
        
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
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.tick_params(axis='x', labelsize=fontsize, rotation=rotation)
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


##### return_data_quantile_wise

def return_data_quantile_wise(data, column_name):
    
    import numpy as np
    
    data = data.copy()
    
    lst_quantile_dataframe = list()  # The list containing dataframe of each quantile of column_name
    lst_quantile_dataframe_column_name = list() # The list containing allele series of each quantile
    lst_quantile_dataframe_column_name_median = list() # The list containing median of allele series of each quantile
    # lst_quantile_accumulated_gain_vs_allele_mean= list()

    count = 0

    for i in range(20):

        if (i==19):
            lst_quantile_dataframe.append(data.iloc[data[column_name][(data[column_name] >= data[column_name].quantile(round((count),2))) & (data[column_name] <=data[column_name].quantile(round((count+0.05),2)))].index.tolist()])
            lst_quantile_dataframe_column_name.append(np.log10(data.iloc[data[column_name][(data[column_name] >= data[column_name].quantile(round((count),2))) & (data[column_name] <=data[column_name].quantile(round((count+0.05),2)))].index.tolist()]["allele count"]))
            lst_quantile_dataframe_column_name_median.append(np.median(np.log10(data.iloc[data[column_name][(data[column_name] >= data[column_name].quantile(round((count),2))) & (data[column_name] <=data[column_name].quantile(round((count+0.05),2)))].index.tolist()]["allele count"])))
    #         count = count + 0.05

        if (i!=19):
            lst_quantile_dataframe.append(data.iloc[data[column_name][(data[column_name] >= data[column_name].quantile(round((count),2))) & (data[column_name] <data[column_name].quantile(round((count+0.05),2)))].index.tolist()])
            lst_quantile_dataframe_column_name.append(np.log10(data.iloc[data[column_name][(data[column_name] >= data[column_name].quantile(round((count),2))) & (data[column_name] <data[column_name].quantile(round((count+0.05),2)))].index.tolist()]["allele count"]))
            lst_quantile_dataframe_column_name_median.append(np.median(np.log10(data.iloc[data[column_name][(data[column_name] >= data[column_name].quantile(round((count),2))) & (data[column_name] <data[column_name].quantile(round((count+0.05),2)))].index.tolist()]["allele count"])))

        count = count + 0.05
        
    return lst_quantile_dataframe, lst_quantile_dataframe_column_name, lst_quantile_dataframe_column_name_median
    
    


##### return a dataframe with four columns  quantile, column_name as provided, accumulate_gain
def return_quantile_wise_allele_count(lst_quantile_dataframe, column_name, ignore_label_column =False):
    
    import pandas as pd
    import numpy as np
    
    lst_quantile_dataframe = lst_quantile_dataframe.copy()
    
    quantile_range = ['0-5%', '5-10%', '10-15%', '15-20%', '20-25%', '25-30%', '30-35%', '35-40%', '40-45%', \
                       '45-50%', '50-55%', '55-60%', '60-65%', '65-70%', '70-75%', '75-80%', '80-85%', '85-90%', \
                       '90-95%', '95-100%']
    
    df_tmp_full = pd.DataFrame(columns = ["quantile of " + str(column_name), "AC"])
    for i in range(len(quantile_range)):
        df_tmp = pd.DataFrame(columns = ["quantile of " + str(column_name), "AC"])
        df_tmp["AC"] = lst_quantile_dataframe[i]["allele count"]
        df_tmp["log10(AC)"] = np.log10(lst_quantile_dataframe[i]["allele count"])
        df_tmp[column_name] = lst_quantile_dataframe[i][column_name]
        if ignore_label_column==False:
            df_tmp["label"] = lst_quantile_dataframe[i]["label"]
        df_tmp.loc[:, "quantile of " + str(column_name) ] = quantile_range[i]
        df_tmp_full = pd.concat([df_tmp_full, df_tmp])
#         print(df_now_new)
    return df_tmp_full









"""Plot three subplots for the column name vs generation number for mutator, non mutator amd all"""
def return_mutator_non_mutator_column_name_wise_graph_together(data, mutator_list, non_mutator_list, column_name, fontsize=10):
    
    
    """
    Plot three subplots for the specified column against generation number for mutator, non-mutator, and the entire population.
    
    Parameters:
        data (pd.DataFrame): The input DataFrame containing the data.
        mutator_list (list): List of labels for mutator data points.
        non_mutator_list (list): List of labels for non-mutator data points.
        column_name (str): The column to plot against generation number.
        fontsize (int, optional): Font size for plot labels. Defaults to 10.
    
    Returns:
        ax (numpy.ndarray): Array of subplot axes.
    
    Notes:
        - Requires 'label' and 'generation_number' columns in the input DataFrame.
    
    """
    
    
    
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter
    from scipy.stats import spearmanr
    import seaborn as sns
    data = data.copy();
    offset = 0

    df_non_mutator_full_mutation_data = data[data['label'].isin(non_mutator_list)]
    df_mutator_full_mutation_data = data[data['label'].isin(mutator_list)]
    
    fig, ax = plt.subplots(3, 1 , dpi = 600, figsize = (2.1, 2.1), sharex=True);
    sns.lineplot(data = df_mutator_full_mutation_data, x = "generation_number", y=column_name, color="red", lw = 0.25, estimator ="median", ax = ax[0]);
    sns.lineplot(data = df_non_mutator_full_mutation_data, x = "generation_number", y=column_name, color="green", lw = 0.25, estimator ="median", ax = ax[1]);
    sns.lineplot(data = data, x = "generation_number", y=column_name, color="blue", lw = 0.25, estimator ="median", ax = ax[2]);

    # mutator_spearman_corr = "{:.3e}".format(spearmanr(df_mutator_full_mutation_data["generation_number"], df_mutator_full_mutation_data[column_name])[0])
    # mutator_spearman_p = "{:.3e}".format(spearmanr(df_mutator_full_mutation_data["generation_number"], df_mutator_full_mutation_data[column_name])[1])
    # ax[0].text((max(ax[0].get_xlim()) -100* offset), (max(ax[0].get_ylim()) - offset), "Corr:" + str(mutator_spearman_corr) + "| p- value: " + str(mutator_spearman_p), fontsize=3, color='black', ha='right', va='top')

    # print(spearmanr(df_mutator_full_mutation_data["generation_number"], df_mutator_full_mutation_data[column_name]))
    # # formatted_corr_all = "{:.2e}".format(df_non_mutator_full_mutation_data["generation_number"].corr(df_non_mutator_full_mutation_data["column_name"]))
    # non_mutator_spearman_corr = "{:.3e}".format(spearmanr(df_non_mutator_full_mutation_data["generation_number"], df_non_mutator_full_mutation_data[column_name])[0])
    # non_mutator_spearman_p = "{:.3e}".format(spearmanr(df_non_mutator_full_mutation_data["generation_number"], df_non_mutator_full_mutation_data[column_name])[1])
    # ax[1].text((max(ax[1].get_xlim()) -100* offset), (max(ax[1].get_ylim()) - offset), "Corr:" + str(non_mutator_spearman_corr) + "| p- value: " + str(non_mutator_spearman_p), fontsize=3, color='black', ha='right', va='top')
    # print(spearmanr(df_non_mutator_full_mutation_data["generation_number"], df_non_mutator_full_mutation_data[column_name]))

    

    # all_spearman_corr = "{:.3e}".format(spearmanr(data["generation_number"], data[column_name])[0])
    # all_spearman_p = "{:.3e}".format(spearmanr(data["generation_number"], data[column_name])[1])
    # ax[2].text((max(ax[2].get_xlim()) -100* offset), (max(ax[2].get_ylim()) - offset), "Corr:" + str(all_spearman_corr) + "| p- value: " + str(all_spearman_p), fontsize=3, color='black', ha='right', va='top')

    # ax[0].set_title(column_name + " vs generation for mutator", fontsize=fontsize);
    # ax[1].set_title(column_name + " vs generation for non mutator", fontsize=fontsize);
    # ax[2].set_title(column_name + " vs generation for all population", fontsize=fontsize);
    # print(all_spearman_corr, all_spearman_p)

    # print(spearmanr(data["generation_number"], data[column_name]))

    # Apply the custom tick formatter
    formatter = FuncFormatter(format_ticks);
    ax[2].xaxis.set_major_formatter(formatter);
    for i in range(3):
        ax[i]=custom_figure_axis(ax[i], fontsize=fontsize, show_ticks = True);
        ax[i].set_ylabel(column_name)  ;
        # ax[i].legend(ncol=2, fontsize=3, frameon=False)

    plt.tight_layout();
    
    return ax




def plot_ac_vs_column(ax, data:pd.DataFrame(), lst_quantile_dataframe_column_name_median:list(), column_name:str, show_ticks:bool = True, fontsize:int=10):
    
    import seaborn as sns
    import utility
    import matplotlib.pyplot as plt
    
    data=data.copy()
    
    quantile_range = ['0-5%', '5-10%', '10-15%', '15-20%', '20-25%', '25-30%', '30-35%', '35-40%', '40-45%', \
                   '45-50%', '50-55%', '55-60%', '60-65%', '65-70%', '70-75%', '75-80%', '80-85%', '85-90%', \
                   '90-95%', '95-100%']

    
    sns.boxplot(data = data, x = "quantile of " + str(column_name), y = "log10(AC)", ax=ax[0], boxprops=dict(edgecolor="black", linewidth=0),flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, dodge=True,  medianprops={"color": "black","linewidth":0.5}, color ="lime")
    sns.lineplot(data = data, x = "quantile of " + str(column_name), y ="log10(AC)", estimator="median", lw= 0.5, color = "red", ax=ax[1])
    # sns.lineplot(data = data, x = "quantile", y = "AC", estimator="median", lw= 1, color = "lime", ax=ax[1])
    ax[1].scatter(x=quantile_range, y=lst_quantile_dataframe_column_name_median, s=0.5, color="red")


    ax[0]=utility.custom_figure_axis(ax[0], fontsize=fontsize, show_ticks = show_ticks)
    ax[1]=utility.custom_figure_axis(ax[1], fontsize=fontsize, show_ticks = show_ticks)
    ax[1].legend().set_visible(False)
    # ax[0].set_title("Ecoli AC vs quantile wise gain score", fontsize=5)
    
    ax[0].set_xlabel(None)
    # plt.subplots_adjust(hspace=0.2)
    # plt.subplots_adjust(wspace=0.2)  # You can adjust the value as needed

    plt.tight_layout()
    # plt.savefig("ecoli log10 AC vs quantile wise AG lineplot.pdf", dpi = 600, bbox_inches="tight")
    
    return ax


"""return one sided greater manwhiteney p value top 50% vs last 50% quantile""" 
    
def return_quantile_wise_pvalue(lst_quantile_dataframe:list()):

    """
    Calculate the Mann-Whitney U test p-value between the first 0-50% quantile and the rest (in log scale).

    Parameters:
    - lst_quantile_dataframe (list of DataFrames): A list of DataFrames representing quantile data.

    Returns:
    - formatted_p (str): The formatted p-value resulting from the Mann-Whitney U test.
    """
    
    import numpy as np
    import utility
    
    lst_first_50_percent_quantile = list()
    lst_last_50_percent_quantile = list()

    for i in range(len(lst_quantile_dataframe)):
        if i<10:
            lst_first_50_percent_quantile.extend(lst_quantile_dataframe[i]["allele count"])
        if i>=10:
            lst_last_50_percent_quantile.extend(lst_quantile_dataframe[i]["allele count"])

    formatted_p = "{:.2e}".format(utility.man_whiteney(np.log10(lst_last_50_percent_quantile), np.log10(lst_first_50_percent_quantile)))
    print("The p value between first 0-50% quantile and rest is : ", formatted_p)
    
    return formatted_p


""" density plot of quantile of a given column"""

def plot_non_overlapping_ag(lst_quantile_accumulated_gain_vs_allele:list(), column_name:str):
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(1, 1 , dpi = 600, figsize = (5, 5))
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[0][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[1][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[2][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[3][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[4][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[5][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[6][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[7][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[8][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[9][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[10][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[11][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[12][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[13][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[14][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[15][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[16][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[17][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[18][column_name], ax=ax)
    sns.kdeplot(lst_quantile_accumulated_gain_vs_allele[19][column_name], ax=ax)
    
    return ax

"""Convert xtick label to the 1K scale"""

def format_ticks(x:int, pos:int):
    """
    Format tick labels for a plot axis.
    
    Parameters:
        x (int): The tick value.
        pos (int): The position of the tick.
    
    Returns:
        str: The formatted tick label.
    
    Behavior:
        If x is greater than or equal to 1000, the function returns the value in 'k' format.
        Otherwise, the function returns the original value.
    
    """
    
    if x >= 1000:
        return f'{x//1000}k'
    else:
        return x
    
    
def return_combined_fitness_esm_data(fitness_dataframe, esm_dataframe):

    """
    Combines fitness data and ESM score data by calculating the median fitness and
    median ESM score for each generation. It then merges the dataframes based on the generation number.

    Parameters:
        fitness_dataframe (pd.DataFrame): DataFrame containing fitness data with columns 'Generation' and 'Fitness'.
        esm_dataframe (pd.DataFrame): DataFrame containing ESM score data with columns 'generation_number' and 'esm_score'.

    Returns:
        pd.DataFrame: Combined DataFrame with columns 'generation_number', 'median fitness', and 'median esm_score'.
    """

    fitness_median_generation_wise = fitness_dataframe.groupby(["Generation"])["Fitness"].agg(['median', 'std'])
    fitness_median_generation_wise = fitness_median_generation_wise.reset_index()
    fitness_median_generation_wise.columns = ["generation_number", "median fitness", "std fitness"]

    esm_median_generation_wise = esm_dataframe.groupby(["generation_number"])["esm_score"].agg(['median', 'std'])
    esm_median_generation_wise = esm_median_generation_wise.reset_index()
    esm_median_generation_wise.columns = ["generation_number", "median esm_score", "std esm_scor"]

    combined_fitness_esm = pd.merge(fitness_median_generation_wise, esm_median_generation_wise)

    return combined_fitness_esm



def return_combined_fitness_ag_data(fitness_dataframe, ag_dataframe):

    """
    Combines fitness data and ag score data by calculating the median fitness and
    median ag score for each generation. It then merges the dataframes based on the generation number.

    Parameters:
        fitness_dataframe (pd.DataFrame): DataFrame containing fitness data with columns 'Generation' and 'Fitness'.
        ag_dataframe (pd.DataFrame): DataFrame containing AG score data with columns 'generation_number' and 'ag_score'.

    Returns:
        pd.DataFrame: Combined DataFrame with columns 'generation_number', 'median fitness', and 'median_ag_score'.
    """

    fitness_median_generation_wise = fitness_dataframe.groupby(["Generation"])["Fitness"].agg(['median', 'std'])
    fitness_median_generation_wise = fitness_median_generation_wise.reset_index()
    fitness_median_generation_wise.columns = ["generation_number", "median fitness", "std fitness"]

    ag_median_generation_wise = ag_dataframe.groupby(["generation_number"])["AG"].agg(['median', 'std'])
    ag_median_generation_wise = ag_median_generation_wise.reset_index()
    ag_median_generation_wise.columns = ["generation_number", "median ag_score", "std ag_score"]

    combined_fitness_ag = pd.merge(fitness_median_generation_wise, ag_median_generation_wise)

    return combined_fitness_ag



def return_generation_grouped_dataframe_with_one_sided_p(data, column_name: str, gen_cut_off:int, cut_off_string1:str, cut_off_string2:str):

    """
    Group and aggregate data by generation and mutator status, and create a new column based on generation cut-offs.

    Parameters:
        data (pd.DataFrame): The input DataFrame containing the data.
        gen_cut_off (int): The generation number used as a cut-off.
        cut_off_string1 (str): Label for data points below or equal to the generation cut-off.
        cut_off_string2 (str): Label for data points above the generation cut-off.
        column_name (str): The column to perform operation.

    Returns:
        pd.DataFrame: A DataFrame with grouped and aggregated data, including the new grouping column.

    Notes:
        - The input DataFrame 'data' should have columns 'generation_number' and 'mutator'.
    """
    
    data = data.copy()
    
    score_median = data.groupby(["generation_number", "mutator"])[column_name].median()
    score_median_reset_index = score_median.reset_index()
    score_median_reset_index['group_gen'] = score_median_reset_index['generation_number'].apply(lambda x: cut_off_string1 if x <= gen_cut_off else cut_off_string2)
    
    mut_gen = score_median_reset_index[score_median_reset_index.loc[:, "mutator"]=="mutator"]
    p_mut = man_whiteney(mut_gen[mut_gen.loc[:, "group_gen"]==cut_off_string2][column_name], 
                                     mut_gen[mut_gen.loc[:, "group_gen"]==cut_off_string1][column_name])
    print("The one sided man_whiteney p value within the mutator groups for data points " + cut_off_string2 +" and " +  cut_off_string1 +" generation w.r.t " + column_name + " score is: ", str(p_mut))
    
    non_mut_gen = score_median_reset_index[score_median_reset_index.loc[:, "mutator"]=="non mutator"]
    p_non_mut = man_whiteney(non_mut_gen[non_mut_gen.loc[:, "group_gen"]==cut_off_string2][column_name], 
                                     non_mut_gen[non_mut_gen.loc[:, "group_gen"]==cut_off_string1][column_name])
    print("The one sided man_whiteney p value within the non mutator groups for data points " + cut_off_string2 + " and " +  cut_off_string1 +" w.r.t  " + column_name + " score is: ", str(p_non_mut))
        
        
    
    
    return score_median_reset_index
    
    
    
def p_value(data_new):
    
    
    df = pd.DataFrame(columns = ["Gene", "esm p", "ag p"])
    
    import scipy
    
    data_new = data_new.copy()

    
    unique_gene = list(set(data_new["Gene"]))
    
    for i in range(len(unique_gene)):
        data = data_new[data_new.loc[:, "Gene"]==unique_gene[i]]
    
        esm_score_median = data.groupby("generation_number")["esm_score"].median()
        ag_score_median = data.groupby("generation_number")["AG"].median()

        if ((len(esm_score_median[esm_score_median.index<=30000])<=15)| (len(esm_score_median[esm_score_median.index>30000])<=15)):            
            continue
        if ((len(ag_score_median[ag_score_median.index<=30000])<=15)| (len(ag_score_median[ag_score_median.index>30000])<=15)):            
            continue

        df.loc[i, "Gene"] = unique_gene[i]
        df.loc[i, "median esm"] = data["esm_score"].median()
        df.loc[i, "median AG"] = data["AG"].median()

        df.loc[i, "mean esm"] = data["esm_score"].mean()
        df.loc[i, "mean AG"] = data["AG"].mean()

        df.loc[i, "esm p"] = scipy.stats.mannwhitneyu(esm_score_median[esm_score_median.index<=30000], esm_score_median[esm_score_median.index>30000], alternative='less')[1]
        df.loc[i, "ag p"] = scipy.stats.mannwhitneyu(ag_score_median[ag_score_median.index<=30000], ag_score_median[ag_score_median.index>30000], alternative='less')[1]

    return df


def return_median_data_gen_wise(data:pd.DataFrame(), x_column: str, y_column: str):
    
    
    """
     Calculate median values of specified columns grouped by 'generation_number'.
    
     Parameters:
     - data (pd.DataFrame): Input DataFrame containing the data.
     - x_column (str): Name of the column for x-axis.
     - y_column (str): Name of the column for y-axis.
    
     Returns:
     - pd.DataFrame: DataFrame containing median values of 'x_column' and 'y_column' grouped by 'generation_number'.
    
     Note:
     - This function calculates the median values of specified columns ('x_column' and 'y_column') grouped by 'generation_number'.
     - The resulting DataFrame has columns 'generation_number', 'median x_column', and 'median y_column'.
     """


    data = data.copy()

    x_column_median = pd.DataFrame(data.groupby(["generation_number"])[x_column].median()).reset_index()
    y_column_median = pd.DataFrame(data.groupby(["generation_number"])[y_column].median()).reset_index()
    x_column_median.columns = ["generation_number", "median "+ x_column]
    y_column_median.columns = ["generation_number", "median "+ y_column]

    return pd.merge(x_column_median, y_column_median, on = "generation_number")




def spearmanr(x1, x2):

    from scipy.stats import spearmanr

    corr, p_value = spearmanr(x1, x2, nan_policy = "omit")

    spearman_corr = "{:.3e}".format(corr)
    spearman_p = "{:.3e}".format(p_value)

    return spearman_corr, spearman_p


def plot_scatter(data:pd.DataFrame(), mutator_list, non_mutator_list, x_column, y_column, ax, logy=True):
    
    """
        Plot scatter plots for mutator, non-mutator, and overall populations.
    
        Parameters:
        - data (pd.DataFrame): Input DataFrame containing the data.
        - mutator_list (list): List of labels for mutator population.
        - non_mutator_list (list): List of labels for non-mutator population.
        - x_column (str): Name of the column for the x-axis.
        - y_column (str): Name of the column for the y-axis.
        - logy (bool, optional): If True, apply log10 transformation to y-axis. Default is True.
        
        Returns:
        - None
        
        Note:
        - The function generates scatter plots for mutator, non-mutator, and overall populations
          based on the specified x_column and y_column.

    """
    
    import matplotlib.pyplot as plt
    import seaborn as sns

    data = data.copy()

    all_population_combined_median  = return_median_data_gen_wise(data, x_column=x_column, y_column=y_column)
    non_mutator_population_combined_median = return_median_data_gen_wise(data[data["label"].isin(non_mutator_list)], x_column=x_column, y_column=y_column)
    mutator_population_combined_median  = return_median_data_gen_wise(data[data["label"].isin(mutator_list)], x_column=x_column, y_column=y_column)

    norm = plt.Normalize(all_population_combined_median['generation_number'].min(), all_population_combined_median['generation_number'].max())
    sm = plt.cm.ScalarMappable(cmap="mako", norm=norm)
    y_updated = "median "+ y_column
    x_updated = "median "+ x_column

    if logy:
        y_updated = "log10(median " + y_column + ")"
        mutator_population_combined_median[y_updated] = np.log10(mutator_population_combined_median["median "+ y_column])
        non_mutator_population_combined_median[y_updated] = np.log10(non_mutator_population_combined_median["median "+ y_column])
        all_population_combined_median[y_updated] = np.log10(all_population_combined_median["median "+ y_column])



    ax[0] = sns.scatterplot(mutator_population_combined_median, x = x_updated, y = y_updated, hue="generation_number",ax= ax[0], legend=True, s = 5, palette='mako')
    ax[1] = sns.scatterplot(non_mutator_population_combined_median, x = x_updated, y =y_updated, hue="generation_number",ax= ax[1], legend=False, s = 5, palette='mako')
    ax[2] = sns.scatterplot(all_population_combined_median, x = x_updated, y = y_updated, hue="generation_number",ax= ax[2], legend=False, s = 5, palette='mako')

    offset = 0.01


    mutator_spearman_corr, mutator_spearman_p = spearmanr(mutator_population_combined_median[x_updated], mutator_population_combined_median[y_updated])
    non_mutator_spearman_corr, non_mutator_spearman_p = spearmanr(non_mutator_population_combined_median[x_updated], non_mutator_population_combined_median[y_updated])
    combined_spearman_corr, combined_spearman_p = spearmanr(all_population_combined_median[x_updated], all_population_combined_median[y_updated])

    # print(mutator_spearman_corr, mutator_spearman_p )

    ax[0].text((max(ax[0].get_xlim()) -100* offset), (max(ax[0].get_ylim()) - offset), "Corr:" + str(mutator_spearman_corr) + "| p- value: " + str(mutator_spearman_p), fontsize=3, color='red', ha='right', va='top')
    ax[1].text((max(ax[1].get_xlim()) -100* offset), (max(ax[1].get_ylim()) - offset), "Corr:" + str(non_mutator_spearman_corr) + "| p- value: " + str(non_mutator_spearman_p), fontsize=3, color='red', ha='right', va='top')
    ax[2].text((max(ax[2].get_xlim()) -100* offset), (max(ax[2].get_ylim()) - offset), "Corr:" + str(combined_spearman_corr) + "| p- value: " + str(combined_spearman_p), fontsize=3, color='red', ha='right', va='top')

    # formatted_correlation_mutator = "{:.2e}".format(mutator_population_combined_median[x_updated].corr(mutator_population_combined_median[y_updated]))
    # ax[0].text((max(ax[0].get_xlim()) - 100*offset), (max(ax[0].get_ylim()) - offset), "Corr:" + formatted_correlation_mutator, fontsize=5, color='red', ha='right', va='top')

    # formatted_correlation_non_mutator = "{:.2e}".format(non_mutator_population_combined_median[x_updated].corr(non_mutator_population_combined_median[y_updated]))
    # ax[1].text((max(ax[1].get_xlim()) - 100*offset), (max(ax[1].get_ylim()) - offset), "Corr:" + formatted_correlation_non_mutator, fontsize=5, color='red', ha='right', va='top')

    # formatted_corr_all = "{:.2e}".format(all_population_combined_median[x_updated].corr(all_population_combined_median[y_updated]))
    # ax[2].text((max(ax[2].get_xlim()) -100* offset), (max(ax[2].get_ylim()) - offset), "Corr:" + formatted_corr_all, fontsize=5, color='red', ha='right', va='top')


    ax[0].set_title("mutator", fontsize=5)
    ax[1].set_title("non mutator", fontsize=5)
    ax[2].set_title("Overall population", fontsize=5)

    for i in range(3):
        ax[i] = custom_figure_axis(ax[i], fontsize=5, show_ticks = True, rotation=0)
        # ax[i].set_ylabel("log10(Median AC)")

    ax[0].get_legend().remove()
    # cbar = ax[0].figure.colorbar(sm, ax=ax[2])
    # ax[0].figure.colorbar(sm)
    # ax[1].figure.colorbar(sm)
    # ax[2].figure.colorbar(sm)


    plt.suptitle(y_updated + " vs " + x_updated + " scatter plots for all type", fontsize=6)

    plt.tight_layout()


    return ax#mutator_population_combined_median, non_mutator_population_combined_median, all_population_combined_median
    
    
    
    
def dba_stat(data, group_by_col,target_col_name, num_permutations, random_state = 42):


    """
        Performs a permutation test using Kendall's Tau to evaluate the correlation 
        between a median wise target column and a grouping variable.

        Parameters:
        - data (pd.DataFrame): Input DataFrame containing the grouping column and target column.
        - group_by_col (str): Column name to group the data by.
        - target_col_name (str): Column name of the target variable.
        - num_permutations (int): Number of permutations to perform (default=1000).
        - random_state (int): Random seed for reproducibility (default=42).

        Returns:
        - float: Permutation-based p-value.
        - float: Observed Kendall's Tau correlation.
    """
    data = data.copy()

    import numpy as np
    from scipy.stats import kendalltau

    np.random.seed(random_state)

    def return_kendall_corr(x, y):
        kendall_corr, kendall_p = kendalltau(x, y)
        return kendall_corr
    
    f = return_kendall_corr
    
    generation_medians_original = data.groupby(group_by_col).agg({
            target_col_name: 'median',    
        }).reset_index()
    
    # print(generation_medians_original)

    d_obs = f(np.array(range(len(generation_medians_original[target_col_name]))), generation_medians_original[target_col_name])

    original_scores = data[target_col_name].values

    # print(original_scores.shape)
    pseudo_corr = []

    for i in range(1, num_permutations + 1):
        random_column_name = f'random_{i}_' + target_col_name
        data.loc[:, random_column_name] = np.random.permutation(original_scores)
        # print(data.loc[:, random_column_name])
        
        generation_medians_random = data.groupby(group_by_col).agg({
            random_column_name: 'median'}).reset_index()

        # print(generation_medians_random)

        pseudo_corr.append(f(np.array(range(len(generation_medians_random[random_column_name]))), generation_medians_random[random_column_name]))

    pseudo_corr = np.array(pseudo_corr)

    freq_p_value = (np.sum(pseudo_corr >= d_obs)+1) / (num_permutations + 1)

    return freq_p_value, d_obs


def return_box_with_p_effect_size(ax, data, x_columns, y_columns, group1_name, group2_name, y_pos, hue_columns=None, fontsize=5):


    import pingouin as pg
    
    # Extract x-tick labels and their positions from the Axes object
    xtick_labels = [tick.get_text() for tick in ax.get_xticklabels()]
    # print(xtick_labels)
    xtick_positions = [tick.get_position()[0] for tick in ax.get_xticklabels()]

    try:

        if not hue_columns is None:
            

            for idx, (label, pos) in enumerate(zip(xtick_labels, xtick_positions)):
                # Filter data for the current x-tick category
                data_tmp = data[data[x_columns] == label]

                

                
                # print("hello")
                # print(data_tmp[data_tmp[hue_columns] == group1_name])
                # Separate groups based on hue column
                group1 = data_tmp[data_tmp[hue_columns] == group1_name][y_columns]
                # print(group1)
                group2 = data_tmp[data_tmp[hue_columns] == group2_name][y_columns]

                # print("hi")

                # Compute p-value using Mann-Whitney U test
                p_value = man_whiteney(group1, group2)
                

                # Compute effect size using Cohen's d
                d_value = pg.compute_effsize(group1, group2, eftype="cohen")

                # Annotate the plot with p-value and effect size
                ax.text(pos, y_pos, f"p: {p_value:.2e} \nd: {d_value:.2e}", ha="center", fontsize=fontsize, color="red")

                print(f"x: {label}, p-value: {p_value:.2e}, Cohen's d: {d_value:.2e}")

            return ax
            
        if ((hue_columns is None) and (len(xtick_labels) == 2)):      
            # for idx, (label, pos) in enumerate(zip(xtick_labels, xtick_positions)): 
            group1 = data[data[x_columns] == group1_name][y_columns]
            group2 = data[data[x_columns] == group2_name][y_columns]

            
     
            # Compute p-value using Mann-Whitney U test
            p_value = man_whiteney(group1, group2)

            # Compute effect size using Cohen's d
            d_value = pg.compute_effsize(group1, group2, eftype="cohen")

            # Annotate the plot with p-value and effect size
            ax.text(0.5, y_pos, f"p: {p_value:.2e} \nd: {d_value:.2e}", ha="center", fontsize=fontsize, color="red")

            # print(f"x: {label}, p-value: {p_value:.2e}, Cohen's d: {d_value:.2e}")
            return ax

    except:

        print("There is some error")


