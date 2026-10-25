#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 26 13:49:20 2023

@author: abhishekh
"""
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

def apply_figure_style(fontsize=6, font_family="Arial", figure_dpi=150):
    """
    Apply a consistent global Matplotlib/Seaborn style across figures.

    Parameters:
        fontsize (int, optional): Default font size for figure text.
        font_family (str, optional): Preferred sans-serif font family.
        figure_dpi (int, optional): Default DPI for generated figures.
    """
    import logging
    import matplotlib as mpl
    import seaborn as sns
    from matplotlib import font_manager

    # Silence repetitive font lookup warnings when a requested font is unavailable.
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

    available_fonts = {f.name for f in font_manager.fontManager.ttflist}
    fallback_chain = [
        font_family,
        "Arial",
        "Liberation Sans",
        "DejaVu Sans",
        "Noto Sans",
    ]
    resolved_fonts = []
    seen = set()
    for name in fallback_chain:
        if name not in seen and name in available_fonts:
            resolved_fonts.append(name)
            seen.add(name)

    if not resolved_fonts:
        resolved_fonts = ["DejaVu Sans"]

    style_rc = {
        "font.family": "sans-serif",
        "font.sans-serif": resolved_fonts,
        "font.size": fontsize,
        "axes.titlesize": fontsize,
        "axes.labelsize": fontsize,
        "xtick.labelsize": fontsize,
        "ytick.labelsize": fontsize,
        "legend.fontsize": fontsize,
        "figure.titlesize": fontsize,
        "axes.linewidth": 1.2,
        "figure.dpi": figure_dpi,
    }
    mpl.rcParams.update(style_rc)
    # Pass RC explicitly so later seaborn styling calls preserve this baseline.
    sns.set_theme(style="white", rc=style_rc)

###### return one sided  mannwhitneyu test p value

def man_whiteney(group1, group2, alternative="greater"):   
    
    """
    Perform a Mann-Whitney U test and return the p-value.
    
    This function compares two independent groups with a configurable alternative.
    
    Parameters:
        group1 (array-like): The data of the first group.
        group2 (array-like): The data of the second group.
    
    Parameters:
        alternative (str, optional): One of `greater`, `less`, or `two-sided`.
            Default is `greater`.

    Returns:
        float: The p-value from the one-sided Mann-Whitney U test.
    
    Notes:
        The null hypothesis depends on `alternative`.
    """

    from scipy.stats import mannwhitneyu
    statistic, p_value = mannwhitneyu(group1, group2, alternative=alternative, nan_policy = "omit")
    return p_value
    


def median_based_distance(group1, group2):
    """

    
    Parameters:
        group1: list or numpy array, first sample
        group2: list or numpy array, second sample
    
    Returns:
        d_median: float, effect size based on median and MAD
    """

    try:

        import numpy as np
        from scipy.stats import median_abs_deviation
        # Compute medians
        median1 = pd.Series(group1).median()
        median2 = pd.Series(group2).median()


        
        # Compute MAD (median absolute deviation)
        mad1 = median_abs_deviation(group1, nan_policy = "omit")  # Scale for consistency with std
        mad2 = median_abs_deviation(group2, nan_policy = "omit")



        # Compute pooled MAD (use NaN-aware counts to match nan_policy="omit" in MAD)
        n1 = pd.Series(group1).count()
        n2 = pd.Series(group2).count()
        pooled_mad = np.sqrt(((n1 - 1) * mad1**2 + (n2 - 1) * mad2**2) / (n1 + n2 - 2))
        
        # Compute Cohen's d using median and MAD
        d_median = (median1 - median2) / pooled_mad if pooled_mad > 0 else np.nan
        
        return d_median

    except:
        return np.nan


def median_based_distance_ci(group1, group2, n_boot=2000, ci=95, random_state=42):
    """
    Returns:
        effect_size: observed median-based distance
        ci_low: lower CI bound
        ci_high: upper CI bound
    """
    import numpy as np
    import pandas as pd

    g1 = pd.Series(group1).dropna().to_numpy()
    g2 = pd.Series(group2).dropna().to_numpy()

    if len(g1) < 2 or len(g2) < 2:
        return np.nan, np.nan, np.nan

    observed = median_based_distance(g1, g2)
    if np.isnan(observed):
        return np.nan, np.nan, np.nan

    rng = np.random.default_rng(random_state)
    boot_stats = []

    for _ in range(n_boot):
        boot_g1 = rng.choice(g1, size=len(g1), replace=True)
        boot_g2 = rng.choice(g2, size=len(g2), replace=True)

        stat = median_based_distance(boot_g1, boot_g2)
        if not np.isnan(stat):
            boot_stats.append(stat)

    if len(boot_stats) == 0:
        return observed, np.nan, np.nan

    alpha = (100 - ci) / 2
    ci_low = np.percentile(boot_stats, alpha)
    ci_high = np.percentile(boot_stats, 100 - alpha)

    return observed, ci_low, ci_high



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
    """
    Split a dataframe into 20 quantile bins of `column_name` and summarize allele counts.

    Parameters:
        data (pd.DataFrame): Input dataframe containing `column_name` and `allele count`.
        column_name (str): Column used to compute quantile bins.

    Returns:
        tuple:
            - list[pd.DataFrame]: Dataframes for each quantile slice.
            - list[pd.Series]: log10 allele-count series for each quantile slice.
            - list[float]: Median log10 allele count for each quantile slice.
    """
    
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
    """
    Build a long-format dataframe of allele counts grouped by quantile labels.

    Parameters:
        lst_quantile_dataframe (list[pd.DataFrame]): Quantile-split dataframes.
        column_name (str): Feature used to generate quantiles.
        ignore_label_column (bool, optional): If True, omit the `label` column.

    Returns:
        pd.DataFrame: Combined dataframe with quantile label, AC, log10(AC),
        feature values, and optionally label.
    """
    
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




def plot_ac_vs_column(ax, data: pd.DataFrame, lst_quantile_dataframe_column_name_median: list, column_name: str, show_ticks: bool = True, fontsize: int = 10):
    """
    Plot quantile-wise allele count distributions and median trends.

    Parameters:
        ax (list[matplotlib.axes.Axes]): Two axes for boxplot and lineplot.
        data (pd.DataFrame): Long-format quantile dataframe from
            `return_quantile_wise_allele_count`.
        lst_quantile_dataframe_column_name_median (list): Median log10(AC) per quantile.
        column_name (str): Feature name used in quantile labels.
        show_ticks (bool, optional): Whether to show tick labels.
        fontsize (int, optional): Axis font size.

    Returns:
        list[matplotlib.axes.Axes]: Updated axes.
    """
    
    import seaborn as sns
    import matplotlib.pyplot as plt
    
    data=data.copy()
    
    quantile_range = ['0-5%', '5-10%', '10-15%', '15-20%', '20-25%', '25-30%', '30-35%', '35-40%', '40-45%', \
                   '45-50%', '50-55%', '55-60%', '60-65%', '65-70%', '70-75%', '75-80%', '80-85%', '85-90%', \
                   '90-95%', '95-100%']

    
    sns.boxplot(data = data, x = "quantile of " + str(column_name), y = "log10(AC)", ax=ax[0], boxprops=dict(edgecolor="black", linewidth=0),flierprops=dict(marker='o', markeredgecolor='black'), width=0.4,linewidth=0.5, fliersize=0.15, dodge=True,  medianprops={"color": "black","linewidth":0.5}, color ="lime")
    sns.lineplot(data = data, x = "quantile of " + str(column_name), y ="log10(AC)", estimator="median", lw= 0.5, color = "red", ax=ax[1])
    # sns.lineplot(data = data, x = "quantile", y = "AC", estimator="median", lw= 1, color = "lime", ax=ax[1])
    ax[1].scatter(x=quantile_range, y=lst_quantile_dataframe_column_name_median, s=0.5, color="red")


    ax[0]=custom_figure_axis(ax[0], fontsize=fontsize, show_ticks = show_ticks)
    ax[1]=custom_figure_axis(ax[1], fontsize=fontsize, show_ticks = show_ticks)
    ax[1].legend().set_visible(False)
    # ax[0].set_title("Ecoli AC vs quantile wise gain score", fontsize=5)
    
    ax[0].set_xlabel(None)
    # plt.subplots_adjust(hspace=0.2)
    # plt.subplots_adjust(wspace=0.2)  # You can adjust the value as needed

    plt.tight_layout()
    # plt.savefig("ecoli log10 AC vs quantile wise AG lineplot.pdf", dpi = 600, bbox_inches="tight")
    
    return ax


"""return one sided greater manwhiteney p value top 50% vs last 50% quantile""" 
    
def return_quantile_wise_pvalue(lst_quantile_dataframe: list):

    """
    Calculate the Mann-Whitney U test p-value between the first 0-50% quantile and the rest (in log scale).

    Parameters:
    - lst_quantile_dataframe (list of DataFrames): A list of DataFrames representing quantile data.

    Returns:
    - formatted_p (str): The formatted p-value resulting from the Mann-Whitney U test.
    """
    
    import numpy as np
    
    lst_first_50_percent_quantile = list()
    lst_last_50_percent_quantile = list()

    for i in range(len(lst_quantile_dataframe)):
        if i<10:
            lst_first_50_percent_quantile.extend(lst_quantile_dataframe[i]["allele count"])
        if i>=10:
            lst_last_50_percent_quantile.extend(lst_quantile_dataframe[i]["allele count"])

    formatted_p = "{:.2e}".format(man_whiteney(np.log10(lst_last_50_percent_quantile), np.log10(lst_first_50_percent_quantile)))
    print("The p value between first 0-50% quantile and rest is : ", formatted_p)
    
    return formatted_p


""" density plot of quantile of a given column"""

def plot_non_overlapping_ag(lst_quantile_accumulated_gain_vs_allele:list(), column_name:str):
    """
    Draw KDE curves for `column_name` across the 20 quantile dataframes.

    Parameters:
        lst_quantile_accumulated_gain_vs_allele (list[pd.DataFrame]): Quantile-split data.
        column_name (str): Column to visualize with KDE.

    Returns:
        matplotlib.axes.Axes: Axis containing overlaid KDE plots.
    """
    
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
        return str(x)
    
    
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
    esm_median_generation_wise.columns = ["generation_number", "median esm_score", "std esm_score"]

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
    score_median_reset_index['group_gen'] = score_median_reset_index['generation_number'].apply(
        lambda x: cut_off_string1 if x <= gen_cut_off else cut_off_string2
    )
    score_median_reset_index["mutator_norm"] = (
        score_median_reset_index["mutator"]
        .astype(str)
        .str.strip()
        .str.lower()
        .str.replace(" ", "_", regex=False)
    )
    
    mut_gen = score_median_reset_index[score_median_reset_index.loc[:, "mutator_norm"]=="mutator"]
    mut_late = mut_gen[mut_gen.loc[:, "group_gen"]==cut_off_string2][column_name]
    mut_early = mut_gen[mut_gen.loc[:, "group_gen"]==cut_off_string1][column_name]
    if len(mut_late) > 0 and len(mut_early) > 0:
        p_mut = man_whiteney(mut_late, mut_early)
        print(
            "The one sided man_whiteney p value within the mutator groups for data points "
            + cut_off_string2 + " and " + cut_off_string1
            + " generation w.r.t " + column_name + " score is: ", str(p_mut)
        )
    else:
        print("Insufficient mutator-group data for one-sided Mann-Whitney test.")
    
    non_mut_gen = score_median_reset_index[score_median_reset_index.loc[:, "mutator_norm"]=="non_mutator"]
    non_mut_late = non_mut_gen[non_mut_gen.loc[:, "group_gen"]==cut_off_string2][column_name]
    non_mut_early = non_mut_gen[non_mut_gen.loc[:, "group_gen"]==cut_off_string1][column_name]
    if len(non_mut_late) > 0 and len(non_mut_early) > 0:
        p_non_mut = man_whiteney(non_mut_late, non_mut_early)
        print(
            "The one sided man_whiteney p value within the non mutator groups for data points "
            + cut_off_string2 + " and " + cut_off_string1
            + " w.r.t  " + column_name + " score is: ", str(p_non_mut)
        )
    else:
        print("Insufficient non-mutator-group data for one-sided Mann-Whitney test.")
        
        
    
    
    return score_median_reset_index
    
    
    
def p_value(data_new):
    """
    Compute per-gene one-sided Mann-Whitney p-values for pre/post generation cutoffs.

    For each gene, this function compares median generation-wise `esm_score` and `AG`
    values before and after generation 30000, and reports summary statistics.

    Parameters:
        data_new (pd.DataFrame): Input dataframe with columns including
            `Gene`, `generation_number`, `esm_score`, and `AG`.

    Returns:
        pd.DataFrame: Per-gene table containing median/mean metrics and p-values.
    """
    
    
    df = pd.DataFrame(columns = ["Gene", "esm p", "ag p"])
    
    import scipy
    
    data_new = data_new.copy()

    
    unique_gene = list(set(data_new["Gene"]))
    
    for i in range(len(unique_gene)):
        data = data_new[data_new.loc[:, "Gene"]==unique_gene[i]]
    
        esm_score_median = data.groupby("generation_number")["esm_score"].median()
        ag_score_median = data.groupby("generation_number")["AG"].median()

        if ((len(esm_score_median[esm_score_median.index<=30000])<=15) or (len(esm_score_median[esm_score_median.index>30000])<=15)):
            continue
        if ((len(ag_score_median[ag_score_median.index<=30000])<=15) or (len(ag_score_median[ag_score_median.index>30000])<=15)):
            continue

        df.loc[i, "Gene"] = unique_gene[i]
        df.loc[i, "median esm"] = data["esm_score"].median()
        df.loc[i, "median AG"] = data["AG"].median()

        df.loc[i, "mean esm"] = data["esm_score"].mean()
        df.loc[i, "mean AG"] = data["AG"].mean()

        df.loc[i, "esm p"] = scipy.stats.mannwhitneyu(esm_score_median[esm_score_median.index<=30000], esm_score_median[esm_score_median.index>30000], alternative='less')[1]
        df.loc[i, "ag p"] = scipy.stats.mannwhitneyu(ag_score_median[ag_score_median.index<=30000], ag_score_median[ag_score_median.index>30000], alternative='less')[1]

    return df


def return_median_data_gen_wise(data: pd.DataFrame, x_column: str, y_column: str):
    
    
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




def spearmanr_formatted(x1, x2):
    """
    Compute Spearman correlation and p-value with NaN handling.

    Parameters:
        x1 (array-like): First variable.
        x2 (array-like): Second variable.

    Returns:
        tuple[str, str]: Correlation and p-value formatted in scientific notation.
    """

    from scipy.stats import spearmanr as _spearmanr

    corr, p_value = _spearmanr(x1, x2, nan_policy = "omit")

    spearman_corr = "{:.3e}".format(corr)
    spearman_p = "{:.3e}".format(p_value)

    return spearman_corr, spearman_p


def plot_scatter(data: pd.DataFrame, mutator_list, non_mutator_list, x_column, y_column, ax, logy=True):
    
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

    rng = np.random.default_rng(random_state)

    def return_kendall_corr(x, y):
        kendall_corr, kendall_p = kendalltau(x, y)
        return kendall_corr
    
    f = return_kendall_corr
    
    generation_medians_original = data.groupby(group_by_col).agg({
            target_col_name: 'median',    
        }).reset_index()
    
    # print(generation_medians_original)

    d_obs = f(np.array(range(len(generation_medians_original[target_col_name]))), generation_medians_original[target_col_name])

    group_labels = data[group_by_col].values
    original_scores = data[target_col_name].values

    # print(original_scores.shape)
    pseudo_corr = []

    for _ in range(num_permutations):
        permuted_scores = rng.permutation(original_scores)
        generation_medians_random = pd.DataFrame(
            {group_by_col: group_labels, "permuted_score": permuted_scores}
        ).groupby(group_by_col)["permuted_score"].median().reset_index()

        pseudo_corr.append(
            f(
                np.array(range(len(generation_medians_random["permuted_score"]))),
                generation_medians_random["permuted_score"],
            )
        )

    pseudo_corr = np.array(pseudo_corr)

    freq_p_value = (np.sum(pseudo_corr >= d_obs)+1) / (num_permutations + 1)

    return freq_p_value, d_obs



def get_significance_star(p):
    """
    Map a p-value to conventional significance-star notation.

    Parameters:
        p (float): P-value.

    Returns:
        str: One of `****`, `***`, `**`, `*`, or `ns`.
    """
    if p < 0.0001:
        return '****'
    elif p < 0.001:
        return '***'
    elif p < 0.01:
        return '**'
    elif p < 0.05:
        return '*'
    else:
        return 'ns'



def return_box_with_p_effect_size(
    ax,
    data,
    x_columns,
    y_columns,
    group1_name,
    group2_name,
    y_pos,
    hue_columns=None,
    fontsize=5,
    p_adjust_method=None,
    alternative="greater",
):
    """
    Annotate boxplots with significance stars, p-values, effect size, and CI.

    The function performs Mann-Whitney testing between two groups and
    computes a median-based effect size with 95% bootstrap CI.

    When `hue_columns` is provided, each x-category gets one test:
    `group1_name` vs `group2_name`. If `p_adjust_method` is set, p-values are
    adjusted across all valid x-categories (same family) using
    `statsmodels.stats.multitest.multipletests`.

    Parameters:
        ax (matplotlib.axes.Axes): Axis containing the boxplot.
        data (pd.DataFrame): Source dataframe.
        x_columns (str): Column mapped to x-axis categories.
        y_columns (str): Numeric response column.
        group1_name (str): First group label.
        group2_name (str): Second group label.
        y_pos (float): Y-position for annotation text.
        hue_columns (str, optional): Hue column for within-category comparisons.
        fontsize (int, optional): Annotation font size.
        p_adjust_method (str, optional): Method for p-value adjustment (for
            example `fdr_bh`, `bonferroni`). Uses raw p-values when None.
        alternative (str, optional): One of `greater`, `less`, or `two-sided`
            for Mann-Whitney testing. Default is `greater`.

    Returns:
        matplotlib.axes.Axes: Annotated axis.
    """
    xtick_labels = [tick.get_text() for tick in ax.get_xticklabels()]
    xtick_positions = [tick.get_position()[0] for tick in ax.get_xticklabels()]

    if hue_columns is not None:
        results = []
        for label, pos in zip(xtick_labels, xtick_positions):
            data_tmp = data[data[x_columns] == label]
            try:
                group1 = data_tmp[data_tmp[hue_columns] == group1_name][y_columns]
                group2 = data_tmp[data_tmp[hue_columns] == group2_name][y_columns]
                if len(group1) == 0 or len(group2) == 0:
                    print(f"Length of one of the groups is zero for {label}")
                    continue
            except Exception:
                print(f"There is some error while grouping wrt hue for label {label}")
                continue

            try:
                p_raw = man_whiteney(group1, group2, alternative=alternative)
                d_value, ci_low, ci_high = median_based_distance_ci(group1, group2)
            except Exception:
                print(f"There is some error while applying stats for label {label}")
                continue

            results.append(
                {
                    "label": label,
                    "pos": pos,
                    "p_raw": p_raw,
                    "p_reported": p_raw,
                    "d": d_value,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                }
            )

        if p_adjust_method and len(results) > 0:
            try:
                from statsmodels.stats.multitest import multipletests

                pvals = [res["p_raw"] for res in results]
                _, p_adj, _, _ = multipletests(pvals, method=p_adjust_method)
                for res, p_val in zip(results, p_adj):
                    res["p_reported"] = float(p_val)
            except Exception as exc:
                print(
                    f"Could not apply p-value correction ({p_adjust_method}); "
                    f"using raw p-values. Error: {exc}"
                )

        for res in results:
            star = get_significance_star(res["p_reported"])
            p_label = (
                f"p_raw: {res['p_raw']:.2e}\np_adjusted: {res['p_reported']:.2e}"
                if p_adjust_method
                else f"p_raw: {res['p_raw']:.2e}"
            )
            ax.text(
                res["pos"],
                y_pos,
                f"{star} \n {p_label} \nd: {res['d']:.2e} \n95% CI [{res['ci_low']:.2e}, {res['ci_high']:.2e}]",
                ha="center",
                fontsize=fontsize,
                color="red",
            )
            print(
                f"Label={res['label']}, p_raw={res['p_raw']:.2e}, "
                f"p_reported={res['p_reported']:.2e}, effect size={res['d']:.2e}, "
                f"95% CI [{res['ci_low']:.2e}, {res['ci_high']:.2e}]"
            )
        return ax

    if (hue_columns is None) and (len(xtick_labels) == 2):
        try:
            group1 = data[data[x_columns] == group1_name][y_columns]
            group2 = data[data[x_columns] == group2_name][y_columns]
            if len(group1) == 0 or len(group2) == 0:
                print("Length of one of the groups is zero")
                return ax
        except Exception:
            print("There is some error while grouping without hue")
            return ax

        try:
            p_value = man_whiteney(group1, group2, alternative=alternative)
            d_value, ci_low, ci_high = median_based_distance_ci(group1, group2)
            star = get_significance_star(p_value)
        except Exception:
            print("There is some error while applying stats without hue")
            return ax

        ax.text(
            0.5,
            y_pos,
            f"{star} \n p: {p_value:.2e} \nd: {d_value:.2e} \n95% CI [{ci_low:.2e}, {ci_high:.2e}]",
            ha="center",
            fontsize=fontsize,
            color="red",
        )
        print(
            f"P-value = {p_value:.2e}, effect size = {d_value:.2e}, "
            f"95% CI [{ci_low:.2e}, {ci_high:.2e}]"
        )
    return ax

    # except:

        # print("There is some error")


def summarize_boxplot_pairwise_stats(
    data,
    x_columns,
    y_columns,
    group1_name,
    group2_name,
    hue_columns=None,
    p_adjust_method=None,
    alternative="greater",
    ci=95,
):
    """
    Build a tidy statistics table for boxplot comparisons used in figures.

    Returns one row per comparison with:
    - category label
    - sample sizes
    - Mann-Whitney raw/adjusted p-values
    - median-based effect size and bootstrap CI
    """
    data = data.copy()
    results = []

    if hue_columns is not None:
        categories = pd.Series(data[x_columns].dropna().unique()).tolist()
        for label in categories:
            subset = data[data[x_columns] == label]
            group1 = subset[subset[hue_columns] == group1_name][y_columns].dropna()
            group2 = subset[subset[hue_columns] == group2_name][y_columns].dropna()
            if len(group1) == 0 or len(group2) == 0:
                continue

            p_raw = man_whiteney(group1, group2, alternative=alternative)
            effect_size, ci_low, ci_high = median_based_distance_ci(group1, group2, ci=ci)
            results.append(
                {
                    "category": label,
                    "x_column": x_columns,
                    "y_column": y_columns,
                    "hue_column": hue_columns,
                    "group1_name": group1_name,
                    "group2_name": group2_name,
                    "n_group1": int(group1.shape[0]),
                    "n_group2": int(group2.shape[0]),
                    "test": "Mann-Whitney U",
                    "alternative": alternative,
                    "effect_size_name": "median_based_distance",
                    "effect_size": effect_size,
                    "ci_level": ci,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "p_raw": float(p_raw),
                    "p_adjust_method": p_adjust_method if p_adjust_method else "none",
                    "p_adjusted": float(p_raw),
                }
            )
    else:
        group1 = data[data[x_columns] == group1_name][y_columns].dropna()
        group2 = data[data[x_columns] == group2_name][y_columns].dropna()
        if len(group1) > 0 and len(group2) > 0:
            p_raw = man_whiteney(group1, group2, alternative=alternative)
            effect_size, ci_low, ci_high = median_based_distance_ci(group1, group2, ci=ci)
            results.append(
                {
                    "category": "overall",
                    "x_column": x_columns,
                    "y_column": y_columns,
                    "hue_column": hue_columns,
                    "group1_name": group1_name,
                    "group2_name": group2_name,
                    "n_group1": int(group1.shape[0]),
                    "n_group2": int(group2.shape[0]),
                    "test": "Mann-Whitney U",
                    "alternative": alternative,
                    "effect_size_name": "median_based_distance",
                    "effect_size": effect_size,
                    "ci_level": ci,
                    "ci_low": ci_low,
                    "ci_high": ci_high,
                    "p_raw": float(p_raw),
                    "p_adjust_method": p_adjust_method if p_adjust_method else "none",
                    "p_adjusted": float(p_raw),
                }
            )

    results_df = pd.DataFrame(results)
    if results_df.empty:
        return results_df

    if p_adjust_method and len(results_df) > 0:
        try:
            from statsmodels.stats.multitest import multipletests

            _, p_adj, _, _ = multipletests(results_df["p_raw"].tolist(), method=p_adjust_method)
            results_df["p_adjusted"] = p_adj.astype(float)
        except Exception as exc:
            print(
                f"Could not apply p-value correction ({p_adjust_method}); "
                f"using raw p-values. Error: {exc}"
            )

    return results_df
