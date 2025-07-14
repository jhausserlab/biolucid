"""
Visualization functions for batch effect analysis results.

This module provides standalone visualization functions that can be used
to plot various aspects of the batch effect analysis results.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy import stats
from typing import Tuple

def results_to_df(results) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Create a DataFrame containing all scores from results dictionary.
    
    Args:
        results: Dictionary containing analysis results with required keys:
                - Global_q_sh_score
                - q_sh_score_per_batch
                - b_score_per_batch
                - residuals: DataFrame with 'residuals' and 'sample_id' columns
    
    Returns:
        Tuple[pd.DataFrame, pd.DataFrame]: Tuple containing global and per-sample results DataFrames
    """
    
    # Create base DataFrame
    global_results_df = pd.DataFrame({
    'Global_b_score': [results.Global_b_score],
    'Global_q_sp_score': [results.Global_q_sp_score], 
    'Global_q_sh_score': [results.Global_q_sh_score],}, index=[0]) 

    per_sample_results_df = pd.DataFrame({
        'q_sh_score_per_batch': results.q_sh_score_per_batch,
        'q_sp_score_per_batch': results.q_sp_score_per_batch,
        'b_score_per_batch': results.b_score_per_batch,
    })

    # add selection suggestions
    per_sample_results_df['recommendation'] = pd.cut(
        per_sample_results_df['b_score_per_batch'],
        bins=[-float('inf'), 0.6, 0.8, float('inf')],
        labels=['Keep', 'Keep with caution', 'Drop'],
        include_lowest=True
    )
    
    return global_results_df, per_sample_results_df

def plot_scatter_analysis(results_df, figsize=(5, 5)):
    """
    Plot scatter plots for analysis visualization based on per-sample results.
    
    Args:
        results_df: DataFrame created by create_analysis_df function
        figsize: Figure size tuple (width, height)
    """
    plt.figure(figsize=figsize)
    
    # Plot 1: q_sh_score_per_batch vs q_sp_score_per_batch
    max_val = max(results_df['q_sh_score_per_batch'].max(), 
                 results_df['q_sp_score_per_batch'].max())
    sns.scatterplot(data=results_df, x='q_sp_score_per_batch', 
                   y='q_sh_score_per_batch')
    plt.plot([0, max_val], [0, max_val], 'k:', label='y=x')
    plt.title('q_sp_vs_q_sh_score_per_batch')
    plt.legend()
    for idx, row in results_df.iterrows():
        plt.annotate(idx, (row['q_sp_score_per_batch'], 
                         row['q_sh_score_per_batch']))
    
    plt.tight_layout()
    plt.show()