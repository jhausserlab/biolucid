"""
Visualization functions for batch effect analysis results.

This module provides standalone visualization functions that can be used
to plot various aspects of the batch effect analysis results.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from typing import Dict, Optional
from scipy import stats

def create_analysis_df(results):
    """
    Create a DataFrame containing all analysis metrics from results dictionary.
    
    Args:
        results: Dictionary containing analysis results with required keys:
                - Biological_batch_scores
                - B_scores_per_batch
                - b_scores_per_batch
                - residuals: DataFrame with 'residuals' and 'sample_id' columns
    
    Returns:
        pandas.DataFrame: DataFrame containing all analysis metrics
    """
    # Calculate residuals sums
    residuals_sums_per_batch = {
        sample: np.mean(group['residuals'] ** 2)
        for sample, group in results['residuals'].groupby('sample_id')
    }
    
    # Create base DataFrame
    results_df = pd.DataFrame({
        'Biological Effect Scores': results['Biological_batch_scores'],
        'Batch Effect Scores': results['B_scores_per_batch'],
        'Normalized Batch Effect Scores': results['b_scores_per_batch'],
        'Residuals Sums': residuals_sums_per_batch
    })
    
    # Add square root columns
    results_df['Batch Effect Scores_sqrt'] = np.sqrt(results_df['Batch Effect Scores'])
    results_df['Normalized Batch Effect Scores_sqrt'] = np.sqrt(results_df['Normalized Batch Effect Scores'])
    results_df['Residuals Sums_sqrt'] = np.sqrt(results_df['Residuals Sums'])
    
    return results_df

def plot_scatter_analysis(results_df, figsize=(15, 5)):
    """
    Plot three scatter plots for analysis visualization.
    
    Args:
        results_df: DataFrame created by create_analysis_df function
        figsize: Figure size tuple (width, height)
    """
    plt.figure(figsize=figsize)
    
    # Plot 1: Normalized batch effect vs Biological Effect Scores
    plt.subplot(1, 3, 1)
    sns.scatterplot(data=results_df, x='Biological Effect Scores', 
                   y='Normalized Batch Effect Scores')
    plt.title('Normalized Batch Effect vs\nBiological Effect Scores')
    for idx, row in results_df.iterrows():
        plt.annotate(idx, (row['Biological Effect Scores'], 
                         row['Normalized Batch Effect Scores']))
    
    # Plot 2: Batch Effect Scores_sqrt vs Residuals Sums_sqrt
    plt.subplot(1, 3, 2)
    max_val = max(results_df['Batch Effect Scores_sqrt'].max(), 
                 results_df['Residuals Sums_sqrt'].max())
    sns.scatterplot(data=results_df, x='Residuals Sums_sqrt', 
                   y='Batch Effect Scores_sqrt')
    plt.plot([0, max_val], [0, max_val], 'k:', label='y=x')
    plt.title('Batch Effect Scores (sqrt) vs\nResiduals Sums (sqrt)')
    plt.legend()
    for idx, row in results_df.iterrows():
        plt.annotate(idx, (row['Residuals Sums_sqrt'], 
                         row['Batch Effect Scores_sqrt']))
    
    # Plot 3: Biological Effect Scores vs Residuals Sums_sqrt with correlation
    plt.subplot(1, 3, 3)
    sns.scatterplot(data=results_df, x='Residuals Sums_sqrt', 
                   y='Biological Effect Scores')
    correlation, p_value = stats.pearsonr(results_df['Residuals Sums_sqrt'], 
                                        results_df['Biological Effect Scores'])
    plt.title(f'Biological Effect Scores vs\nResiduals Sums (sqrt)\nr={correlation:.2f}, p={p_value:.3f}')
    for idx, row in results_df.iterrows():
        plt.annotate(idx, (row['Residuals Sums_sqrt'], 
                         row['Biological Effect Scores']))
    
    plt.tight_layout()
    plt.show()

def plot_scores_bar(results_df, score_type='Batch Effect Scores', sqrt=False, 
                   color='skyblue', figsize=(10, 6), rotation=45, **kwargs):
    """
    Plot bar chart for any score type from the analysis DataFrame.
    
    Args:
        results_df: DataFrame created by create_analysis_df function
        score_type: Type of score to plot, options:
                   - 'Batch Effect Scores'
                   - 'Normalized Batch Effect Scores'
                   - 'Biological Effect Scores'
        sqrt: Whether to use square root values
        color: Bar color
        figsize: Figure size tuple (width, height)
        rotation: X-axis label rotation angle
        **kwargs: Additional arguments passed to plt.bar()
    """
    plt.figure(figsize=figsize)
    
    # Select data to plot
    column = f'{score_type}_sqrt' if sqrt else score_type
    if column not in results_df.columns:
        raise ValueError(f"Column {column} not found in DataFrame")
    
    # Create bar plot
    bars = plt.bar(results_df.index, results_df[column], color=color, **kwargs)
    
    # Add value labels
    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height,
                f'{height:.2e}',
                ha='center', va='bottom')
    
    # Customize plot
    plt.xlabel('Batch')
    ylabel = f'sqrt({score_type})' if sqrt else score_type
    plt.ylabel(ylabel)
    plt.title(f'{ylabel} per Batch')
    plt.xticks(rotation=rotation)
    
    plt.tight_layout()
    plt.show() 