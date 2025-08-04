"""
Core calculation functions for batch effect and biological effect analysis.

This module implements the mathematical calculations for both batch effect analysis, including:
- Likelihood estimation formula for each parameter
- Calculation of q_sh and q_sp scores
"""

import logging
from typing import Dict, Tuple
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData
from .models import bioLUCIDResult

logger = logging.getLogger(__name__)

def calculate_per_sample_cluster_means(data: pd.DataFrame) -> pd.Series:
    """
    Calculate group-specific means.
    
    Args:
        data: DataFrame with columns [sample_id, cluster_id, gene_id, expression_value]
        
    Returns:
        X_st: Mean expression per sample-cluster-gene combination
    """
    X_st = data.groupby(['sample_id', 'cluster_id', 'gene_id'])['expression_value'].mean()
    
    return X_st

def calculate_μ(data: pd.Series) -> pd.Series:
    """
    Calculate μ (global mean expression per gene).
    Likelihood estimation formula: μ = 1/st * ΣstX_st

    Args:
        data: X_st
        
    Returns:
        μ: Global mean expression per gene
    """
    
    μ = data.groupby('gene_id').mean()

    return μ

def calculate_cluster_effect(
    X_st: pd.Series,
    μ: pd.Series
) -> pd.Series:
    """
    Calculate cluster effect (tau_t).
    Likelihood estimation formula: tau_t = 1/s * Σs(X_st - μ)
    
    Args:
        X_st: Sample-cluster means
        μ: Global means
        
    Returns:
        Cluster effect per cluster-gene combination
    """
    # Prepare data
    X_st_df = X_st.reset_index()
    
    # Merge and calculate differences
    merged = X_st_df.merge(
        μ.reset_index(),
        on=['gene_id'],
        suffixes=('_st', '_mean')
    )
    
    # Calculate cluster effect
    tau_t = (merged
              .groupby(['cluster_id', 'gene_id'])
              .apply(lambda x: (x['expression_value_st'] - x['expression_value_mean']).mean()))
    
    return tau_t

def calculate_sample_effect(
    X_st: pd.Series,
    μ: pd.Series,
) -> pd.Series:
    """
    Calculate sample effect (sigma_s).
    Likelihood estimation formula: sigma_s = 1/t * Σt(X_st - μ)
    
    Args:
        X_st: Sample-cluster means
        μ: Global means
        
    Returns:
        Sample effect per sample-gene combination
    """
    # Prepare data
    X_st_df = X_st.reset_index()
    merged = X_st_df.merge(
        μ.reset_index(),
        on=['gene_id'],
        suffixes=('_st', '_mean')
    )
    
    # Calculate sample effect
    sigma_s = (merged
              .groupby(['sample_id', 'gene_id'])
              .apply(lambda x: (x['expression_value_st'] - 
                              x['expression_value_mean']).mean()))
    
    return sigma_s

def calculate_residuals(
    X_st: pd.Series,
    μ: pd.Series,
    tau_t: pd.Series,
    sigma_s: pd.Series
) -> pd.DataFrame:
    """
    Calculate residuals (ϵst).
    Formula: ϵst = X_st - μ - sigma_s - tau_t
    
    Args:
        X_st: Sample-cluster means
        μ: Global means
        tau_t: Cluster effects
        sigma_s: Sample effects
        
    Returns:
        Residuals DataFrame
    """
    # Prepare data
    X_st_df = X_st.reset_index()
    μ_df = μ.reset_index()
    tau_t_df = tau_t.reset_index().rename(columns={0: 'tau_t_value'})
    sigma_s_df = sigma_s.reset_index().rename(columns={0: 'sigma_s_value'})
    
    # Merge - the first merge will create _x and _y suffixes for expression_value
    merged = (X_st_df
             .merge(μ_df, on=['gene_id'])  # Creates expression_value_x (X_st) and expression_value_y (μ)
             .merge(tau_t_df, on=['cluster_id', 'gene_id'])
             .merge(sigma_s_df, on=['sample_id', 'gene_id']))
    
    # Calculate residuals using the pandas-generated column names
    residuals = (merged['expression_value_x'] -    # X_st values
                merged['expression_value_y'] -     # μ values
                merged['tau_t_value'] -             # τ_t values
                merged['sigma_s_value'])            # σ_s values
    
    residuals_df = pd.DataFrame({
        'sample_id': merged['sample_id'],
        'cluster_id': merged['cluster_id'],
        'gene_id': merged['gene_id'],
        'residuals': residuals
    })

    return residuals_df

def calculate_q_sh_per_sample(sigma_s: pd.Series) -> Dict[str, float]:
    """
    Calculate amount of gene expression variation that is shared across celltypes per sample.
    Formula: q_sh = (1/g * Σg(sigma_s)²)^(1/2)
    
    Args:
        sigma_s: batch effects vector (sample-gene level)
        
    Returns:
        Dict mapping sample_id to q_sh score
    """
    sigma_s_df = sigma_s.reset_index()
    
    # Calculate q_sh per sample
    q_sh_per_sample = {
        sample: np.sqrt(np.mean(group[0] ** 2))
        for sample, group in sigma_s_df.groupby('sample_id')
    }
    
    return q_sh_per_sample

def calculate_q_sp_per_sample(residuals_df: pd.DataFrame) -> Dict[str, float]:
    """
    Calculate amount of gene expression variation that is specific to a cell type per sample.
    Formula: q_sp = (1/gt * Σgt(ϵst)²)^(1/2)
    
    Args:
        residuals_df: DataFrame containing residuals
        
    Returns:
        Dict mapping sample_id to q_sp score
    """
    # Calculate q_sp per sample (average over genes and cell types within each sample)
    q_sp_per_sample = {
        sample: np.sqrt(np.mean(group['residuals'] ** 2))
        for sample, group in residuals_df.groupby('sample_id')
    }
    
    return q_sp_per_sample

def calculate_b_scores_per_sample(
    q_sh_per_sample: Dict[str, float],
    q_sp_per_sample: Dict[str, float]
) -> Dict[str, float]:
    """
    Calculate proportion of batch effect per sample (b score).
    Formula: b_s = q_sh² / (q_sh² + q_sp²)
    
    Args:
        q_sh_per_sample: Dict of q_sh scores per sample
        q_sp_per_sample: Dict of q_sp scores per sample
        
    Returns:
        Dict of b scores per sample
    """
    b_scores_per_sample = {}
    
    for sample in q_sh_per_sample.keys():
        q_sh_squared = q_sh_per_sample[sample] ** 2
        q_sp_squared = q_sp_per_sample[sample] ** 2
        b_scores_per_sample[sample] = q_sh_squared / (q_sh_squared + q_sp_squared)
    
    return b_scores_per_sample

def calculate_global_scores(
    q_sh_per_sample: Dict[str, float],
    q_sp_per_sample: Dict[str, float],
    b_scores_per_sample: Dict[str, float]
) -> Tuple[float, float, float]:
    """
    Calculate global scores by averaging across samples.
    
    Returns:
        Tuple of (global_q_sh, global_q_sp, global_b_score)
    """
    global_q_sh = np.mean(list(q_sh_per_sample.values()))
    global_q_sp = np.mean(list(q_sp_per_sample.values()))
    global_b_score = np.mean(list(b_scores_per_sample.values()))
    
    return global_q_sh, global_q_sp, global_b_score

def bioLUCID_calculation(
    adata: AnnData,
    params: Dict
) -> bioLUCIDResult:
    """
    Calculate batch and biological effect scores and store them in a bioLUCIDResult object.
    
    Args:
        adata: Input data
        params: Analysis parameters
        
    Returns:
        bioLUCIDResult containing scores and components
    """
    
    # Prepare data
    if 'logTPM' not in adata.layers:
        raise ValueError("logTPM layer not found. Run preprocessing first.")
    
    # Create DataFrame for calculations
    cells_df = pd.DataFrame(
        adata.layers['logTPM'],
        columns=adata.var_names
    )
    cells_df['sample_id'] = adata.obs[params['batch_key']].values
    cells_df['cluster_id'] = adata.obs[params['celltype_key']].values
    
    # Melt data
    data = cells_df.melt(
        id_vars=['sample_id', 'cluster_id'],
        var_name='gene_id',
        value_name='expression_value'
    )
    
    # Calculate components according to mentor's formulas
    X_st = calculate_per_sample_cluster_means(data)
    μ = calculate_μ(X_st)
    tau_t = calculate_cluster_effect(X_st, μ)
    sigma_s = calculate_sample_effect(X_st, μ)
    
    # Calculate residuals
    residuals_df = calculate_residuals(X_st, μ, tau_t, sigma_s)

    # Calculate batch effect scores per sample
    q_sh_per_sample = calculate_q_sh_per_sample(sigma_s)
    q_sp_per_sample = calculate_q_sp_per_sample(residuals_df)
    b_scores_per_sample = calculate_b_scores_per_sample(q_sh_per_sample, q_sp_per_sample)
    
    # Calculate global scores
    global_q_sh, global_q_sp, global_b_score = calculate_global_scores(
        q_sh_per_sample, q_sp_per_sample, b_scores_per_sample
    )
    
    return bioLUCIDResult(

        # Global scores
        Global_b_score=global_b_score,
        Global_q_sh_score=global_q_sh,
        Global_q_sp_score=global_q_sp,

        # Per-sample scores
        b_score_per_batch=b_scores_per_sample,
        q_sh_score_per_batch=q_sh_per_sample,
        q_sp_score_per_batch=q_sp_per_sample,

        # Components
        components={
            'μ': μ,
            'X_st': X_st,
            'tau_t': tau_t,
            'sigma_s': sigma_s,
            'residuals':residuals_df
        }
    )