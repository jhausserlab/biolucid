"""
Core calculation functions for batch effect and biological effect analysis.

This module implements the mathematical calculations for both batch effect
and biological effect analysis, including:
- Batch effect component calculation
- Differential expression analysis
- Pathway enrichment analysis
- Score calculations
"""

import logging
from typing import Dict, Tuple, List, Optional
import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc
import gseapy as gp
from anndata import AnnData
from .models import BatchEffectResult, BiologicalEffectResult
from ..config.settings import get_pathway_dbs

logger = logging.getLogger(__name__)

def calculate_means(data: pd.DataFrame) -> Tuple[pd.Series, pd.Series]:
    """
    Calculate global and group-specific means.
    
    Args:
        data: DataFrame with columns [sample_id, cluster_id, gene_id, expression_value]
        
    Returns:
        Tuple containing:
        - X_mean: Global mean expression per gene
        - X_st: Mean expression per sample-cluster-gene combination
    """
    # Calculate global mean per gene
    X_mean = data.groupby('gene_id')['expression_value'].mean()
    
    # Calculate mean per sample-cluster combination
    X_st = data.groupby(['sample_id', 'cluster_id', 'gene_id'])['expression_value'].mean()
    
    return X_mean, X_st

def calculate_cluster_effect(
    X_st: pd.Series,
    X_mean: pd.Series
) -> pd.Series:
    """
    Calculate cluster effect (∆t).
    Formula: ∆t = 1/s * Σs(X_st - X_mean)
    
    Args:
        X_st: Sample-cluster means
        X_mean: Global means
        
    Returns:
        Cluster effect per cluster-gene combination
    """
    # Prepare data
    X_st_df = X_st.reset_index()
    
    # Merge and calculate differences
    merged = X_st_df.merge(
        X_mean.reset_index(),
        on=['gene_id'],
        suffixes=('_st', '_mean')
    )
    
    # Calculate cluster effect
    delta_t = (merged
              .groupby(['cluster_id', 'gene_id'])
              .apply(lambda x: (x['expression_value_st'] - x['expression_value_mean']).mean()))
    
    return delta_t

def calculate_sample_effect(
    X_st: pd.Series,
    X_mean: pd.Series,
    delta_t: pd.Series
) -> pd.Series:
    """
    Calculate sample effect (∆s).
    Formula: ∆s = 1/t * Σt(X_st - X_mean - ∆t)
    
    Args:
        X_st: Sample-cluster means
        X_mean: Global means
        delta_t: Cluster effects
        
    Returns:
        Sample effect per sample-gene combination
    """
    # Prepare data
    X_st_df = X_st.reset_index()
    merged = (X_st_df
             .merge(X_mean.reset_index(),
                   on=['gene_id'],
                   suffixes=('_st', '_mean'))
             .merge(delta_t.reset_index(),
                   on=['cluster_id', 'gene_id']))
    
    # Calculate sample effect
    delta_s = (merged
              .groupby(['sample_id', 'gene_id'])
              .apply(lambda x: (x['expression_value_st'] - 
                              x['expression_value_mean'] - 
                              x[0]).mean()))
    
    return delta_s

def calculate_batch_effect(
    adata: AnnData,
    params: Dict
) -> BatchEffectResult:
    """
    Calculate batch effect scores.
    
    Args:
        adata: Input data
        params: Analysis parameters
        
    Returns:
        BatchEffectResult containing scores and components
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
    
    # Calculate components
    X_mean, X_st = calculate_means(data)
    delta_t = calculate_cluster_effect(X_st, X_mean)
    delta_s = calculate_sample_effect(X_st, X_mean, delta_t)
    
    # Calculate residuals
    residuals = calculate_residuals(X_st, X_mean, delta_t, delta_s)
    
    # Calculate scores
    B_scores = calculate_B_scores(delta_s)
    b_scores = calculate_b_scores(B_scores, residuals)
    
    return BatchEffectResult(
        b_score=b_scores['global'],
        B_score=B_scores['global'],
        b_scores_per_batch=b_scores['per_batch'],
        B_scores_per_batch=B_scores['per_batch'],
        components={
            'X_mean': X_mean,
            'X_st': X_st,
            'delta_t': delta_t,
            'delta_s': delta_s
        },
        residuals=residuals
    )

def calculate_residuals(
    X_st: pd.Series,
    X_mean: pd.Series,
    delta_t: pd.Series,
    delta_s: pd.Series
) -> pd.Series:
    """
    Calculate residuals (ϵst).
    Formula: ϵst = X_st - X_mean - ∆s - ∆t
    
    Args:
        X_st: Sample-cluster means
        X_mean: Global means
        delta_t: Cluster effects
        delta_s: Sample effects
        
    Returns:
        Residual values
    """
    # Prepare data
    X_st_df = X_st.reset_index()
    
    # Merge all components
    merged = (X_st_df
             .merge(X_mean.reset_index(),
                   on=['gene_id'],
                   suffixes=('_st', '_mean'))
             .merge(delta_t.reset_index(),
                   on=['cluster_id', 'gene_id'])
             .merge(delta_s.reset_index(),
                   on=['sample_id', 'gene_id']))
    
    # Calculate residuals
    residuals = (merged['expression_value_st'] - 
                merged['expression_value_mean'] - 
                merged['0_x'] -  # delta_t
                merged['0_y'])   # delta_s
    
    residuals_df = pd.DataFrame({
        'sample_id': X_st_df['sample_id'],
        'cluster_id': X_st_df['cluster_id'],
        'gene_id': X_st_df['gene_id'],
        'residuals': residuals
    })

    return residuals_df

def calculate_B_scores(delta_s: pd.Series) -> Dict[str, float]:
    """
    Calculate batch effect magnitude (B score).
    Formula: B = 1/sg * Σsg(∆sg)²
    
    Args:
        delta_s: Sample effects
        
    Returns:
        Dictionary with global and per-batch B scores
    """
    delta_s_df = delta_s.reset_index()
    
    # Calculate global B score
    B_score = np.mean(delta_s_df[0] ** 2)
    
    # Calculate per-batch B scores
    B_scores_per_batch = {
        sample: np.mean(group[0] ** 2)
        for sample, group in delta_s_df.groupby('sample_id')
    }
    
    return {
        'global': B_score,
        'per_batch': B_scores_per_batch
    }

def calculate_b_scores(
    B_scores: Dict[str, float],
    residuals_df: pd.DataFrame
) -> Dict[str, float]:
    """
    Calculate proportion of batch effect (b score).
    Formula: b = B / (B + 1/sgt·Σϵstg²)
    
    Args:
        B_scores: B score values
        residuals_df: DataFrame containing residuals and batch information
        
    Returns:
        Dictionary with global and per-batch b scores
    """
    # Calculate global epsilon sum
    epsilon_sum = np.mean(residuals_df['residuals'] ** 2)
    
    epsilon_sums_per_batch = {
        sample: np.mean(group['residuals'] ** 2)
        for sample, group in residuals_df.groupby('sample_id')
    }

    # Calculate global b score
    b_score = B_scores['global'] / (B_scores['global'] + epsilon_sum)
    
    # Calculate per-batch b scores
    b_scores_per_batch = {
        sample: B_s / (B_s + eps_sum)
        for sample, (B_s, eps_sum) in 
        zip(B_scores['per_batch'].keys(),
            zip(B_scores['per_batch'].values(),
                epsilon_sums_per_batch.values()))
    }
    
    return {
        'global': b_score,
        'per_batch': b_scores_per_batch
    }

def calculate_expression_profile(
    adata: AnnData,
    batch: str,
    celltype: str,
    params: Dict
) -> pd.Series:
    """
    Calculate average gene expression profile for a specific batch-celltype combination.
    
    Args:
        adata: Input data
        batch: Batch identifier
        celltype: Cell type identifier
        params: Analysis parameters
        
    Returns:
        Average expression profile as pandas Series
    """
    # Get specific batch and celltype cells
    mask = (adata.obs[params['batch_key']] == batch) & \
           (adata.obs[params['celltype_key']] == celltype)
    subset = adata[mask]

    if 'logTPM' not in subset.layers:
        raise ValueError("logTPM layer not found. Run preprocessing first.")
    
    # Calculate mean expression
    if scipy.sparse.issparse(subset.layers['logTPM']):
        profile = pd.Series(
            np.array(subset.layers['logTPM'].mean(axis=0)).flatten(),
            index=adata.var_names
        )
    else:
        profile = pd.Series(
            subset.layers['logTPM'].mean(axis=0),
            index=adata.var_names
        )
    
    return profile

def get_ranked_genes(
    adata: AnnData,
    batch: str,
    celltype: str,
    params: Dict
) -> pd.Series:
    """
    Get ranked gene list based on expression profile difference.
    
    Args:
        adata: Input data
        batch: Batch identifier
        celltype: Cell type identifier
        params: Analysis parameters
        
    Returns:
        Ranked gene list with differences as values
    """
    # Get batch-specific profile
    batch_profile = calculate_expression_profile(adata, batch, celltype, params)
    
    # Get overall celltype profile (including the batch)
    celltype_mask = adata.obs[params['celltype_key']] == celltype
    celltype_adata = adata[celltype_mask]
    
    if scipy.sparse.issparse(celltype_adata.layers['logTPM']):
        overall_profile = pd.Series(
            np.array(celltype_adata.layers['logTPM'].mean(axis=0)).flatten(),
            index=adata.var_names
        )
    else:
        overall_profile = pd.Series(
            celltype_adata.layers['logTPM'].mean(axis=0),
            index=adata.var_names
        )
    
    # Calculate difference
    diff_profile = batch_profile - overall_profile
    
    # Sort by from highest to lowest difference
    ranked_genes = diff_profile.sort_values(ascending=False)
    
    return ranked_genes

def run_enrichment(
    ranked_genes: pd.Series,
    species: str,
    params: Dict
) -> Optional[pd.DataFrame]:
    """
    Run GSEA enrichment analysis.
    
    Args:
        ranked_genes: Ranked gene list
        species: Species name
        params: Analysis parameters
        
    Returns:
        Enrichment results or None if analysis fails
    """
    if ranked_genes.empty:
        return None
        
    try:
        # Get appropriate databases
        databases = get_pathway_dbs(species)
        
        # Convert series to rnk format
        rnk = pd.Series(
            ranked_genes.values,
            index=ranked_genes.index
        )
        
        # Run GSEA
        enr = gp.prerank(
            rnk=rnk,
            gene_sets=databases,
            organism=species,
            processes=params['processes'],
            permutation_num=params['permutation_num'],
            min_size=params['min_size'],
            max_size=params['max_size'],
            seed=params['seed']
        )
        
        return enr.res2d
        
    except Exception as e:
        logger.warning(f"Enrichment analysis failed: {str(e)}")
        return None

def calculate_pathway_score(
    enrichment_results: pd.DataFrame
) -> float:
    """
    Calculate pathway-based biological effect score using FDR q-values.
    """
    if enrichment_results is None or len(enrichment_results) == 0:
        return 0.0
    
    # Use FDR q-val instead of adjusted p-value
    score = sum(1 - enrichment_results['FDR q-val']) / len(enrichment_results)
    return score

def calculate_biological_effect(adata: AnnData, params: Dict) -> BiologicalEffectResult:
    """Calculate biological effect scores."""
    # Verify required data layers
    for layer in ['counts', 'logTPM']:
        if layer not in adata.layers:
            raise ValueError(f"{layer} layer not found. Run preprocessing first.")
    
    batches = adata.obs[params['batch_key']].unique()
    celltypes = adata.obs[params['celltype_key']].unique()
    
    logger.info(f"Analyzing {len(batches)} batches and {len(celltypes)} cell types")
    
    celltype_scores = {}
    batch_scores = {}
    enrichment_results = {}
    
    for batch in batches:
        celltype_scores[batch] = {}
        enrichment_results[batch] = {}
        
        for celltype in celltypes:
            logger.info(f"Processing batch {batch}, celltype {celltype}")
            
            # Get ranked genes
            ranked_genes = get_ranked_genes(adata, batch, celltype, params)
            logger.info(f"Found and Rank {len(ranked_genes)} differentially expressed genes")
            
            # Run enrichment and calculate score
            if len(ranked_genes) > 0:
                enrichment = run_enrichment(ranked_genes, params['species'], params)
                if enrichment is not None and len(enrichment) > 0:
                    score = calculate_pathway_score(enrichment)
                    logger.info(f"Enrichment analysis successful, score: {score:.3f}")
                else:
                    score = 0.0
                    logger.warning("No significant enrichment found")
            else:
                score = 0.0
                logger.warning("No differentially expressed genes found")
            
            celltype_scores[batch][celltype] = score
            enrichment_results[batch][celltype] = enrichment
        
        # Calculate batch score
        batch_scores[batch] = np.mean(list(celltype_scores[batch].values()))
        logger.info(f"Batch {batch} score: {batch_scores[batch]:.3f}")
    
    # Calculate summary score
    summary_score = np.mean(list(batch_scores.values()))
    logger.info(f"Overall biological effect score: {summary_score:.3f}")
    
    return BiologicalEffectResult(
        Biological_summary_score=summary_score,
        Biological_batch_scores=batch_scores,
        Biological_celltype_scores=celltype_scores,
        Biological_enrichment_results=enrichment_results
    ) 