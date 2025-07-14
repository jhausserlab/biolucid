"""
Utility functions for data preprocessing and analysis.

This module contains functions for data preprocessing, validation,
and general utility functions used throughout the analysis pipeline.
"""

import logging
from typing import Tuple, Dict
import numpy as np
import pandas as pd
import scipy.sparse
import scanpy as sc
from anndata import AnnData

# Set up logging
logger = logging.getLogger(__name__)

def setup_logging(verbose: bool = True):
    """
    Configure logging settings.
    
    Args:
        verbose: Whether to show INFO level logs
    """
    level = logging.INFO if verbose else logging.WARNING
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )

def validate_and_filter_celltypes(adata: AnnData, params: Dict) -> AnnData:
    """
    Validate batch and cell type information and filter cell types.
    Ensures each cell type has sufficient cells in EACH batch.
    
    Args:
        adata: Input data
        params: Analysis parameters
        
    Returns:
        Filtered AnnData with valid cell types
        
    Raises:
        ValueError: If validation fails
    """
    # Check batch information
    n_batches = adata.obs[params['batch_key']].nunique()
    if n_batches < 2:
        raise ValueError(f"Found only {n_batches} batch(es). At least 2 batches are required for analysis")
    
    # Calculate cell counts for each batch-celltype combination
    batch_celltype_counts = pd.crosstab(
        adata.obs[params['batch_key']], 
        adata.obs[params['celltype_key']]
    )
    
    # Find cell types that have sufficient cells in ALL batches
    min_cells = params.get('min_cells', 20)
    valid_celltypes = []
    
    for celltype in batch_celltype_counts.columns:
        counts = batch_celltype_counts[celltype]
        if all(counts >= min_cells):
            valid_celltypes.append(celltype)
    
    if len(valid_celltypes) < 2:
        # Prepare detailed error message
        insufficient_counts = batch_celltype_counts.loc[:, ~batch_celltype_counts.columns.isin(valid_celltypes)]
        error_msg = (
            f"Found only {len(valid_celltypes)} cell type(s) with >= {min_cells} cells in all batches. "
            f"At least 2 valid cell types are required.\n"
            f"Cell counts per batch:\n{insufficient_counts}"
        )
        raise ValueError(error_msg)
    
    # Filter data to keep only valid cell types
    mask = adata.obs[params['celltype_key']].isin(valid_celltypes)
    adata_ct_selection = adata[mask].copy()
    
    # Log detailed information
    logger.info(f"Retained {len(adata_ct_selection)} / {len(adata)} cells after cell type filtering")
    logger.info(f"Retained cell types: {valid_celltypes}")
    logger.info("Cell counts per batch for retained cell types:")
    retained_counts = batch_celltype_counts.loc[:, valid_celltypes]
    for batch in retained_counts.index:
        logger.info(f"Batch {batch}:")
        for celltype in valid_celltypes:
            count = retained_counts.loc[batch, celltype]
            logger.info(f"  {celltype}: {count} cells")
    
    return adata_ct_selection

def select_abundant_genes(adata: AnnData, params: Dict) -> Tuple[AnnData, np.ndarray]:
    """
    Select genes with sufficient expression across all cells.
    
    Args:
        adata: Input data after cell type filtering
        params: Analysis parameters
        
    Returns:
        Tuple of (filtered AnnData, abundant gene indices)
        
    Raises:
        ValueError: If insufficient abundant genes found
    """
    # Get raw counts
    counts = get_counts_matrix(adata)
    
    # Calculate mean UMI per cell for each gene
    gene_means = np.array(counts.mean(axis=0)).flatten()
    abundant_genes = np.where(gene_means >= params.get('abundant_gene_threshold', 1))[0]
    
    min_genes = params.get('min_abundant_genes', 300)
    if len(abundant_genes) < min_genes:
        raise ValueError(
            f"Found only {len(abundant_genes)} abundant genes (>= {params['abundant_gene_threshold']} UMI/cell). "
            f"At least {min_genes} genes are required for meaningful analysis"
        )
    
    # Create filtered AnnData
    adata_ct_genes_selection = adata[:, abundant_genes].copy()
    
    logger.info(f"Selected {len(abundant_genes)} abundant genes")
    return adata_ct_genes_selection, abundant_genes

def get_counts_matrix(adata: AnnData) -> np.ndarray:
    """
    Locate and return raw counts matrix.
    
    Args:
        adata: Input data
        
    Returns:
        Raw counts matrix
        
    Raises:
        ValueError: If no counts matrix is found
    """
    def is_counts(X) -> bool:
        """Check if matrix contains count data."""
        if scipy.sparse.issparse(X):
            return (np.issubdtype(X.data.dtype, np.integer) or 
                    np.all(np.mod(X.data, 1) == 0))
        else:
            return (np.issubdtype(X.dtype, np.integer) or 
                    np.all(np.mod(X, 1) == 0))
    
    # Check possible locations
    locations = {
        'layers["counts"]': (
            'counts' in adata.layers,
            lambda: adata.layers['counts']
        ),
        'raw.X': (
            adata.raw is not None and is_counts(adata.raw.X),
            lambda: adata.raw.X
        ),
        'X': (
            is_counts(adata.X),
            lambda: adata.X
        )
    }
    
    for loc_name, (exists, getter) in locations.items():
        if exists:
            logger.info(f"Found counts in {loc_name}")
            counts = getter()
            return scipy.sparse.csr_matrix(counts).toarray()
            
    raise ValueError("No raw counts found in data")

def normalize_data(adata: AnnData) -> AnnData:
    """
    Normalize data using median counts and log transform.
    
    Args:
        adata: Input data
        
    Returns:
        Normalized data with 'logTPM' layer
    """
    adata = adata.copy()
    
    # Get counts
    adata.X = get_counts_matrix(adata)

    # Store in layers['counts']
    logger.info("Storing raw counts in 'layers[\"counts\"]'")
    adata.layers['counts'] = adata.X.copy()

    logger.debug(f"Data range after getting counts: [{adata.X.min()}, {adata.X.max()}]")
    
    # Normalize to median counts
    total_counts = np.sum(adata.X, axis=1)
    median_counts = np.median(total_counts)
    logger.debug(f"Median counts: {median_counts}")
    
    sc.pp.normalize_total(adata, target_sum=median_counts)
    logger.debug(f"Data range after normalization: [{adata.X.min()}, {adata.X.max()}]")
    
    # Log transform
    adata.layers['logTPM'] = np.log1p(adata.X)
    logger.debug(f"Data range after log1p: [{adata.layers['logTPM'].min()}, {adata.layers['logTPM'].max()}]")
    
    return adata

def preprocess_data(adata: AnnData, params: Dict) -> AnnData:
    """
    Complete preprocessing pipeline.
    
    Workflow:
    1. Validate batches and filter cell types
    2. Select abundant genes
    3. Create and normalize one dataset:
       - Abundant genes dataset for batch effects analysis
    
    Args:
        adata: Input data
        params: Analysis parameters
        
    Returns:
        normalized data with abundant genes for batch scores,
    """
    # 1. Validate and filter cell types
    adata_ct_selection = validate_and_filter_celltypes(adata, params)
    
    # 2. Select abundant genes
    adata_ct_genes_selection, abundant_genes = select_abundant_genes(adata_ct_selection, params)
    
    # 3. Normalize both datasets
    adata_batch = normalize_data(adata_ct_genes_selection)  # abundant genes for batch scores
    
    return adata_batch