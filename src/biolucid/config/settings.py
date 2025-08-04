"""
Configuration parameters for batch effect analysis.

This module contains all the parameters used in the analysis,
including data processing parameters, analysis thresholds.
"""

from typing import Dict

# Basic analysis parameters
DEFAULT_PARAMS = {
    # Column names in AnnData object
    'batch_key': 'batch',           # Column containing batch information in adata.obs
    'celltype_key': 'celltype',     # Column containing cell type information in adata.obs
    
    # Filtering parameters
    'min_cells': 20,                # Minimum cells per cell type
    'abundant_gene_threshold': 1,    # Minimum UMI/cell for abundant genes
    'min_abundant_genes': 100,      # Minimum number of abundant genes required
}


# Validation ranges for parameters
VALID_RANGES = {
    'min_cells': (20, float('inf')),         # At least 20 cells per cell type
    'abundant_gene_threshold': (0.1, float('inf')),  # Reasonable UMI threshold
    'min_abundant_genes': (100, float('inf'))  # At least 100 abundant genes
}


def validate_params(params: Dict) -> bool:
    """
    Validate configuration parameters.
    
    Args:
        params: Dictionary containing configuration parameters
        
    Returns:
        bool: True if all parameters are valid
         
    Raises:
        ValueError: If any parameter is invalid
    """
    # Check required keys
    required_keys = {'batch_key', 'celltype_key', 'min_cells'}

    missing_keys = required_keys - set(params.keys())

    if missing_keys:
        raise ValueError(f"Missing required parameters: {missing_keys}")
    
    # Validate numeric parameters
    for param, (min_val, max_val) in VALID_RANGES.items():
        if params[param] < min_val or params[param] > max_val:
            raise ValueError(
                f"Parameter {param} must be between {min_val} and {max_val}"
            )
    
    return True