"""
Configuration parameters for batch effect analysis.

This module contains all the parameters used in the analysis,
including data processing parameters, analysis thresholds, and pathway databases.
"""

from typing import Dict, List

# Basic analysis parameters
DEFAULT_PARAMS = {
    # Column names in AnnData object
    'batch_key': 'batch',           # Column containing batch information in adata.obs
    'celltype_key': 'celltype',     # Column containing cell type information in adata.obs
    
    # Filtering parameters
    'min_cells': 20,                # Minimum cells per cell type
    'abundant_gene_threshold': 1,    # Minimum UMI/cell for abundant genes
    'min_abundant_genes': 300,      # Minimum number of abundant genes required

    # GSEA parameters
    'permutation_num': 1000,        # Number of permutations for GSEA
    'processes': 4,                 # Number of processes for parallel computing
    'min_size': 5,                 # Minimum gene set size
    'max_size': 1000,                # Maximum gene set size
    'seed': 817,                     # Random seed for reproducibility
    
    # Species information (used for pathway analysis)
    'species': 'human',             # Species ('human' or 'mouse')
}

# Pathway databases for different species used in enrichment analysis
PATHWAY_DBS = {
    'human': [
        'GO_Molecular_Function_2023',    
        'GO_Biological_Process_2023',    
        'MSigDB_Hallmark_2020',         
        'Reactome_2022',                
        'WikiPathway_2023_Human',      
        'KEGG_2021_Human',              
        'BioPlanet_2019'                
    ],
    'mouse': [
        'GO_Molecular_Function_2023',    
        'GO_Biological_Process_2023',    
        'MSigDB_Hallmark_2020',         
        'Reactome_2022',                
        'BioPlanet_2019'
    ]
}

# Validation ranges for parameters
VALID_RANGES = {
    'min_cells': (20, float('inf')),         # At least 20 cells per cell type
    'abundant_gene_threshold': (0.1, float('inf')),  # Reasonable UMI threshold
    'min_abundant_genes': (100, float('inf')),  # At least 100 abundant genes
    'permutation_num': (10, 10000),          # Reasonable range for permutations
    'processes': (1, 32),                    # Reasonable range for processes
    'min_size': (5, 100),                    # Reasonable range for gene set size
    'max_size': (50, 1000),                  # Reasonable range for gene set size
}

# Valid values for species 
VALID_SPECIES = {'human', 'mouse'}

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
    required_keys = {'batch_key', 'celltype_key', 'min_cells', 
                    'permutation_num', 'processes', 'min_size',
                    'max_size', 'species', 'seed'}
    missing_keys = required_keys - set(params.keys())
    if missing_keys:
        raise ValueError(f"Missing required parameters: {missing_keys}")
    
    # Validate numeric parameters
    for param, (min_val, max_val) in VALID_RANGES.items():
        if params[param] < min_val or params[param] > max_val:
            raise ValueError(
                f"Parameter {param} must be between {min_val} and {max_val}"
            )
    
    # Validate species
    if params['species'].lower() not in VALID_SPECIES:
        raise ValueError(
            f"Species must be one of {VALID_SPECIES}, got {params['species']}"
        )
    
    return True

def get_pathway_dbs(species: str) -> List[str]:
    """
    Get pathway databases for a specific species.
    
    Args:
        species: Species name ('human' or 'mouse')
        
    Returns:
        List of pathway database names
        
    Raises:
        ValueError: If species is not supported
    """
    species = species.lower()
    if species not in PATHWAY_DBS:
        raise ValueError(f"Unsupported species: {species}")
    return PATHWAY_DBS[species] 