"""
Main analyzer class for batch effect analysis.

This module provides the main interface for running batch effect analysis,
coordinating data preprocessing, calculation, and result collection.
"""

import logging
from typing import Optional, Dict, Tuple
from anndata import AnnData
import time

from ..config.settings import DEFAULT_PARAMS, validate_params
from ..utils.preprocessing import setup_logging, preprocess_data, get_data_summary
from .calculator import calculate_batch_effect, calculate_biological_effect
from .models import AnalysisResult, BatchEffectResult, BiologicalEffectResult

logger = logging.getLogger(__name__)

class Analyzer:
    """
    Main class for batch effect analysis.
    
    This class coordinates the complete analysis pipeline, including:
    - Data preprocessing
    - Batch effect calculation
    - Biological effect calculation
    - Result management
    
    Attributes:
        adata_raw: Original input data
        params: Analysis parameters
        adata_bio: Processed data with all genes for biological effect analysis
        adata_batch: Processed data with abundant genes for batch effect analysis
        abundant_genes: List of indices for abundant genes
        results: Analysis results
    """
    
    def __init__(
        self,
        adata: AnnData,
        params: Optional[Dict] = None,
        verbose: bool = True
    ):
        """
        Initialize the analyzer.
        
        Args:
            adata: Input data
            params: Custom parameters (optional)
            verbose: Whether to show detailed progress
        """
        # Store original data
        self.adata_raw = adata
        
        # Initialize parameters
        self.params = DEFAULT_PARAMS.copy()
        if params is not None:
            self.params.update(params)
        
        # Initialize data containers
        self.adata_bio = None      # For biological effect analysis
        self.adata_batch = None    # For batch effect analysis
        self.abundant_genes = None # Indices of abundant genes
        self.results = None        # Analysis results
        
        # Set up logging
        setup_logging(verbose)
        
        # Validate parameters
        validate_params(self.params)
        logger.info("Initialized analyzer with valid parameters")

        # Print parameters
        logger.info("=== Analysis Parameters ===")
        for key, value in self.params.items():
            logger.info(f"{key}: {value}")

    
    def preprocess_data(self) -> Tuple[AnnData, AnnData]:
        """
        Run data preprocessing.
        
        Returns:
            Tuple of (biological analysis data, batch effect analysis data)
        """
        logger.info("Starting data preprocessing")
        
        # Record original dimensions
        original_dims = get_data_summary(self.adata_raw, self.params)
        
        # Run preprocessing pipeline
        self.adata_bio, self.adata_batch, self.abundant_genes = preprocess_data(
            self.adata_raw, 
            self.params
        )
        
        # Record processed dimensions
        bio_dims = get_data_summary(self.adata_bio, self.params)
        batch_dims = get_data_summary(self.adata_batch, self.params)
        
        logger.info("Data preprocessing completed")
        
        # Biological effect related logs
        logger.info("=== Biological Effect Analysis Data ===")
        logger.info(f"Retained {len(self.adata_bio)}/{len(self.adata_raw)} cells")
        logger.info(f"Original dimensions: {original_dims}")
        logger.info(f"Biological analysis dimensions: {bio_dims}")
        
        # Batch effect related logs
        logger.info("=== Batch Effect Analysis Data ===")
        logger.info(f"Selected {len(self.abundant_genes)} abundant genes")
        logger.info(f"Batch effect analysis dimensions: {batch_dims}")
        
        
        return self.adata_bio, self.adata_batch
    
    def run_analysis(self) -> AnalysisResult:
        """
        Run complete analysis pipeline.
        
        Workflow:
        1. Preprocess data and create two datasets:
           - Full dataset for biological scores
           - Abundant genes dataset for batch effects
        2. Calculate batch effects using abundant genes
        3. Calculate biological effects using all genes
        
        Returns:
            Complete analysis results
            
        Raises:
            ValueError: If analysis fails
        """
        logger.info("Starting analysis pipeline")
        start_time = time.time()
        
        try:
            # 1. Preprocess data if not done
            if self.adata_bio is None or self.adata_batch is None:
                self.preprocess_data()
            
            # 2. Calculate batch effect using abundant genes
            logger.info(f"Calculating batch effect scores using {len(self.abundant_genes)} abundant genes")
            batch_start = time.time()
            batch_effect = calculate_batch_effect(self.adata_batch, self.params)
            batch_time = time.time() - batch_start
            
            # 3. Calculate biological effect using all genes
            logger.info(f"Calculating biological effect scores using all {self.adata_bio.n_vars} genes")
            bio_start = time.time()
            bio_effect = calculate_biological_effect(self.adata_bio, self.params)
            bio_time = time.time() - bio_start
            
            # Compile results
            self.results = AnalysisResult(
                batch_effect=batch_effect,
                bio_effect=bio_effect,
                metadata={
                    'parameters': self.params,
                    'dimensions': {
                        'original': get_data_summary(self.adata_raw, self.params),
                        'biological': get_data_summary(self.adata_bio, self.params),
                        'batch_effect': get_data_summary(self.adata_batch, self.params),
                        'abundant_genes': len(self.abundant_genes)
                    },
                    'computation_time': {
                        'total': time.time() - start_time,
                        'batch_effect': batch_time,
                        'bio_effect': bio_time
                    }
                }
            )
            
            logger.info("Analysis pipeline completed")
            return self.results
            
        except Exception as e:
            logger.error(f"Analysis failed: {str(e)}")
            raise
    
    def get_summary(self) -> Dict:
        """
        Get analysis summary.
        
        Returns:
            Dictionary containing summary statistics
            
        Raises:
            ValueError: If analysis hasn't been run
        """
        if self.results is None:
            raise ValueError("No results available. Run analysis first.")
        return self.results.get_summary()
    
    def get_batch_effect_scores(self) -> Dict[str, float]:
        """
        Get batch effect scores.
        
        Returns:
            Dictionary with global and per-batch scores
            
        Raises:
            ValueError: If analysis hasn't been run
        """
        if self.results is None:
            raise ValueError("No results available. Run analysis first.")
            
        return {
            'global_score': self.results.batch_effect.b_score,
            'per_batch': self.results.batch_effect.b_scores_per_batch
        }
    
    def get_biological_effect_scores(self) -> Dict[str, float]:
        """
        Get biological effect scores.
        
        Returns:
            Dictionary with global and per-batch scores
            
        Raises:
            ValueError: If analysis hasn't been run
        """
        if self.results is None:
            raise ValueError("No results available. Run analysis first.")
            
        return {
            'global_score': self.results.bio_effect.summary_score,
            'per_batch': self.results.bio_effect.batch_scores
        }
    
    def get_celltype_scores(self) -> Dict[str, Dict[str, float]]:
        """
        Get per-celltype scores.
        
        Returns:
            Nested dictionary with scores for each batch-celltype pair
            
        Raises:
            ValueError: If analysis hasn't been run
        """
        if self.results is None:
            raise ValueError("No results available. Run analysis first.")
        return self.results.bio_effect.celltype_scores
    
    def get_key_results(self) -> Dict:
        """
        Get key analysis results.
        """
        if self.results is None:
            raise ValueError("No results available. Run analysis first.")
    
        return {
            # Batch effect results
            'b_score': self.results.batch_effect.b_score,
            'B_score': self.results.batch_effect.B_score,
            'b_scores_per_batch': self.results.batch_effect.b_scores_per_batch,
            'B_scores_per_batch': self.results.batch_effect.B_scores_per_batch,
            'components': self.results.batch_effect.components,
            'residuals': self.results.batch_effect.residuals,
    
            # Biological effect results
            'Biological_summary_score': self.results.bio_effect.Biological_summary_score,
            'Biological_batch_scores': self.results.bio_effect.Biological_batch_scores,
            'Biological_celltype_scores': self.results.bio_effect.Biological_celltype_scores,
            'Biological_enrichment_results': self.results.bio_effect.Biological_enrichment_results
        } 