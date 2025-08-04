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
from ..utils.preprocessing import setup_logging, preprocess_data
from .calculator import bioLUCID_calculation
from .models import bioLUCIDResult

logger = logging.getLogger(__name__)

class BatchEffectAnalyzer:
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
        adata_batch: Processed data with abundant genes for batch effect analysis
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
        self.adata_batch = None    # For batch effect analysis
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

    
    def preprocess_data(self) -> AnnData:
        """
        Run data preprocessing.
        
        Returns:
            Tuple of (biological analysis data, batch effect analysis data)
        """
        logger.info("Starting data preprocessing")
        
        # Run preprocessing pipeline
        self.adata_batch= preprocess_data(
            self.adata_raw, 
            self.params
        )
        logger.info("Data preprocessing completed")
        # Batch effect related logs
        logger.info("=== Batch Effect Analysis Data ===")
        logger.info(f"Raw data: {self.adata_raw}")
        logger.info(f"Processed data: {self.adata_batch}")
        
        return self.adata_batch
    
    def run_analysis(self) -> bioLUCIDResult:
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
        
        try:
            # 1. Preprocess data if not done
            if self.adata_batch is None:
                logger.info("Preprocessing data")
                self.preprocess_data()
            
            # 2. Run bioLUCID using abundant genes
            bioLUCID_result = bioLUCID_calculation(self.adata_batch, self.params)

            # 3. Compile results
            self.results = bioLUCID_result
            
            logger.info("Analysis pipeline completed")
            return self.results
            
        except Exception as e:
            logger.error(f"Analysis failed: {str(e)}")
            raise