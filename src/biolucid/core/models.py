"""
Data models for storing and managing analysis results.

This module defines the data structures used to store both batch effect
and biological effect analysis results, along with utility methods for
result manipulation and summary generation.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional, List
import numpy as np
import pandas as pd

@dataclass
class BatchEffectResult:
    """
    Results from batch effect analysis.
    
    Attributes:
        b_score: Global proportion of batch effect in total variance (0-1)
        B_score: Global batch effect (non-negative)
        b_scores_per_batch: Per-batch proportion of batch effect
        B_scores_per_batch: Per-batch effect
        components: Component calculations Xmean, Xst, deltaS, deltaT
        residuals: Calculation residuals
    """
    b_score: float
    B_score: float
    b_scores_per_batch: Dict[str, float]
    B_scores_per_batch: Dict[str, float]
    components: Dict[str, pd.Series]
    residuals: pd.Series
    
    def __post_init__(self):
        """Validate results after initialization."""
        self._validate_scores()
    
    def _validate_scores(self):
        """Validate score values."""
        # Global score validation
        if not 0 <= self.b_score <= 1:
            raise ValueError("b_score must be between 0 and 1")
        if self.B_score < 0:
            raise ValueError("B_score must be non-negative")
            
        # Per-batch score validation
        for batch, score in self.b_scores_per_batch.items():
            if not 0 <= score <= 1:
                raise ValueError(f"b_score for batch {batch} must be between 0 and 1")
        for batch, score in self.B_scores_per_batch.items():
            if score < 0:
                raise ValueError(f"B_score for batch {batch} must be non-negative")
    
    def to_dict(self) -> Dict:
        """
        Convert results to dictionary format.
        
        Returns:
            Dictionary containing all results
        """
        return {
            'b_score': self.b_score,
            'B_score': self.B_score,
            'b_scores_per_batch': self.b_scores_per_batch,
            'B_scores_per_batch': self.B_scores_per_batch,
            'components': {
                k: v.to_dict() for k, v in self.components.items()
            },
            'residuals': self.residuals.to_dict()
        }

@dataclass
class BiologicalEffectResult:
    """
    Results from biological effect analysis.
    
    Attributes:
        Biological_summary_score: Global biological effect score
        Biological_batch_scores: Biological effect score per batch
        Biological_celltype_scores: Score per batch-celltype combination
        Biological_enrichment_results: Detailed enrichment results for each batch-celltype pair
    """
    Biological_summary_score: float
    Biological_batch_scores: Dict[str, float]
    Biological_celltype_scores: Dict[str, Dict[str, float]]
    Biological_enrichment_results: Dict[str, Dict[str, Optional[pd.DataFrame]]]
    
    def __post_init__(self):
        """Validate results after initialization."""
        if not 0 <= self.Biological_summary_score <= 1:
            raise ValueError("Biological summary score must be between 0 and 1")
            
        for batch, score in self.Biological_batch_scores.items():
            if not 0 <= score <= 1:
                raise ValueError(f"Score for batch {batch} must be between 0 and 1")
    
    def to_dict(self) -> Dict:
        """Convert results to dictionary format."""
        return {
            'Biological_summary_score': self.Biological_summary_score,
            'Biological_batch_scores': self.Biological_batch_scores,
            'Biological_celltype_scores': self.Biological_celltype_scores,
            'Biological_enrichment_results': {
                batch: {
                    celltype: df.to_dict() if df is not None else None
                    for celltype, df in celltype_results.items()
                }
                for batch, celltype_results in self.Biological_enrichment_results.items()
            }
        }

@dataclass
class AnalysisResult:
    """
    Complete analysis results including both batch and biological effects.
    
    Attributes:
        batch_effect: Batch effect analysis results
        bio_effect: Biological effect analysis results
        metadata: Analysis metadata (parameters, statistics)
    """
    batch_effect: BatchEffectResult
    bio_effect: BiologicalEffectResult
    metadata: Dict = field(default_factory=dict)
    
    def to_dict(self) -> Dict:
        """
        Convert complete results to dictionary format.
        
        Returns:
            Dictionary containing all results
        """
        return {
            'batch_effect': self.batch_effect.to_dict(),
            'bio_effect': self.bio_effect.to_dict(),
            'metadata': self.metadata
        } 