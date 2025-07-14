"""
Data models for storing and managing bioLUCID analysis results.
"""

from dataclasses import dataclass
from typing import Dict, Union
import pandas as pd

@dataclass
class bioLUCIDResult:
    """
    Results from bioLUCID analysis.
    
    Attributes:
        Global_b_score: Global proportion of batch effect in total variance (0-1)
        Global_q_sh_score: Global amount of shared gene expression variation
        Global_q_sp_score: Global amount of cell-type-specific gene expression variation
        b_score_per_batch: Per-batch proportion of batch effect
        q_sh_score_per_batch: Per-batch shared variation scores
        q_sp_score_per_batch: Per-batch cell-type-specific variation scores
        components: Component calculations μ, X_st, tau_t, sigma_s, residuals
    """
    # Global scores
    Global_b_score: float
    Global_q_sh_score: float
    Global_q_sp_score: float
    
    # Per-batch scores
    b_score_per_batch: Dict[str, float]
    q_sh_score_per_batch: Dict[str, float]
    q_sp_score_per_batch: Dict[str, float]
    
    # Components
    components: Dict[str, Union[pd.Series, pd.DataFrame]]
    
    def __post_init__(self):
        """Validate results after initialization."""
        self._validate_scores()
    
    def _validate_scores(self):
        """Validate score values."""
        # Global score validation
        if not 0 <= self.Global_b_score <= 1:
            raise ValueError("Global_b_score must be between 0 and 1")
        if self.Global_q_sh_score < 0:
            raise ValueError("Global_q_sh_score must be non-negative")
        if self.Global_q_sp_score < 0:
            raise ValueError("Global_q_sp_score must be non-negative")
            
        # Per-batch score validation
        for batch, score in self.b_score_per_batch.items():
            if not 0 <= score <= 1:
                raise ValueError(f"b_score for batch {batch} must be between 0 and 1")
        for batch, score in self.q_sh_score_per_batch.items():
            if score < 0:
                raise ValueError(f"q_sh_score for batch {batch} must be non-negative")
        for batch, score in self.q_sp_score_per_batch.items():
            if score < 0:
                raise ValueError(f"q_sp_score for batch {batch} must be non-negative")
    
    def to_dict(self) -> Dict:
        return {
            'Global_b_score': self.Global_b_score,
            'Global_q_sh_score': self.Global_q_sh_score,
            'Global_q_sp_score': self.Global_q_sp_score,
            'b_score_per_batch': self.b_score_per_batch,
            'q_sh_score_per_batch': self.q_sh_score_per_batch,
            'q_sp_score_per_batch': self.q_sp_score_per_batch,
            'components': {
                k: v.to_dict() if hasattr(v, 'to_dict') else v 
                for k, v in self.components.items()
            }
        }