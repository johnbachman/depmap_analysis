"""
This file contains functions to extract data and data models to carry the
data for the depmap script
"""
from typing import List, Tuple
from pydantic import BaseModel


class Results(BaseModel):
    """The results data model"""
    all_x_corrs: List[float] = []
    avg_x_corrs: List[float] = []
    top_x_corrs: List[Tuple[str, str, float]] = []
    all_azb_corrs: List[float] = []
    azb_avg_corrs: List[float] = []
    all_azfb_corrs: List[float] = []  # Background for filtered A,B
    azfb_avg_corrs: List[float] = []
    all_reactome_corrs: List[float] = []
    reactome_avg_corrs: List[float] = []
    all_x_filtered_corrs: List[float] = []
    avg_x_filtered_corrs: List[float] = []
