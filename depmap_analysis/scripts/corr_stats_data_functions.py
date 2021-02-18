"""This file contains functions to extract data and data models to carry the
data """
import numpy as np
from indra_db.util.s3_path import S3Path
from numpy.core.multiarray import ndarray
from typing import Optional, List, Tuple
from pydantic import BaseModel


class HistData(BaseModel):
    """Holds a histogram and its associated bin edges"""
    histogram: ndarray
    bin_edges: ndarray

    def bin_witdh(self) -> float:
        return self.bin_edges[1] - self.bin_edges[0]

    def bar_placements(self) -> ndarray:
        return self.bin_edges[:-1] + self.bin_witdh()/2


class Histogram(BaseModel):
    """Holds histogram data and its cold storage location on S3"""
    name: str
    histogram_data: HistData


class IndividualHistograms(BaseModel):
    """These histograms are produced with:

    np.histogram(a=data, bins='auto')
    """
    all_x_corrs: Histogram
    avg_x_corrs: Histogram
    top_x_corrs: Histogram
    all_azb_corrs: Histogram
    azb_avg_corrs: Histogram
    all_azfb_corrs: Histogram
    azfb_avg_corrs: Histogram
    all_x_filtered_corrs: Histogram
    avg_x_filtered_corrs: Histogram
    all_reactome_corrs: Optional[Histogram] = None
    reactome_avg_corrs: Optional[Histogram] = None


class CombinedHistograms(BaseModel):
    """These histograms are produced with 'density=True' using np.histogram

    np.histogram(a=data, bins='auto', density=True)
    """
    azfb_avg_corrs: Histogram
    avg_x_filtered_corrs: Histogram

    azb_avg_corrs: Histogram
    avg_x_corrs: Histogram
    reactome_avg_corrs: Optional[Histogram] = None


class HistogramDirectory(BaseModel):
    """Root for all stored histograms and the s3 location for the json repr"""
    individual: IndividualHistograms
    combined: CombinedHistograms
    s3_location: str

    def get_s3_path(self) -> S3Path:
        """Return an S3Path object of the s3 cold storage location

        Returns
        -------
        S3Path
        """
        return S3Path.from_string(self.s3_location)


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


def get_combined_hists(results: Results) -> CombinedHistograms:
    res_dict = results.dict()
    ch_in = {}
    names = ['azfb_avg_corrs', 'avg_x_filtered_corrs', 'azb_avg_corrs',
             'avg_x_corrs', 'reactome_avg_corrs']

    for name, data in res_dict.items():
        if name in names and len(data) > 0:
            _data = [t[0] for t in data] if isinstance(data[0], tuple) \
                else data
            be, h = np.histogram(a=_data, bins='auto', density=True)
            ch_in[name] = Histogram(name=name,
                                    histogram_data=HistData(histogram=h,
                                                            bin_edges=be))
    return CombinedHistograms(**ch_in)


def get_individial_hists(results: Results) -> IndividualHistograms:
    res_dict = results.dict()

    ih_in = {}
    for name, data in res_dict.items():
        _data = [t[0] for t in data] if isinstance(data[0], tuple) else data
        be, h = np.histogram(a=_data, bins='auto')
        ih_in[name] = Histogram(name=name,
                                histogram_data=HistData(histogram=h,
                                                        bin_edges=be))
    return IndividualHistograms(**ih_in)


def get_hist_dir(results: Results):
    comb_hists: CombinedHistograms = get_combined_hists(results)
    ind_hists: IndividualHistograms = get_individial_hists(results)
    return HistogramDirectory(individual=ind_hists, combined=comb_hists)
