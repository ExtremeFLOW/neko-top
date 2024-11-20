"""
Compute a metric for a given case.

This script is used to compute a metric for a given case. The metric is
specified by the user and the case is provided as a dictionary. The script
returns the computed metric.
"""

# Global imports
import os

# Local imports
from .separation_angle import inflection_benchmark


def compute_metric(case: dict, metric: str, **kwargs) -> float:
    """
    Compute a metric for a given case.

    Parameters
    ----------
    case : dict
        The case dictionary.
    metric : str
        The metric to compute.
    kwargs : dict
        The keyword arguments.

    Returns
    -------
    float
        The computed metric.
    """

    if metric == "Frequency of Inflection Points":

        # Extract keyword arguments
        cache_dir = kwargs.get("cache_dir", None)
        cases_dir = kwargs.get("cases_dir", None)

        case_dir = os.path.join(cases_dir, case["name"])

        circ_radius = kwargs.get("circ_radius", None)
        if circ_radius is None:
            raise ValueError(f"circ_radius must be provided for {metric}")

        file_name = os.path.join(case_dir, 'circ_' + circ_radius + '.csv')

        benchmark_results = inflection_benchmark(file_name, cache_dir)

        match kwargs.get("which", "mean"):
            case "0":
                return benchmark_results["max_freq"][0]
            case "1":
                return benchmark_results["max_freq"][1]
            case "2":
                return benchmark_results["max_freq"][2]

    else:
        print(f"Metric {metric} not recognized.")
        exit(1)
