import os
import numpy as np
import matplotlib.pyplot as plt
from .compute_metric import compute_metric

__all__ = ["plot_metric"]

# ============================================================================ #


def plot_metric(experiment: list[dict],
                metric: str,
                variable: str,
                variant_key: str = None,
                variant_list: list[str] = None,
                **kwargs):
    """
    Plots a given metric as a function of a given variable set for an experiment.

    Parameters
    ----------
    experiment : list[dict]
        The experiment dictionary.
    metric : str
        The metric to plot.
    variable : str
        The variable to plot against.
    variant_key : str
        The key in the case dictionary that corresponds to the variant.
    variant_list : list[str]
        The list of variants to plot.
    """

    # set up the axis for a test plot
    fig, ax = plt.subplots()
    ax.set_title(f"{metric} vs {variable}")
    ax.set_xlabel(variable)
    ax.set_ylabel(metric)

    if variant_key is None:
        plot_metric_all(experiment, metric, variable, ax, **kwargs)
    else:
        plot_metric_variant(experiment, metric, variable, variant_key,
                            variant_list, ax, **kwargs)
        ax.legend(**kwargs.get("legend_options", {}))

    # Apply plot_style options
    ax.set(**kwargs.get("axes_options", {}))

    # Save the figure
    if kwargs.get("save_fig", False):
        fig_dir = kwargs.get("fig_dir", None)
        fig_name = kwargs.get("fig_name", None)

        if fig_dir is None or fig_name is None:
            raise ValueError(
                "fig_dir and fig_name must be provided to save the figure.")

        fig.savefig(os.path.join(fig_dir, fig_name))
        plt.close(fig)


# ============================================================================ #
# Plotting functions
# ============================================================================ #


def plot_metric_all(experiment: list[dict], metric: str, variable: str,
                    ax: plt.Axes, **kwargs):
    """
    Plots a given metric as a function of a given variable set for an
    experiment.

    Parameters
    ----------
    experiment : list[dict]
        The experiment dictionary.
    metric : str
        The metric to plot.
    variable : str
        The variable to plot against.
    """

    # Construct the x and y data
    x = np.array(
        [case[variable] for case in experiment if case[variable] is not None])
    y = np.zeros(len(x), dtype=float)

    # Loop through each case/case in the experiment
    idx = 0
    for case in experiment:

        # ---------------------------------------------------------------- #
        # Set up the data for the current case

        # Determine if the current case should be plotted
        if case[variable] is None:
            continue

        # ---------------------------------------------------------------- #
        # Construct the y value for the metric.

        try:
            y[idx] = compute_metric(case, metric, **kwargs)
        except Exception as e:
            print(f"Skipping case {case['name']}:\n\t {e}")
            y[idx] = np.nan

        # ---------------------------------------------------------------- #
        # Increment the index and continue to the next case
        idx += 1

    # ------------------------------------------------------------------------ #
    # Plot the data

    ax.plot(x, y, **kwargs.get("plot_options", {}))


def plot_metric_variant(experiment: list[dict], metric: str, variable: str,
                        variant_key: str, variant_list: list[str],
                        ax: plt.Axes, **kwargs):
    """
    Plots a given metric as a function of a given variable set for an
    experiment.

    Parameters
    ----------
    experiment : list[dict]
        The experiment dictionary.
    metric : str
        The metric to plot.
    variable : str
        The variable to plot against.
    variant_key : str
        The key in the case dictionary that corresponds to the variant.
    variant_list : list[str]
        The list of variants to plot.
    """

    for variant in variant_list:

        x = np.array([
            case[variable] for case in experiment
            if case[variable] is not None and case[variant_key] == variant
        ])
        y = np.zeros(len(x), dtype=float)

        # Loop through each case/case in the experiment
        idx = 0
        for case in experiment:

            # ---------------------------------------------------------------- #
            # Set up the data for the current case

            # Determine if the current case should be plotted
            if case[variable] is None:
                continue
            if not case[variant_key] == variant:
                continue

            # ---------------------------------------------------------------- #
            # Construct the y value for the metric.

            try:
                y[idx] = compute_metric(case, metric, **kwargs)
            except Exception as e:
                print(f"Skipping case {case['name']}:\n\t {e}")
                y[idx] = np.nan

            # ---------------------------------------------------------------- #
            # Increment the index and continue to the next case
            idx += 1

        # -------------------------------------------------------------------- #
        # Plot the data

        # Add to the existing plot
        ax.plot(x, y, label=variant, **kwargs.get("plot_options", {}))
