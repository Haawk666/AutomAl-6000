# Internal imports:
import core
import utils
# External imports:
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import copy
import time
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def column_detection(project, plot=False):
    """Experimental column detection"""
    pass


def plot_gamma(project):

    peak_gamma = []
    avg_gamma = []

    for vertex in project.graph.vertices:
        if not vertex.void and not vertex.is_edge_column:
            peak_gamma.append(vertex.peak_gamma)
            avg_gamma.append(vertex.avg_gamma)

    peak_gamma, avg_gamma = zip(*(sorted(zip(peak_gamma, avg_gamma), reverse=True)))

    fig = plt.figure(constrained_layout=True)
    gs = GridSpec(1, 1, figure=fig)
    ax_values = fig.add_subplot(gs[0, 0])

    ax_values.plot(
        range(0, len(peak_gamma)),
        peak_gamma,
        c='k',
        label='Peak gamma'
    )
    ax_values.plot(
        range(0, len(peak_gamma)),
        avg_gamma,
        c='r',
        label='Avg gamma'
    )
    ax_values.plot(
        range(0, len(peak_gamma)),
        [project.pixel_average] * len(peak_gamma),
        c='b',
        label='Average image pixel value'
    )
    ax_values.set_title('Column detection summary')
    ax_values.set_xlabel('# Column')
    ax_values.set_ylabel('Pixel intensity')
    ax_values.legend()

    plt.show()












