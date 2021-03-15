
# Internals
import core
# Externals
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import copy
import time
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
plt.rcParams.update({'font.size': 12})

base_filenames = [
    'test_set/prepared/0_Medium_Qprime',
    'test_set/prepared/0_multi_phase',
    'test_set/prepared/0_Small_L',
    'test_set/prepared/0_Smart_aligned_Qprime',
    'test_set/prepared/023',
    'test_set/prepared/c_Kenji_No1_exMgCu_64min_250C_009_IFFT'
]

# base_filenames = [
#     'test_set/prepared/0_Medium_Qprime',
#     'test_set/prepared/0_multi_phase',
#     'test_set/prepared/0_Small_L',
# ]

# base_filenames = [
#     'test_set/prepared/0_Small_L',
# ]

results_filename = 'test_set/test_results'

sequences = [
    [16, 24, 6, 9, 10, 11, 12, 13, 8, 9, 10, 11, 12, 13, 11, 9, 10, 20],
    [16, 11, 7, 6, 11, 9, 10, 12, 13, 11, 8, 6, 9, 10, 12, 13, 11, 20],
    [16, 11, 7, 6, 23, 12, 13, 9, 10, 11, 20]
]

poor_sequences = [
    [16, 11, 7, 6, 11, 12, 12, 13, 9, 10, 8, 6, 11, 12, 13, 6, 9, 10, 11, 20]
]

colors = [
    (0.9, 0.1, 0.1, 1.0),
    (0.1, 0.9, 0.1, 1.0),
    (0.1, 0.1, 0.9, 1.0),
    (0.9, 0.9, 0.1, 1.0),
    (0.9, 0.1, 0.9, 1.0),
    (0.1, 0.9, 0.9, 1.0)
]


def plot_results(results):

    grid_size = int(np.ceil(np.abs(np.sqrt(len(results)))))
    fig = plt.figure(constrained_layout=True, figsize=(6 * grid_size, 5 * grid_size), dpi=300)
    gs = GridSpec(grid_size, grid_size, figure=fig)
    ax = []
    row = 0
    column = 0
    if grid_size == 1:
        ax.append(fig.add_subplot(gs[0, 0]))
    elif grid_size == 2:
        ax.append(fig.add_subplot(gs[0, 0]))
        ax.append(fig.add_subplot(gs[0, 1]))
        ax.append(fig.add_subplot(gs[1, 0]))
        ax.append(fig.add_subplot(gs[1, 1]))
    else:
        for grid in range(grid_size ** 2):
            ax.append(fig.add_subplot(gs[row, column]))
            if np.mod(grid, grid_size - 1) == 0 and grid > 0:
                row += 1
                column = 0
            else:
                column += 1
    ax = ax[0:len(results)]
    ax_2 = []
    for axis in ax:
        ax_2.append(axis.twinx())

    for im, im_title in enumerate(base_filenames):

        for v, version in enumerate(results[im]):

            ax[im].plot(
                [i for i in range(len(sequences[v]))],
                [version[i][3] for i in range(len(sequences[v]))],
                c=colors[v],
                linestyle='solid',
                label='v{} (Error (%))'.format(v)
            )

            ax_2[im].plot(
                [i for i in range(len(sequences[v]))],
                [version[i][2] for i in range(len(sequences[v]))],
                c=colors[v],
                linestyle='dotted',
            )
            ax[im].plot(
                np.nan,
                c=colors[v],
                linestyle='dotted',
                label='v{} (chi)'.format(v)
            )

        ax[im].set_title(im_title)
        ax[im].set_ylabel('Error (%)')
        ax_2[im].set_ylabel('Chi (1)')
        ax[im].set_xlabel('Step (#)')
        ax[im].legend()

    fig.suptitle('Algorithm test results')

    plt.savefig(results_filename, bbox_inches='tight')


def compare(instance, control_instance):
    num_errors = 0
    num_precipitate_errors = 0
    num_columns = 0
    num_precipitate_columns = 0

    for vertex in instance.graph.vertices:
        if not vertex.is_edge_column:
            if not vertex.atomic_species == control_instance.graph.vertices[vertex.i].atomic_species:
                num_errors += 1
                if vertex.is_in_precipitate:
                    num_precipitate_errors += 1
                    num_precipitate_columns += 1
            else:
                if vertex.is_in_precipitate:
                    num_precipitate_columns += 1
            num_columns += 1

    if not num_columns == 0:
        percent_error = 100 * (num_errors / num_columns)
    else:
        percent_error = 100

    if not num_precipitate_columns == 0:
        percent_precipitate_error = 100 * (num_precipitate_errors / num_precipitate_columns)
    else:
        percent_precipitate_error = 100

    return [num_precipitate_errors, percent_precipitate_error, num_errors, percent_error]


def test_algorithm():

    logger.info('Testing algorithm versions...')

    time_all_1 = time.time()

    results = []

    for b, base_filename in enumerate(base_filenames):

        logger.info('    Analyzing image {}:'.format(base_filename))

        version_results = []

        init = core.Project.load(base_filename + '_init', supress_logging=True)
        control = core.Project.load(base_filename + '_control', supress_logging=True)

        for s, sequence in enumerate(sequences):

            logger.info('        With algorithm version {}'.format(s))

            version_result = []

            this_init = copy.deepcopy(init)
            total_time = 0

            for i, step in enumerate(sequence):

                time_1 = time.time()
                this_init.column_characterization_2(starting_index=this_init.starting_index, sub_method=step, supress_logging=True)
                time_2 = time.time()
                total_time += time_2 - time_1
                errors = compare(this_init, control)
                this_init.graph.calc_chi()
                version_result.append([step, total_time, this_init.graph.chi, errors[1], errors[3]])

            version_results.append(version_result)

            if results_filename:
                this_init.save(base_filename + '_result_v{}'.format(s), supress_logging=True)

        results.append(version_results)

    time_all_2 = time.time()
    summary = 'Testing complete in {:.1f} seconds. Summary:\n'.format(time_all_2 - time_all_1)
    for im, image_result in enumerate(results):
        summary += '    Image {}:\n'.format(im)
        summary += '                     Seconds      chi           total error      precipitate error\n'
        summary += '        -----------------------------------------------------------------------------------------\n'
        for v, version in enumerate(image_result):
            summary += '        Version {}: {:5.0f}         {:7.4f}      {:7.1f}          {:7.1f}\n'.format(v, version[-1][1], version[-1][2], version[-1][3], version[-1][4])

    logger.info(summary)

    if results_filename:
        plot_results(results)


