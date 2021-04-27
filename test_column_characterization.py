
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
    'test_set/0_Medium_Qprime',
    'test_set/0_multi_phase',
    'test_set/0_Small_L',
    'test_set/008',
    'test_set/012a',
    'test_set/018_Ba',
    'test_set/023',
    'test_set/030',
    'test_set/a_ARM_FICAL_6005_185C_2h_WQ_045_IFFT',
    'test_set/b_ARM_FICAL_150C_20h_WQ_006_IFFT',
    'test_set/c_Kenji_No1_exMgCu_64min_250C_009_IFFT',
    'test_set/haakon_016'
]

results_filename = 'test_set/results/test_results'

sequences = [
    [11, 7, 12, 13, 9, 10, 11, 8, 12, 13, 9, 10, 11, 8, 12, 13, 9, 10, 11, 8, 12, 13, 9, 10, 17, 11, 20]
]

colors = [
    (0.9, 0.1, 0.1, 1.0),
    (0.1, 0.9, 0.1, 1.0),
    (0.2, 0.2, 0.9, 1.0),
    (0.9, 0.9, 0.1, 1.0),
    (0.9, 0.1, 0.9, 1.0),
    (0.1, 0.9, 0.9, 1.0),
    (0.1, 0.1, 0.1, 1.0),
    (0.5, 0.5, 0.5, 1.0),
    (0.9, 0.5, 0.1, 1.0),
    (0.6, 0.5, 0.7, 1.0),
    (0.8, 0.9, 0.3, 1.0),
    (0.1, 0.3, 0.9, 1.0)
]


def plot_results(results):

    grid_size = int(np.ceil(np.abs(np.sqrt(len(results)))))
    fig = plt.figure(constrained_layout=True, figsize=(6 * grid_size, 5 * grid_size), dpi=400)
    gs = GridSpec(grid_size, grid_size, figure=fig)
    ax = []
    row = 0
    column = 0
    for grid in range(grid_size ** 2):
        ax.append(fig.add_subplot(gs[row, column]))
        if column == grid_size - 1:
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

        ax[im].vlines([0, 6, 12, 18, 24], 0, 100, linestyles='dashed', colors='k')
        ax[im].set_title(im_title)
        ax[im].set_ylim([0, 100])
        ax[im].set_ylabel('Error (%)')
        ax_2[im].set_ylabel('Chi (1)')
        ax_2[im].set_ylim([0, 0.15])
        ax[im].set_xlabel('Step (#)')
        ax[im].legend()

    fig.suptitle('Algorithm test - results pr image')

    plt.savefig(results_filename + '_pr_image', bbox_inches='tight')

    grid_size = int(np.ceil(np.abs(np.sqrt(len(sequences)))))
    fig = plt.figure(constrained_layout=True, figsize=(6 * grid_size, 5 * grid_size), dpi=400)
    gs = GridSpec(grid_size, grid_size, figure=fig)
    ax = []
    row = 0
    column = 0
    for grid in range(grid_size ** 2):
        ax.append(fig.add_subplot(gs[row, column]))
        if column == grid_size - 1:
            row += 1
            column = 0
        else:
            column += 1
    ax = ax[0:len(sequences)]
    ax_2 = []
    for axis in ax:
        ax_2.append(axis.twinx())

    for v, version in enumerate(sequences):

        for im, im_title in enumerate(base_filenames):

            ax[v].plot(
                [i for i in range(len(version))],
                [results[im][v][i][3] for i in range(len(version))],
                c=colors[im],
                linestyle='solid',
                label='image {}'.format(im)
            )

            ax_2[v].plot(
                [i for i in range(len(version))],
                [results[im][v][i][2] for i in range(len(version))],
                c=colors[im],
                linestyle='dotted',
            )
            # ax[v].plot(
            #     np.nan,
            #     c=colors[im],
            #     linestyle='dotted',
            #     label='im{} (chi)'.format(im)
            # )

        ax[v].vlines([0, 6, 12, 18, 24], 0, 100, linestyles='dashed', colors='k')
        ax[v].set_title('Version {}'.format(v))
        ax[v].set_ylabel('Error (%)')
        ax[v].set_ylim([0, 100])
        ax_2[v].set_ylabel('Chi (1)')
        ax_2[v].set_ylim([0, 0.15])
        ax[v].set_xlabel('Step (#)')
        ax[v].legend()

    plt.savefig(results_filename + '_pr_version', bbox_inches='tight')

    grid_size = int(np.ceil(np.abs(np.sqrt(len(sequences)))))
    fig = plt.figure(constrained_layout=True, figsize=(6 * grid_size, 5 * grid_size), dpi=400)
    gs = GridSpec(grid_size, grid_size, figure=fig)
    ax = []
    row = 0
    column = 0
    for grid in range(grid_size ** 2):
        ax.append(fig.add_subplot(gs[row, column]))
        if column == grid_size - 1:
            row += 1
            column = 0
        else:
            column += 1
    ax = ax[0:len(sequences)]

    for v, version in enumerate(sequences):

        for im, im_title in enumerate(base_filenames):
            ax[v].plot(
                [i for i in range(len(version))],
                [results[im][v][i][3] for i in range(len(version))],
                c=colors[im],
                linestyle='solid',
                label='im{} (Error (%))'.format(im)
            )

        ax[v].vlines([0, 6, 12, 18, 24], 0, 100, linestyles='dashed', colors='k')
        ax[v].set_ylim([0, 50])
        ax[v].set_xticks([0, 6, 12, 18, 24], False)
        ax[v].set_xticks([1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 25, 26], True)

    plt.savefig(results_filename + '_pr_version_error_only', bbox_inches='tight')


def compare(instance, control_instance):
    num_errors = 0
    num_precipitate_errors = 0
    num_columns = 0
    num_precipitate_columns = 0
    num_zeta_errors = 0
    num_precipitate_zeta_errors = 0

    for vertex in instance.graph.vertices:
        if not vertex.is_edge_column:
            if not vertex.zeta == control_instance.graph.vertices[vertex.i].zeta:
                num_zeta_errors += 1
                if vertex.is_in_precipitate:
                    num_precipitate_zeta_errors += 1
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
        percent_zeta_error = 100 * (num_zeta_errors / num_columns)
        if percent_zeta_error > 50:
            percent_zeta_error = 100 - percent_zeta_error
            num_zeta_errors = num_columns - num_zeta_errors
    else:
        percent_error = 100
        percent_zeta_error = 100

    if not num_precipitate_columns == 0:
        percent_precipitate_error = 100 * (num_precipitate_errors / num_precipitate_columns)
        percent_precipitate_zeta_error = 100 * (num_precipitate_zeta_errors / num_precipitate_columns)
        if percent_precipitate_zeta_error > 50:
            percent_precipitate_zeta_error = 100 - percent_precipitate_zeta_error
            num_precipitate_zeta_errors = num_precipitate_columns - num_precipitate_zeta_errors
    else:
        percent_precipitate_error = 100
        percent_precipitate_zeta_error = 100

    return [num_precipitate_errors, percent_precipitate_error, num_errors, percent_error, num_precipitate_zeta_errors,
            percent_precipitate_zeta_error, num_zeta_errors, percent_zeta_error]


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
                version_result.append([step, total_time, this_init.graph.chi, errors[1], errors[3], errors[5], errors[7]])

            version_results.append(version_result)

            if results_filename:
                pass
                # this_init.save(base_filename + '_result_v{}'.format(s), supress_logging=True)

        results.append(version_results)

    time_all_2 = time.time()
    avg_precipitate_error_percent = [0] * len(sequences)
    var_precipitate_error_percent = [0] * len(sequences)
    for image_result in results:
        for v, version in enumerate(image_result):
            avg_precipitate_error_percent[v] += version[-1][3]
    for image_result in results:
        for v, version in enumerate(image_result):
            var_precipitate_error_percent[v] += (version[-1][3] - avg_precipitate_error_percent[v] / len(results)) ** 2
    summary = 'Testing complete in {:.1f} seconds.\n'.format(time_all_2 - time_all_1)
    summary += '\n    Average algorithm performance (average precipitate error percent):\n'
    for v, version in enumerate(sequences):
        summary += '        Version {}: avg: {:7.1f}    var: {:7.1f}\n'.format(v, avg_precipitate_error_percent[v] / len(results), var_precipitate_error_percent[v] / (len(results) - 1))
    summary += '\n    All results:\n'
    for im, image_result in enumerate(results):
        summary += '    Image {}:\n'.format(base_filenames[im])
        summary += '                     Seconds      chi             total error      precipitate error       zeta error\n'
        summary += '        -------------------------------------------------------------------------------------------------------\n'
        for v, version in enumerate(image_result):
            summary += '        Version {}: {:5.0f}         {:7.4f}      {:7.1f}          {:7.1f}                 {:7.1f}\n'.format(v, version[-1][1], version[-1][2], version[-1][4], version[-1][3], version[-1][6])

    logger.info(summary)

    if results_filename:
        with open(results_filename + '.txt', 'w') as f:
            f.write(summary)
        plot_results(results)



