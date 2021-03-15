

# Internal imports
import core
import utils
import data_module
# External imports
import csv
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
import numpy as np
import pickle
import copy
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def test_vertex_data_manager(data=None):
    if data is None:
        data = [
            {'x': 2.5, 'y': 2.4, 'advanced_species': 'Si_1'},
            {'x': 0.5, 'y': 0.7, 'advanced_species': 'Si_1'},
            {'x': 2.2, 'y': 2.9, 'advanced_species': 'Si_1'},
            {'x': 1.9, 'y': 2.2, 'advanced_species': 'Si_1'},
            {'x': 3.1, 'y': 3.0, 'advanced_species': 'Si_1'},
            {'x': 2.3, 'y': 2.7, 'advanced_species': 'Si_1'},
            {'x': 2.0, 'y': 1.6, 'advanced_species': 'Si_1'},
            {'x': 1.0, 'y': 1.1, 'advanced_species': 'Si_1'},
            {'x': 1.5, 'y': 1.6, 'advanced_species': 'Si_1'},
            {'x': 1.1, 'y': 0.9, 'advanced_species': 'Si_1'}
        ]

    data_manager_object = data_module.VertexDataManager(
        '',
        attr_keys=['x', 'y'],
        category_key='advanced_species',
        save_filename='test'
    )

    data_manager_object.original_dict_data = data
    data_manager_object.category_list = ['Si_1']
    data_manager_object.n = len(data_manager_object.original_dict_data)

    data_manager_object.num_data_categories = len(data_manager_object.category_list)

    data_manager_object.matrix_data = data_manager_object.vectorize_data()
    data_manager_object.composite_model = []
    for c, category_data in enumerate(data_manager_object.matrix_data):
        data_manager_object.composite_model.append(
            data_module.MultivariateNormalDist(category_data, data_manager_object.category_list[c], data_manager_object.attribute_keys))

    data_manager_object.concatenated_matrix_data = data_manager_object.concatenate_categories()
    data_manager_object.uncategorized_normal_dist = data_module.MultivariateNormalDist(data_manager_object.concatenated_matrix_data, 'All categories',
                                                                                       data_manager_object.attribute_keys)

    data_manager_object.normalized_concatenated_matrix_data = np.array(data_manager_object.concatenated_matrix_data)
    data_manager_object.normalized_matrix_data = []
    for category_data in data_manager_object.matrix_data:
        data_manager_object.normalized_matrix_data.append(np.array(category_data))
    data_manager_object.norm_data()
    data_manager_object.normalized_uncategorized_normal_dist = data_module.MultivariateNormalDist(data_manager_object.normalized_concatenated_matrix_data,
                                                                       'All categories', data_manager_object.attribute_keys)
    data_manager_object.pca_feature_vector = data_manager_object.normalized_uncategorized_normal_dist.covar_matrix_eigenvectors.T
    data_manager_object.composed_uncategorized_data = np.matmul(data_manager_object.pca_feature_vector, data_manager_object.normalized_concatenated_matrix_data)
    data_manager_object.composed_data = []
    for normalized_category_data in data_manager_object.normalized_matrix_data:
        data_manager_object.composed_data.append(np.matmul(data_manager_object.pca_feature_vector, normalized_category_data))
    data_manager_object.composed_uncategorized_normal_dist = data_module.MultivariateNormalDist(data_manager_object.composed_uncategorized_data, 'Column',
                                                                                                data_manager_object.pc_keys)
    data_manager_object.composed_normal_dist = []
    for c, composed_category_data in enumerate(data_manager_object.composed_data):
        data_manager_object.composed_normal_dist.append(
            data_module.MultivariateNormalDist(composed_category_data, data_manager_object.category_list[c], data_manager_object.pc_keys))

    data_manager_object.dual_plot('x', 'y')

    data_manager_object.composite_model[0].get_max_probability()

    print('Covariance matrix: {}\n'.format(data_manager_object.uncategorized_normal_dist.covar_matrix))
    print('Covariance determinant: {}\n'.format(data_manager_object.uncategorized_normal_dist.covar_matrix_determinant))
    print('Inverse covariance matrix: {}\n'.format(data_manager_object.uncategorized_normal_dist.inverse_covar_matrix))
    print('Means: {}\n'.format(data_manager_object.uncategorized_normal_dist.means))
    print('covar_matrix_eigenvalues: {}\n'.format(data_manager_object.uncategorized_normal_dist.covar_matrix_eigenvalues))
    print('covar_matrix_eigenvectors: {}\n'.format(data_manager_object.uncategorized_normal_dist.covar_matrix_eigenvectors))

    data_manager_object.save()

    print(data_manager_object.calc_prediction({'x': 1.0, 'y': 8.2}))

    print(data_manager_object.report())


def test_statistical_basics():
    x = [1.5, 1.1, 2.4, 2.1, 3.7, 4.6, 4.2, 6.4, 6.2, 7.5]
    y = [0.9, 2.7, 2.0, 2.6, 2.5, 2.1, 4.6, 4.5, 6.4, 6.3]

    means = [sum(x) / len(x), sum(y) / len(y)]
    print(means)
    vars_ = [sum([(a - means[0]) ** 2 for a in x]) / 9, sum([(a - means[1]) ** 2 for a in y]) / 9]
    print(vars_)
    covars = [
        [vars_[0], sum([(a - means[0]) * (b - means[1]) for a, b in zip(x, y)]) / 9],
        [sum([(a - means[1]) * (b - means[0]) for a, b in zip(y, x)]) / 9, vars_[1]]
    ]
    print(covars)

    determinant = covars[0][0] * covars[1][1] - covars[0][1] * covars[1][0]
    print(determinant)

    inverse = [
        [covars[1][1] / determinant, -covars[0][1] / determinant],
        [-covars[1][0] / determinant, covars[0][0] / determinant]
    ]
    print(inverse)
    print('theoretical max: {}'.format(1/np.sqrt(((2 * np.pi) ** 2) * determinant)))
    print('calculated max: {}'.format(utils.multivariate_normal_dist(means, means, determinant, inverse)))

    data = []
    for a, b in zip(x, y):
        data.append({
            'x': a, 'y': b, 'advanced_species': 'Si_1'
        })
    test_vertex_data_manager(data)


def test_grid_population(grid_length):

    grid_size = int(np.ceil(np.abs(np.sqrt(grid_length))))
    print('grid length {} | grid size {}'.format(grid_length, grid_size))

    row = 0
    column = 0

    for grid in range(grid_size ** 2):
        print('grid {} | row {} | column {}'.format(grid, row, column))
        if column == grid_size - 1:
            row += 1
            column = 0
        else:
            column += 1



