
"""General utility/convenience functions"""

# External imports:
import numpy as np
from copy import deepcopy
from PIL import Image
import os.path
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def make_int_from_list(list_):
    """Concatenate the menbers of a list of integers into a integer."""
    return int('0'.join([str(i) for i in list_]))


def cyclic_sort(list_):
    """Find the lowest member of list and perform cyclic permuting so that the lowest member is first.

    Assumes unique members.

    """

    min_ = min(list_)
    min_index = -1
    for k, item in enumerate(list_):
        if item == min_:
            min_index = k
            break

    if min_index == -1:
        return list_
    elif min_index == 0:
        return list_
    else:
        tmp_lst = []
        for j in range(min_index, len(list_)):
            tmp_lst.append(list_[j])
        for j in range(0, min_index):
            tmp_lst.append(list_[j])
        return tmp_lst


def is_circularly_identical(list_1, list_2):
    """Returns a boolean indicating whether list_1 and list_2 are the same in any cyclic perturbation of each other.

    """
    if not len(list_1) == len(list_2):
        return False
    else:
        return ' '.join(map(str, list_2)) in ' '.join(map(str, list_1 * 2))


def circularize_next_index(i, i_max):
    """Make an index cyclic.

    parameters
    ----------
    i : int
        Integer index to be converted to a cyclic index
    i_max : int
        The maximum valid value of i.

    returns
    ----------
    i: int
        If input i is equal to i_max + 1, return 0, otherwise return the input i.

    note
    ----------
    In python, list[-1] will automatically give you the last element of a list, so checking for i == -1 is redundant.
    """

    if i == i_max + 1:
        return 0
    else:
        return i


def y_line(slope, b, x):

    return slope * x + b


def z_score(x, mean, var):
    z = (x - mean) / np.sqrt(var)
    return z


def normal_dist(x, mean, var):
    dist = (x - mean)**2
    if not var == 0:
        dist = dist / (2 * var)
    else:
        dist = 0
    dist = -dist
    dist = np.exp(dist)
    if not var == 0:
        dist = (1 / np.sqrt(2 * np.pi * var)) * dist
    else:
        dist = 0
    return dist


def multivariate_normal_dist(x, means, covar_determinant, inverse_covar_matrix):
    """All arguments assumed to be python lists"""
    k = len(x)
    deviance_vector = np.array([a - b for a, b in zip(x, means)])
    result = np.matmul(deviance_vector.T, np.array(inverse_covar_matrix))
    result = np.matmul(result, deviance_vector)
    result = np.exp(-0.5 * result)
    result = result / np.sqrt(((2 * np.pi) ** k) * covar_determinant)
    return float(result)


def normalize_list(in_list, norm_sum=1):
    factor = norm_sum / sum(in_list)
    norm_list = [a * factor for a in in_list]
    return norm_list


def normalize_dict(in_dict, norm_sum=1):
    sum = 0
    for value in in_dict.values():
        sum += value
    if sum == 0:
        new_dict = {}
        for key, value in in_dict.items():
            new_dict[key] = 1 / len(in_dict)
    else:
        factor = norm_sum / sum
        new_dict = {}
        for key, value in in_dict.items():
            new_dict[key] = value * factor
    return new_dict


def normalize_array(in_array, norm_sum=1):

    factor = in_array.sum() / norm_sum
    norm_array = in_array / factor
    return norm_array


def find_angle(a_1, a_2, b_1, b_2):
    # Returns the smallest angle between two vectors (a_1, b_1) and (a_2, b_2)
    alpha = a_1 * a_2 + b_1 * b_2
    alpha = alpha / np.sqrt((a_1**2 + b_1**2) * (a_2**2 + b_2**2))
    alpha = np.arccos(alpha)

    return alpha


def mean_val(data):
    """Return the mean of the given data

    :param data: list or 1-d np.array with data values
    :type data: [float] or np.array
    :returns mean:
    :rtype float:
    """

    if type(data) is type([1, 2, 3]):
        if not len(data) == 0:
            return sum(data) / len(data)
        else:
            return 0
    elif type(data) is type(np.array([1, 2, 3])):
        if not data.shape[0] == 0:
            return np.sum(data) / data.shape[0]
        else:
            return 0
    else:
        logger.warning('Unknown data format, returning 0..')
        return 0


def variance(data):
    """Return the var of the given data

    :param data: list or 1-d np.array with data values
    :type data: [float] or np.array
    :returns mean:
    :rtype float:
    """

    mean = mean_val(data)

    if type(data) is type([1, 2, 3]):
        sum_ = 0
        for item in data:
            sum_ += (item - mean) ** 2
        if len(data) > 1:
            var = sum_ / (len(data) - 1)
        else:
            var = 1
        return var
    elif type(data) is type(np.array([1, 2, 3])):
        var = 0
        for data_item in data:
            var += (data_item - mean) ** 2
        if not data.shape[0] < 2:
            var = var / (data.shape[0] - 1)
        else:
            var = 1
        return var
    else:
        logger.warning('Unknown data format, returning 1..')
        return 1


def covariance(data_1, data_2):
    """Two-pass numerically stable covariance calculation"""

    mean_1 = mean_val(data_1)
    mean_2 = mean_val(data_2)

    sum_ = 0
    for item_1, item_2 in zip(data_1, data_2):
        sum_ += (item_1 - mean_1) * (item_2 - mean_2)

    if len(data_1) > 2:
        covar = sum_ / (len(data_1) - 1)
    else:
        covar = 0

    return covar


def find_angle_from_points(p1, p2, pivot):

    # Always returns the anti-clockwise angle from p1 to p2 around the pivot (radians)

    vec_1 = (p1[0] - pivot[0], p1[1] - pivot[1])
    vec_2 = (p2[0] - pivot[0], p2[1] - pivot[1])
    alpha = find_angle(vec_1[0], vec_2[0], vec_1[1], vec_2[1])

    if vector_cross_product_magnitude(vec_1[0], vec_2[0], vec_1[1], vec_2[1]) > 0:
        alpha = 2 * np.pi - alpha

    return alpha


def vector_magnitude(vec):

    length = vec[0] ** 2 + vec[1] ** 2
    length = np.sqrt(length)

    return length


def vector_cross_product_magnitude(a_1, a_2, b_1, b_2):

    return a_1 * b_2 - b_1 * a_2


def side(a, b, c):
    """ Returns a position of the point c relative to the line going through a and b
        Points a, b are expected to be different
    """
    d = (c[1]-a[1])*(b[0]-a[0]) - (b[1]-a[1])*(c[0]-a[0])
    return 1 if d > 0 else (-1 if d < 0 else 0)


def is_point_in_closed_segment(a, b, c):
    """ Returns True if c is inside closed segment, False otherwise.
        a, b, c are expected to be collinear
    """
    if a[0] < b[0]:
        return a[0] <= c[0] and c[0] <= b[0]
    if b[0] < a[0]:
        return b[0] <= c[0] and c[0] <= a[0]

    if a[1] < b[1]:
        return a[1] <= c[1] and c[1] <= b[1]
    if b[1] < a[1]:
        return b[1] <= c[1] and c[1] <= a[1]

    return a[0] == c[0] and a[1] == c[1]


def closed_segment_intersect(a, b, c, d):
    """ Verifies if closed segments a, b, c, d do intersect.
    """
    if a == b or c == d or a == c or a == d or b == c or b == d:
        return False

    s1 = side(a, b, c)
    s2 = side(a, b, d)

    # All points are collinear
    if s1 == 0 and s2 == 0:
        return \
            is_point_in_closed_segment(a, b, c) or is_point_in_closed_segment(a, b, d) or \
            is_point_in_closed_segment(c, d, a) or is_point_in_closed_segment(c, d, b)

    # No touching and on the same side
    if s1 and s1 == s2:
        return False

    s1 = side(c, d, a)
    s2 = side(c, d, b)

    # No touching and on the same side
    if s1 and s1 == s2:
        return False

    return True


def dual_sort(list_1, list_2):

    temp_list_1 = np.ndarray([list_1.shape[0]], dtype=type(list_1))
    temp_list_2 = np.ndarray([list_1.shape[0]], dtype=type(list_2))

    for x in range(0, list_1.shape[0]):

        index = list_1.argmin()
        temp_list_1[x] = deepcopy(list_1[index])
        temp_list_2[x] = deepcopy(list_2[index])
        list_1[index] = list_1.max() + 1

    return temp_list_1, temp_list_2


def sort_neighbours(indices, distances, n):
    temp_max = distances.max()
    temp_indices = np.ndarray([n], dtype=np.int)
    temp_distances = np.ndarray([n], dtype=np.float64)

    for j in range(0, n):
        k = distances.argmin()
        temp_distances[j] = distances[k]
        temp_indices[j] = indices[k]
        distances[k] = temp_max + j + 1

    return temp_indices, temp_distances


def cm_fit(mat, x_0, y_0, r):
    # "Centre of mass" fit. Calculates the centre of mass of the pixels in the circle centred at (x_0, y_0) with radius
    # r.
    (y_range, x_range) = mat.shape
    total_mass = 0
    weighted_x_sum = 0
    weighted_y_sum = 0
    for x_i in range(x_0 - r, x_0 + r):
        for y_i in range(y_0 - r, y_0 + r):
            distance = np.sqrt((x_0 - x_i) ** 2 + (y_0 - y_i) ** 2)
            if distance <= r:
                if -1 < x_i < x_range and -1 < y_i < y_range:
                    total_mass += mat[y_i, x_i]
                    weighted_x_sum += mat[y_i, x_i]*x_i
                    weighted_y_sum += mat[y_i, x_i]*y_i
                else:
                    total_mass += 0.1
                    weighted_x_sum += 0.1*x_i
                    weighted_y_sum += 0.1*y_i
    x_fit = weighted_x_sum / total_mass
    y_fit = weighted_y_sum / total_mass
    return x_fit, y_fit


def gen_fft(mat):
    """Generate the fast Fourier transform of a matrix."""
    fourier = np.fft.fft2(mat)
    fourier = np.fft.fftshift(fourier)
    fourier = abs(fourier)
    fourier = np.log10(fourier)
    lowest = np.nanmin(fourier[np.isfinite(fourier)])
    highest = np.nanmax(fourier[np.isfinite(fourier)])
    original_range = highest - lowest
    fourier = (fourier - lowest) / original_range * 255
    return fourier


def normalize_static(mat):
    """Normalize matrix such that its elements are \in [0, 1]"""
    # General matrix norm to interval (0, 1). Assumed to be float types
    if mat.min() == mat.max() or mat.max() == 0.0:
        pass
    else:
        mat = mat - mat.min()
        mat = (mat / mat.max())
    return mat


def gen_framed_mat(mat, r):
    """Return input matrix padded with a frame of zeros and width r.

    :param mat: Matrix to be padded.
    :type mat: np.array
    :param r: Width of frame
    :type r: int
    :returns Framed matrix:
    :rtype np.array:

    """
    hor_frame = np.zeros((r, mat.shape[1]), dtype=type(mat))
    ver_frame = np.zeros((mat.shape[0] + 2 * r, r), dtype=type(mat))
    framed_mat = np.concatenate((hor_frame, mat, hor_frame), axis=0)
    framed_mat = np.concatenate((ver_frame, framed_mat, ver_frame), axis=1)
    return framed_mat


def gen_de_framed_mat(temp_mat, r):
    """Return input matrix de-framed by width r.

    :param temp_mat: Matrix to be de-padded.
    :type temp_mat: np.array
    :param r: Width of frame
    :type r: int
    :returns De-framed matrix:
    :rtype np.array:

    """
    (height, width) = temp_mat.shape
    return temp_mat[r:height - r, r:width - r]


def delete_pixels(mat, x_0, y_0, r):
    # Deletes (or sets elements to 0) in the (discrete)-circular area defined by it's centre (x_0, y_0) and radius r.
    (y_range, x_range) = mat.shape
    for x_i in range(x_0 - r, x_0 + r):
        for y_i in range(y_0 - r, y_0 + r):
            if -1 < x_i < x_range and -1 < y_i < y_range:
                if np.sqrt((x_0 - x_i)**2 + (y_0 - y_i)**2) <= r:
                    mat[y_i, x_i] = 0.0
    return mat


def draw_circle(mat, x_0, y_0, r):
    # Draws a simple pixel-wide ring in the input matrix centered at (x_0, y_0) with radius r. The input matrix is
    # assumed to be sparse, that is, mostly contains 0's.
    for x_i in range(x_0 - r, x_0 + r):
        for y_i in range(y_0 - r, y_0 + r):
            distance = np.sqrt((x_0 - x_i) ** 2 + (y_0 - y_i) ** 2)
            if np.floor(distance) == r:
                mat[y_i, x_i] = 1
    return mat


def circular_average(mat, x_0, y_0, r):
    # Calculate the average value of the elements in the circle centered at (x_0, y_0) with radius r.
    (y_range, x_range) = mat.shape
    element_sum = 0.0
    peak_element = 0.0
    counter = 0
    for x_i in range(x_0 - r, x_0 + r):
        for y_i in range(y_0 - r, y_0 + r):
            distance = np.sqrt((x_0 - x_i) ** 2 + (y_0 - y_i) ** 2)
            if distance <= r:
                if -1 < x_i < x_range and -1 < y_i < y_range:
                    element_sum = element_sum + mat[y_i, x_i]
                    counter = counter + 1
                    if mat[y_i, x_i] > peak_element:
                        peak_element = mat[y_i, x_i]
    element_average = element_sum / counter
    return element_average, peak_element


def find_index_from_coor(mat, x, y):
    return mat[y, x, 1]


def find_coor_from_index(mat, i):

    (M, N) = mat.shape()
    for x in range(0, N):
        for y in range(0, M):
            if mat[y, x, 1] == i:
                return x, y
    print('Error?')


def im_out_static(im_mat, filename_full):
    im_mat = (im_mat - np.min(im_mat)) / (np.max(im_mat) - np.min(im_mat))
    im_mat = np.uint8(np.round(im_mat * 255))
    im_dsp = Image.fromarray(im_mat)
    im_dsp.save(os.path.join(filename_full))


def norm_rgb_tuple(rgb):
    return rgb[0] / 255, rgb[1] / 255, rgb[2] / 255


def rgb_to_hex(rgb):
    return '#{0:02x}{1:02x}{2:02x}'.format(max(0, min(rgb[0], 255)), max(0, min(rgb[1], 255)), max(0, min(rgb[2], 255)))


def hex_to_rgb(hex):
    pass


