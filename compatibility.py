"""Module to handle forwards compatibility between versions.

.. note::

    Only forwards compatibility is maintained. Opening a project file that was saved in a more recent version is
    generally not guarantied to be possible, but opening old project files in newer version should always be possible.

"""
# Internal imports
import graph_2
import core
import utils
# External imports
import logging
# Instantiate logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def convert(obj, old_version):
    """Return project instance in a compatible state.

    :param obj: The core.SuchSoftware instance to be upgraded.
    :param old_version: The version that 'obj' was saved with.

    :type obj: core.SuchSoftware() instance
    :type old_version: list(<int>)

    :returns project instance in a compatible state.
    :rtype core.SuchSoftware() instance

    """

    # Return obj in a compatible state
    # Return None if not possible!

    if old_version == [0, 1, 0]:
        fresh_obj = None

    else:
        try:
            new_project = core.Project('empty')
            new_project.im_mat = obj.im_mat
            (new_project.im_height, new_project.im_width) = new_project.im_mat.shape
            new_project.scale = obj.scale
            new_project.fft_im_mat = utils.gen_fft(new_project.im_mat)
            new_project.filename_full = obj.filename_full
            new_project.alloy = obj.alloy
            new_project.set_alloy_mat()
            new_project.search_mat = obj.search_mat
            new_project.column_centre_mat = obj.column_centre_mat
            new_project.threshold = obj.threshold
            new_project.r = obj.r
            new_project.overhead = obj.overhead
            new_project.num_columns = 0

            new_graph = graph_2.AtomicGraph(new_project.scale, None)
            for vertex in obj.graph.vertices:
                try:
                    new_vertex = graph_2.Vertex(vertex.i, vertex.im_coor_x, vertex.im_coor_y, vertex.r, vertex.scale, zeta=vertex.level, species_index=vertex.species_index, void=vertex.void)
                except:
                    try:
                        new_vertex = graph_2.Vertex(vertex.i, vertex.im_coor_x, vertex.im_coor_y, vertex.r, vertex.scale, zeta=vertex.zeta, species_index=vertex.species_index, void=vertex.void)
                    except:
                        new_vertex = None
                else:
                    if new_vertex is not None:
                        new_vertex.district = vertex.district
                        new_vertex.peak_gamma = vertex.peak_gamma
                        new_vertex.avg_gamma = vertex.avg_gamma
                        new_vertex.normalized_peak_gamma = vertex.normalized_peak_gamma
                        new_vertex.normalized_avg_gamma = vertex.normalized_avg_gamma

                        new_graph.add_vertex(new_vertex)
                        new_project.num_columns += 1
            new_graph.refresh_graph()

            new_project.graph = new_graph

            fresh_obj = new_project
        except:
            fresh_obj = None

    return fresh_obj

