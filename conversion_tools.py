# Program imports:
import utils
import core
# External imports:
import logging
import h5py
import copy
# Instantiate logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


def make_svg(project, filter_, species_type='atomic_species', graph=None, image=True, radii_factor=0.5):
    """Convert an AutomAl 6000 project file into a Scalable Vector Graphics (svg) file.

    The filter dictionary needs the following fields:



    :param project: An AutomAl 6000 project file.
    :type project: core.Project
    :param filter_: A dictionary with filter settings
    :type filter_: dict
    :param species_type: (Optional, default = 'atomic_species') Keyword indicating grouping by atomic species or advanced species
    :type species_type: str
    :param graph: (Optional, default = None) Keyword indicating which, if any, graph to include
    :type graph: str
    :param image: (Optional, default=True) Keyword indicating if a rectangle should be included in place of the raw image
    :type image: Bool
    :param radii_factor: (Optional, default = 0.5) Scale the atomic radii of the columns with this number
    :type radii_factor: float

    :returns A string that can be saved as an svg file:
    :rtype string:

    """

    xml_string = '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n'
    xml_string += '<!-- Created with AutomAl 6000 (http://www.automal.org/) -->\n\n'

    xml_string += '<svg\n   xmlns:xlink="http://www.w3.org/1999/xlink">\n'

    if graph is not None:

        if graph == 'atomic_graph':

            xml_string += '  <g\n'
            xml_string += '     inkscape:groupmode="layer"\n'
            xml_string += '     inkscape:label="Atomic graph [arcs]"\n'
            xml_string += '     id="Atomic graph [arcs]"\n'
            xml_string += '     style="display:inline" >\n'

            project.graph.map_arcs()

            for arc in project.graph.arcs:

                if arc.dual_arc:

                    xml_string += '    <path\n'
                    xml_string += '       id="arc{}"\n'.format(arc.j)
                    xml_string += '       d="M {},{} {},{}"\n'.format(arc.vertex_a.im_coor_x, arc.vertex_a.im_coor_y, arc.vertex_b.im_coor_x, arc.vertex_b.im_coor_y)

                    if arc.co_planar:

                        xml_string += '       stroke="{}"\n'.format(utils.rgb_to_hex((20, 20, 220)))
                        xml_string += '       stroke-width="4" />\n'

                    else:

                        xml_string += '       stroke="{}"\n'.format(utils.rgb_to_hex((0, 0, 0)))
                        xml_string += '       stroke-width="1" />\n'

                else:

                    xml_string += '    <path\n'
                    xml_string += '       id="arc{}"\n'.format(arc.j)
                    xml_string += '       d="M {},{} {},{}"\n'.format(arc.vertex_a.im_coor_x, arc.vertex_a.im_coor_y, arc.vertex_b.im_coor_x, arc.vertex_b.im_coor_y)
                    xml_string += '       stroke="{}"\n'.format(utils.rgb_to_hex((220, 20, 20)))
                    xml_string += '       stroke-width="4" />\n'

            xml_string += '  </g>\n'

            xml_string += '  <g\n'
            xml_string += '     inkscape:groupmode="layer"\n'
            xml_string += '     inkscape:label="Atomic graph [vertices]"\n'
            xml_string += '     id="Atomic graph [vertices]"\n'
            xml_string += '     style="display:inline" >\n'

            for vertex in project.graph.vertices:

                xml_string += '    <circle\n'
                xml_string += '       id="graph_vertex_{}"\n'.format(vertex.i)
                xml_string += '       cx="{}"\n'.format(vertex.im_coor_x)
                xml_string += '       cy="{}"\n'.format(vertex.im_coor_y)
                xml_string += '       r="{}"\n'.format(vertex.r / 2 - 8 / project.scale)

                if vertex.zeta == 0:
                    xml_string += '       style="fill:{};fill-opacity:1;stroke:{};stroke-width:{};stroke-opacity:1"\n'.format(
                        utils.rgb_to_hex((0, 0, 0)),
                        utils.rgb_to_hex((0, 0, 0)),
                        15 / project.scale
                    )
                else:
                    xml_string += '       style="fill:{};fill-opacity:1;stroke:{};stroke-width:{};stroke-opacity:1"\n'.format(
                        utils.rgb_to_hex((255, 255, 255)),
                        utils.rgb_to_hex((0, 0, 0)),
                        15 / project.scale
                    )

            xml_string += '  </g>\n'

    if image:

        # im = Image.fromarray(project.im_mat).tobytes()
        # im_b_string = base64.b64encode(im)
        # im_b_string = '{}'.format(im_b_string)
        # im_b_string = im_b_string[2:]

        xml_string += '  <g\n'
        xml_string += '     inkscape:groupmode="layer"\n'
        xml_string += '     inkscape:label="Image"\n'
        xml_string += '     id="Image"\n'
        xml_string += '     style="display:inline" >\n'

        # xml_string += '    <image\n'
        # xml_string += '       x="0"\n'
        # xml_string += '       y="0"\n'
        # xml_string += '       width="{}"\n'.format(project.im_width)
        # xml_string += '       height="{}"\n'.format(project.im_height)
        # xml_string += '       xlink:href="data:image/png;base64,{}" />\n'.format(im_b_string)

        xml_string += '    <rect\n'
        xml_string += '       id="image_frame"\n'
        xml_string += '       x="0"\n'
        xml_string += '       y="0"\n'
        xml_string += '       width="{}"\n'.format(project.im_width)
        xml_string += '       height="{}"\n'.format(project.im_height)
        xml_string += '       fill="none"\n'
        xml_string += '       stroke="#000000"\n'
        xml_string += '       stroke-width="1" />\n'

        xml_string += '  </g>\n'

    xml_string += '  <g\n'
    xml_string += '     inkscape:groupmode="layer"\n'
    xml_string += '     inkscape:label="matrix"\n'
    xml_string += '     id="matrix"\n'
    xml_string += '     style="display:inline" >\n'

    if not filter_['exclude_matrix_columns']:
        if species_type == 'advanced_species':
            color = project.species_dict[species_type]['Al_1']['color']
        else:
            color = project.species_dict[species_type]['Al']['color']
        pixel_radii = radii_factor * (project.species_dict['atomic_species']['Al']['atomic_radii'] - 15) / project.scale
        for vertex in project.graph.vertices:
            if not vertex.is_in_precipitate:
                if vertex.zeta == 0:
                    in_color = color
                else:
                    in_color = (0, 0, 0)
                xml_string += '    <circle\n'
                xml_string += '       id="vertex{}"\n'.format(vertex.i)
                xml_string += '       style="fill:{};fill-opacity:1;stroke:{};stroke-width:{};stroke-opacity:1"\n'.format(
                    utils.rgb_to_hex(in_color),
                    utils.rgb_to_hex(color),
                    30 / project.scale
                )
                xml_string += '       cx="{}"\n'.format(vertex.im_coor_x)
                xml_string += '       cy="{}"\n'.format(vertex.im_coor_y)
                xml_string += '       r="{}" />\n'.format(pixel_radii)

    xml_string += '  </g>\n'

    if species_type == 'advanced_species':
        species_population = project.graph.get_advanced_species_list()
    else:
        species_population = project.graph.get_atomic_species_list()

    for species in species_population:

        xml_string += '  <g\n'
        xml_string += '     inkscape:groupmode="layer"\n'
        xml_string += '     inkscape:label="{}"\n'.format(species)
        xml_string += '     id="{}"\n'.format(species)
        xml_string += '     style="display:inline" >\n'

        if not filter_['exclude_particle_columns']:
            color = project.species_dict[species_type][species]['color']
            if species_type == 'advanced_species':
                pixel_radii = radii_factor * (project.species_dict['atomic_species'][project.species_dict[species_type][species]['atomic_species']]['atomic_radii'] - 15) / project.scale
            else:
                pixel_radii = radii_factor * (project.species_dict['atomic_species'][species]['atomic_radii'] - 15) / project.scale
            for vertex in project.graph.vertices:
                if vertex.is_in_precipitate and eval('vertex.{}'.format(species_type)) == species:
                    if vertex.zeta == 0:
                        in_color = color
                    else:
                        in_color = (0, 0, 0)
                    xml_string += '    <circle\n'
                    xml_string += '       id="vertex{}"\n'.format(vertex.i)
                    xml_string += '       style="fill:{};fill-opacity:1;stroke:{};stroke-width:{};stroke-opacity:1"\n'.format(
                        utils.rgb_to_hex(in_color),
                        utils.rgb_to_hex(color),
                        30 / project.scale
                    )
                    xml_string += '       cx="{}"\n'.format(vertex.im_coor_x)
                    xml_string += '       cy="{}"\n'.format(vertex.im_coor_y)
                    xml_string += '       r="{}" />\n'.format(pixel_radii)

        xml_string += '  </g>\n'

    xml_string += '</svg>\n'

    return xml_string


def import_from_atomap(filename, scale, debug_obj=None):

    column_file = h5py.File(filename, 'r+')

    project = core.Project('empty', debug_obj=debug_obj, species_dict=core.Project.default_species_dict)
    project.filename_full = filename

    for group_name in column_file:
        if ('atom_lattice' in group_name) or ('sublattice' in group_name):

            data = column_file[group_name]
            original_image_data = data['original_image_data'][:]
            atom_position_array = data['atom_positions'][:]

            project.im_mat = original_image_data
            project.scale = scale
            if not scale == 1000 * data.attrs['pixel_size']:
                logger.warning('The indicated scale ({}), does not agree with the scale indication in the file ({}).'.format(
                    scale,
                    data.attrs['pixel_size']
                ))
                if not data.attrs['pixel_size'] == 1.0:
                    project.scale = 1000 * data.attrs['pixel_size']
            (project.im_height, project.im_width) = project.im_mat.shape
            project.im_mat = utils.normalize_static(project.im_mat)
            project.fft_im_mat = utils.gen_fft(project.im_mat)
            project.search_mat = copy.deepcopy(project.im_mat)
            project.calc_avg_pixel_value()
            project.r = int(100 / project.scale)
            project.overhead = int(7 * (project.r / 10))
            project.graph = core.graph_2.AtomicGraph(
                project.scale,
                active_model=None,
                species_dict=project.species_dict,
                district_size=project.district_size
            )

            for i, column in enumerate(atom_position_array):
                new_vertex = core.graph_2.Vertex(
                    i,
                    column[0],
                    column[1],
                    project.r,
                    project.scale,
                    parent_graph=project.graph
                )
                new_vertex.avg_gamma, new_vertex.peak_gamma = utils.circular_average(
                    project.im_mat,
                    int(column[0]),
                    int(column[1]),
                    project.r
                )
                project.graph.add_vertex(new_vertex)
                project.num_columns += 1

            project.find_edge_columns()

            logger.info('Generated instance from {}'.format(project.filename_full))

    return project





