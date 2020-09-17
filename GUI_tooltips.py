# By Haakon Tvedt @ NTNU
# Contributors:
"""Module container for tooltip strings"""

tooltips = {
    'project': {
        'lbl_filename': 'Project file.',
        'lbl_location': 'Project location.',
        'sbtn_active_model':
            'The currently active model. Changing the active model will influence the results of\n'
            'the column characterization.',
        'btn_import': 'Import image.',
        'btn_open_project': 'Open project.',
        'btn_save_project': 'Save current project.',
        'btn_close_project': 'Close the current project.',
        'btn_species_dict': 'Edit the species dictionary of the project.',
        'btn_export': 'Export project data with the export wizard.',
        'btn_show_stats': 'Print project statistics to the terminal window.'
    }, 'image': {
        'lbl_image_width': 'Image width in pixels.',
        'lbl_image_height': 'Image height in pixels.',
        'lbl_upsampled': 'Indicates whether the image has been upsampled with the bilinear method.',
        'chb_lock_views': 'Lock the different tabs to always display the same area (pluss/minus 2 pixels).',
        'btn_show_source': 'Print the path of the original image file to the terminal window.',
        'btn_align_views': 'Align all tabs with the current area of the current tab. [E]'
    }, 'column_detection': {
        'lbl_atomic_radii':
            'The \'approximate atomic radii, r\' in pixels. Provides a rough estimate of atomic radii, which is used\n'
            'during column detection. This value is calculated from the scale of the image.',
        'lbl_overhead':
            'The \'overhead\' parameter is a parameter used in column detection. This parameter is\n'
            'calculated from the scale of the image.',
        'lbl_pixel_average': 'The average pixel intensity of the image.',
        'lbl_search_matrix_peak': 'The currently brightest pixel in the search matrix.',
        'lbl_num_detected_columns': 'The number of columns currently detected in the image.',
        'chb_plot_results': 'If checked, a summary of the relevant values will be printed after column detection.',
        'chb_automatic_threshold': 'Not in use.',
        'chb_toggle_positions':
            'Toggle whether to display position overlay. Can be used to toggle between the raw\n'
            'image and the position overlay.',
        'sbtn_threshold':
            'Set the threshold value (between 0 and 1). If the maximum pixel value of the search\n'
            'matrix is below this value, column detection will stop searching. A good starting\n'
            'value to try for most images is for instance 0,3.',
        'sbtn_search_size':
            'Set the maximum number of columns that column detection will try to find.',
        'sbtn_scale':
            'Set the scale of the image in (pm/pixel). This should be automatically set from the .dm3\n'
            'metadata, but can be overridden here if needed.',
        'btn_start_column_detection': 'Begin column detection with current settings.',
        'btn_reset_column_detection': 'Delete all columns.',
        'btn_redraw_search_mat':
            'Redraw the search matrix. This must be done if columns have been moved manually,\n'
            'and the column detection is to be continued.'
    }, 'column_characterization': {
        'chb_show_graphics': 'Show graphical updates while the algorithm is running. (This will be slightly slower).',
        'lbl_alloy': 'The elements present according to the current species dictionary of the image.',
        'btn_start_column_characterization': 'Start column characterization.',
        'btn_reset_column_characterization': 'Reset column characterization',
        'btn_invert_zeta': 'Invert all zeta values.'
    }, 'selected_column': {
        'lbl_im_coor_x': 'The x-coordinate of the selected column in image coordinates.',
        'lbl_im_coor_y': 'The y-coordinate of the selected column in image coordinates.',
        'lbl_spatial_coor_x': 'The x-coordinate of the selected column in spatial coordinates.',
        'lbl_spatial_coor_y': 'The y-coordinate of the selected column in spatial coordinates.',
        'lbl_peak_gamma': 'The brightest pixel intensity of the selected column.',
        'lbl_avg_gamma': 'The average pixel intensity of the selected column.',
        'lbl_area_gamma': 'The area weighted average pixel intensity of the selected column.',
        'lbl_normalized_peak_gamma':
            'The brightest pixel intensity of the selected column, normalized against the matrix.',
        'lbl_normalized_avg_gamma':
            'The average pixel intensity of the selected column, normalized against the matrix.',
        'lbl_normalized_area_gamma':
            'The area weighted average pixel intensity of the selected column, normalized against the matrix.',
        'lbl_alpha_min': 'Minimum alpha angle of the selected column.',
        'lbl_alpha_max': 'Maximum alpha angle of the selected column.',
        'lbl_theta_min': 'Minimum theta angle of the selected column.',
        'lbl_theta_max': 'Maximum theta angle of the selected column.',
        'lbl_theta_mean': 'The mean theta angle of the selected column.',
        'lbl_theta_angle_variance': 'The theta variance of the selected column.',
        'chb_precipitate_column': 'If checked, the selected column is considered to be within the precipitate.',
        'chb_show_in_overlay': 'If un-checked, the selected column will not be shown in the atomic overlay.',
        'mchb_move': 'Enable manual repositioning of the selected column.',
        'sbtn_find_column': 'Unique column index.',
        'sbtn_atomic_species': 'The \'atomic species\' of the selected column.',
        'sbtn_advanced_species': 'The \'advanced species\' of the selected column.',
        'sbtn_zeta': 'The zeta value of the selected column.',
        'btn_delete_column': 'Delete the selected column.',
        'btn_print_details': 'Print additional details of the selected column to the terminal window.',
        'btn_snap': 'Zoom in all tabs on the currently selected column.',
        'btn_deselect': 'Deselect the column.',
        'btn_new_column':
            'Insert a new column. (This produces a new column in the middle of the image and enables move-mode.\n'
            'Move the column into position manually and hit enter to accept the position.)'
    }, 'atomic_graph': {
        'lbl_order': 'The order of a graph is the number of vertices in the graph.',
        'lbl_size': 'The size of a graph is the number of arcs in the graph.',
        'lbl_volume': 'In this application, the volume of a graph is the number of meshes in the graph.',
        'lbl_chi':
            'Chi is the number of un-symmetric arcs divided by the total number of arcs (the graph size).\n' 
            'Chi gives a indication of the symmetry density of the graph (chi=0 => Graph is symmetric,\n' 
            'chi=1 => is fully un-symmetric).',
        'lbl_zeta_chi': 'The chi of the zeta graph.',
        'lbl_average_degree': 'The average degree of the graph.',
        'lbl_average_zeta_degree': 'The average degree of the zeta graph.',
        'chb_permute_mode':
            'If checked, selecting three vertices (v_a, v_b, v_c) in sequence while in the atomic graph tab, will\n'
            'perform the permutation P_a(b, c) if possible.',
        'chb_enable_ruler':
            'If checked, selecting two vertices in sequence, will print the distance between the columns to the\n'
            'terminal window.',
        'chb_show_mesh_vertex_order':
            'Show the mesh orders in the atomic graph. Meshes that has order other than 4, will stand out in red print',
        'chb_show_zeta_0': 'Toggle vertices with zeta=0 in the anti-graph.',
        'chb_show_zeta_1': 'Toggle vertices with zeta=1 in the anti-graph.',
        'btn_refresh_graph': 'Refresh graph properties.',
        'btn_refresh_mesh': 'Refresh the mesh list of the atomic graph.',
        'btn_print_hsd':
            'For reference, print all the hard-sphere distances (HSD) between the elements present in the species\n'
            'dictionary.',
        'btn_gen_sub': 'Generate subgraphs with the sub-graph wizard.',
        'btn_build_anti_graph': 'Build the anti-graph of the atomic graph.',
        'btn_match_zeta': 'Permute the arcs of the atomic graph so as to match the zeta graph.'
    }, 'heat_maps': {
        'sbtn_current_map': 'Set the heat map display.'
    }, 'data': {
        'btn_project_plots': 'Produce data plots from current project.',
        'btn_model_plots': 'Produce data plots from available models.',
        'btn_calc_models': 'Calculate a new model with the model wizard.',
        'btn_calc_heat': 'Calculate heat-maps.'
    }, 'overlay': {
        'btn_configure': 'Configure overlay settings.',
        'btn_show_all': 'Show all columns in the overlay.',
        'btn_hide_all': 'Hide all columns from the overlay.',
        'chb_raw_image': 'Toggle HAADF-STEM image.',
        'chb_black_background': 'Toggle black background',
        'chb_edge_columns':
            'Toggle columns that are on the edge of the image (These are not analysed by the column characterization).',
        'chb_precipitate': 'Toggle precipitate columns',
        'chb_matrix': 'Toggle Al matrix columns.',
        'chb_legend': 'Toggle movable legend',
        'chb_scalebar': 'Toggle movable scalebar',
        'chb_zeta_0': 'Toggle columns with zeta=0',
        'chb_zeta_1': 'Toggle columns with zeta=1'
    }, 'debug': {
        'sbtn_starting_index':
            'Set a default starting column for the column characterization. Used to study deterministic\n'
            'results when testing the algorithms',
        'btn_test': '',
        'btn_crash': '',
        'btn_test_data_manager': '',
        'btn_show_prediction': ''
    }
}

