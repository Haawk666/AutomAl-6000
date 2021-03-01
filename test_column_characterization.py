import core
import copy
import time
import logging
# Instantiate logger:
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class TestColumnCharacterizationAlgorithm:

    def __init__(self, column_ready_instances, control_instances):
        logger.info('Testing algorithms...')
        results = []
        for project, control in zip(column_ready_instances, control_instances):
            result = []
            control = core.Project.load(control)
            project = core.Project.load(project)

            thesis_git_result = copy.deepcopy(project)
            time_1 = time.time()
            thesis_git_result = self.thesis_git_characterization(thesis_git_result)
            time_2 = time.time()
            num_error = 0
            num_precipitate = 0
            for vertex in thesis_git_result.graph.vertices:
                if vertex.is_in_precipitate:
                    num_precipitate += 1
                    if not vertex.atomic_species == control.graph.vertices[vertex.i].atomic_species:
                        num_error += 1
            percent_error = 100 * (num_error / num_precipitate)
            result.append([thesis_git_result.graph.chi, time_2 - time_1, percent_error])

            thesis_result = copy.deepcopy(project)
            time_1 = time.time()
            thesis_result = self.thesis_characterization(thesis_result)
            time_2 = time.time()
            num_error = 0
            num_precipitate = 0
            for vertex in thesis_result.graph.vertices:
                if vertex.is_in_precipitate:
                    num_precipitate += 1
                    if not vertex.atomic_species == control.graph.vertices[vertex.i].atomic_species:
                        num_error += 1
            percent_error = 100 * (num_error / num_precipitate)
            result.append([thesis_result.graph.chi, time_2 - time_1, percent_error])

            pre_paper_result = copy.deepcopy(project)
            time_1 = time.time()
            pre_paper_result = self.pre_paper_characterization(pre_paper_result)
            time_2 = time.time()
            num_error = 0
            num_precipitate = 0
            for vertex in pre_paper_result.graph.vertices:
                if vertex.is_in_precipitate:
                    num_precipitate += 1
                    if not vertex.atomic_species == control.graph.vertices[vertex.i].atomic_species:
                        num_error += 1
            percent_error = 100 * (num_error / num_precipitate)
            result.append([pre_paper_result.graph.chi, time_2 - time_1, percent_error])

            paper_result = copy.deepcopy(project)
            time_1 = time.time()
            paper_result = self.paper_characterization(paper_result)
            time_2 = time.time()
            num_error = 0
            num_precipitate = 0
            for vertex in paper_result.graph.vertices:
                if vertex.is_in_precipitate:
                    num_precipitate += 1
                    if not vertex.atomic_species == control.graph.vertices[vertex.i].atomic_species:
                        num_error += 1
            percent_error = 100 * (num_error / num_precipitate)
            result.append([paper_result.graph.chi, time_2 - time_1, percent_error])

            new_result = copy.deepcopy(project)
            time_1 = time.time()
            new_result = self.new_characterization(new_result)
            time_2 = time.time()
            num_error = 0
            num_precipitate = 0
            for vertex in new_result.graph.vertices:
                if vertex.is_in_precipitate:
                    num_precipitate += 1
                    if not vertex.atomic_species == control.graph.vertices[vertex.i].atomic_species:
                        num_error += 1
            percent_error = 100 * (num_error / num_precipitate)
            result.append([new_result.graph.chi, time_2 - time_1, percent_error])

            results.append(result)

        logger.info('testing complete.')
        string = 'Results:\n'
        for index, result in enumerate(results):
            string += '    Project {}:\n'.format(index)
            for version in [0, 1, 2, 3, 4]:
                string += '        v{}: {:.4f}, {:.4f}, {:.1f}\n'.format(version, result[version][0], result[version][1], result[version][2])
        logger.info(string)
        print(string)

    @staticmethod
    def thesis_git_characterization(project):
        logger.info('Running thesis version...')
        sub_methods = [4, 5, 6, 11, 7, 6, 11, 9, 10, 12, 13, 11, 9, 10, 12, 13, 11, 8, 20]
        if project.starting_index is not None:
            starting_index = project.starting_index
        else:
            starting_index = 0
        for sub_method in sub_methods:
            project.column_characterization_2(starting_index, sub_method=sub_method, indent='    ')
        logger.info('Thesis (git) characterization complete.')
        return project

    @staticmethod
    def thesis_characterization(project):
        logger.info('Running thesis version...')
        sub_methods = [4, 5, 6, 11, 7, 6, 11, 9, 10, 12, 13, 11, 8, 6, 9, 10, 12, 13, 11, 20]
        if project.starting_index is not None:
            starting_index = project.starting_index
        else:
            starting_index = 0
        for sub_method in sub_methods:
            project.column_characterization_2(starting_index, sub_method=sub_method, indent='    ')
        logger.info('Thesis characterization complete.')
        return project

    @staticmethod
    def pre_paper_characterization(project):
        logger.info('Running thesis version...')
        sub_methods = [4, 5, 6, 7, 6, 9, 10, 11, 12, 13, 8, 9, 10, 11, 12, 13, 11, 20]
        if project.starting_index is not None:
            starting_index = project.starting_index
        else:
            starting_index = 0
        for sub_method in sub_methods:
            project.column_characterization_2(starting_index, sub_method=sub_method, indent='    ')
        logger.info('Pre-paper characterization complete.')
        return project

    @staticmethod
    def paper_characterization(project):
        logger.info('Running thesis version...')
        sub_methods = [4, 5, 6, 7, 9, 10, 11, 12, 13, 8, 9, 10, 11, 12, 13, 8, 9, 10, 11, 12, 13, 8, 20]
        if project.starting_index is not None:
            starting_index = project.starting_index
        else:
            starting_index = 0
        for sub_method in sub_methods:
            project.column_characterization_2(starting_index, sub_method=sub_method, indent='    ')
        logger.info('Paper characterization complete.')
        return project

    @staticmethod
    def new_characterization(project):
        logger.info('Running thesis version...')
        sub_methods = [4, 5, 6, 11, 7, 6, 11, 9, 10, 11, 12, 13, 8, 9, 10, 11, 12, 13, 8, 9, 10, 11, 12, 13, 8, 17, 9, 10, 11, 20]
        if project.starting_index is not None:
            starting_index = project.starting_index
        else:
            starting_index = 0
        for sub_method in sub_methods:
            project.column_characterization_2(starting_index, sub_method=sub_method, indent='    ')
        logger.info('New characterization complete.')
        return project

