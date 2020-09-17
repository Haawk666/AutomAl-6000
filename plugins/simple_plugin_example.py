import logging

# Instantiate logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Bridge:

    def __init__(self, ui_obj):
        self.ui_obj = ui_obj

    def trigger(self):
        if self.ui_obj.project_instance is not None:
            logger.info('Num columns: {}'.format(self.ui_obj.project_instance.num_columns))
        else:
            logger.info('Not instantiated yet!')


