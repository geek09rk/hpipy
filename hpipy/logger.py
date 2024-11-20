import os
import logging

def logHPIpy(logdir = None, mode = None):
    """
    Initialize a logger for HPIpy with specified directory and mode.

    :param logdir: Directory where the log file will be saved
    :param mode: Mode in which the log file is opened. Defaults to None. If not provided, the default mode is 'a' (append).
    
    :type logdir: str
    :type mode: str

    :return: Logger object configured with the specified settings.
    :rtype: str
    """

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    file_handler = logging.FileHandler(os.path.join(logdir, "HPIpy.log"), mode = mode)
    formatter = logging.Formatter('%(asctime)s :: %(levelname)s :: %(message)s', datefmt = '%b-%d-%Y %H:%M:%S')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    return logger