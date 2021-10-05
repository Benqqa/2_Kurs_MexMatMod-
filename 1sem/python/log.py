import logging
import functools


def _generate_log(path):

    logger = logging.getLogger('LogError')
    logger.setLevel(logging.ERROR)
    file_handler = logging.FileHandler(path)
    log_format = 'Error time: %(asctime)s %(message)s'
    formatter = logging.Formatter(log_format)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    return logger


def log_error(path='Log.log'):
    def error_log(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                logger = _generate_log(path)
                error_msg = 'and it has occurred at ' + func.__name__ + '\n-------------------------\n'
                logger.exception(error_msg)
                return e
        return wrapper
    return error_log
