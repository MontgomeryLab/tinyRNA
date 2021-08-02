import functools
import time
import os


def report_execution_time(step_name: str):
    def timer(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            start = time.time()
            return_val = func(*args, **kwargs)
            end = time.time()
            print("%s took %.2f seconds" % (step_name, end - start))
            return return_val
        return wrapper
    return timer


def from_here(config_file, input_file):
    """Calculates paths relative to the config file which contains them"""
    config_file, input_file = (os.path.expanduser(p) for p in [config_file, input_file])

    if not os.path.isabs(input_file):
        from_here = os.path.dirname(config_file)
        input_file = os.path.normpath(os.path.join(from_here, input_file))

    return input_file
