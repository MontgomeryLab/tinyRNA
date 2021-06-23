import functools
import time


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
