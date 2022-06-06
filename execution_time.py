import functools
import time
import logging
from inspect import getmembers, isfunction, ismethod
import sys


# Built by Siddhant Chhabra(@siddhant-curious) under MIT License

class ExecutionTime:
    logFormatter = logging.Formatter("%(asctime)s [%(levelname)-5.5s]  %(message)s")
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.DEBUG)

    def __init__(self, console=False, module_name=None):
        self.console = console
        self.module_name = module_name
        self.logtime_data = {}

        if self.console:
            self.enable_console()

        if self.module_name is not None:
            self.auto_decorate()

    def timeit(self, method):
        @functools.wraps(method)
        def wrapper(*args, **kwargs):
            start_time = time.perf_counter()
            result = method(*args, **kwargs)
            end_time = time.perf_counter()
            total_time = round((end_time - start_time) * 1000, 4)  # time in milliseconds

            if method.__name__ in self.logtime_data:
                curr = self.logtime_data[method.__name__]
                tt = curr["total_time"] + total_time
                count = curr["times_called"] + 1
                avg_time = round(tt / count, 4)
                self.logtime_data[method.__name__] = {'times_called': count, "total_time": tt, "average_time": avg_time}
            else:
                self.logtime_data[method.__name__] = {'times_called': 1, "total_time": total_time,
                                                      "average_time": total_time}

            if self.console is True:
                ExecutionTime.rootLogger.info(f'Time take by method : {method.__name__} is {total_time} ms')
            return result

        return wrapper

    def enable_console(self):
        consoleHandler = logging.StreamHandler(sys.stdout)
        consoleHandler.setFormatter(ExecutionTime.logFormatter)
        ExecutionTime.rootLogger.addHandler(consoleHandler)

    def auto_decorate(self):
        try:
            module = sys.modules[self.module_name]
            items = getmembers(module, isfunction)
            for name, addr in items:
                setattr(module, name, self.timeit(addr))
        except KeyError as e:
            raise Exception(f'Error Occured, No module by name {self.module_name}')
