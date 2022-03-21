import numpy as np

class c(np.ndarray):
    def __new__(cls, *args):
        return np.array(args)
