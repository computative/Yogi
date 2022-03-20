import numpy as np

class demo1:
    _instance = None
    a = 3

    def __new__(cls):
        if cls._instance == None:
            cls._instance = super().__new__(cls)\
        super().__init__()
        return cls._instance

    def __init__(self):
        self.hello = "hello"

class demo2:
    _instance = None
    a = 3

    def __new__(cls):
        if cls._instance == None:
            cls._instance = super().__new__(cls)
            return np.arange(3)
        return np.arange(4,7)


D = demo1()
E = demo1()
F = demo1()
F = 5


a=3
