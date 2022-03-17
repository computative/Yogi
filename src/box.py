class Box:

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def __contains__(self, item):
        a,b = item
        condx = self.a[0] < a[0] < self.b[0]
        condy = self.a[1] < a[1] < self.b[1]
        condz = self.a[2] < a[2] < self.b[2]
        if condx and condy and condz:
            return True
        return False