__author__ = 'Chaoren'
import copy

class test:
    def __init__(self, a, b):
        self.a = a
        self.b = b
    def __add__(self, other):
        t = copy.copy(self)
        t.a = self.a - other.a
        t.b = self.b - other.b
        return t

o1 = test(10,2)
o2 = test(1,4)

print (o1+o2).a
