# import sys
# print(sys._getframe().f_back.f_code.co_name)

class measured_state():
    def __init__(self):
        self.m = {}

    def __repr__(self, key):
        return self.m[0]

    def update(self, **kwargs):
        self.m = {}
        for arg, value in kwargs.items():
            self.m.update({arg:value})


m = measured_state()

m.update(a=1, b=2, c=3)

print(m.m['a'])
print(m.m['b'])
print(m.m['c'])

m.update(a=3, b=4, c=5)

print(m.m['a'])
print(m.m['b'])
print(m.m['c'])