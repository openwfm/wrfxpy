class a(object):
    def __init__(self):
        self.v=0
        print(self.j)
    def m(self):
        print('a.m')
class b(a):
    def n(self):
        print('b.n')
class c(b):
    def m(self):
        print('c.m')
    def n(self):
        print('c.n')
i=[0,1]
j=i[-1]
      
