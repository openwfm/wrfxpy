class a(object):
    def __init__(self):
        self.v=0
    def m(self):
        print 'a.m'
    def x(self):
        print 'a.x'
class b(a):
    def n(self):
        print 'b.n'
class c(b):
    def m(self):
        print 'c.m'
    def n(self):
        print 'c.n'
      
