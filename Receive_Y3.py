class Point(object):
    def __init__(self, _x, _y, _order=None):
        self.x, self.y, self.order = _x, _y, _order
 
    def calc(self, top, bottom, other_x):
        l = (top * inverse_mod(bottom)) % p
        x3 = (l * l - self.x - other_x) % p
        return Point(x3, (l * (self.x - x3) - self.y) % p)
 
    def double(self):
        if self == INFINITY: return INFINITY
        return self.calc(3 * self.x * self.x, 2 * self.y, self.x)
 
    def __add__(self, other):
        if other == INFINITY: return self
        if self == INFINITY: return other
        if self.x == other.x:
            if (self.y + other.y) % p == 0: return INFINITY
            return self.double()
        return self.calc(other.y - self.y, other.x - self.x, other.x)
 
    def __mul__(self, e):
        if self.order: e %= self.order
        if e == 0 or self == INFINITY: return INFINITY
        result, q = INFINITY, self
        while e:
            if e & 1: result += q
            e, q = e >> 1, q.double()
        return str(result).split()[0] + ' ' + str(result).split()[1]
        
 
    def __str__(self):
        if self == INFINITY: return "infinity"
        return "= %x %x" % (self.x, self.y)
 
 
def inverse_mod(a):
    if a < 0 or a >= p: a = a % p
    c, d, uc, vc, ud, vd = a, p, 1, 0, 0, 1
    while c:
        q, c, d = divmod(d, c) + (c,)
        uc, vc, ud, vd = ud - q * uc, vd - q * vc, uc, vc
    if ud > 0: return ud
    return ud + p
 
 
p, INFINITY = 115792089237316195423570985008687907853269984665640564039457584007908834671663, Point(None, None)
g = Point(55066263022277343669578718895168534326250603453777594175500187360389116729240,
          32670510020758816978083085130507043184471273380659243275938904335757337482424,
          115792089237316195423570985008687907852837564279074904382605163141518161494337)
Value = 3073472814
 
result = '  X3:    %x\n  Y3:  %s' % (Value, g*Value)
f = open('Value_X3_Y3.txt', 'w')
f.write(result)
f.close()

# Hex to Decimal converter:  https://www.rapidtables.com/convert/number/hex-to-decimal.html
# Python 2.7
