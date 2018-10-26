#! /usr/bin/env python

""" Test inverse log10 factorial function """

import math

#Global constants
LR2P = .5 * math.log(2.*math.pi)  #Log Root 2 Pi
L10 = math.log(10.)               #Log 10

def invlfact(n):
    """ Inverse Stirling's approx for log n factorial, using Newton's method """
    x = y = n * L10 - LR2P
    for i in xrange(3):
        x = (y + x) / math.log(x)
    return int(round(x))


def main():
    N = 30
    f = 6
    for i in xrange(4, N+1):
        f *= i
        x = invlfact(math.log10(f))
        print i, f, x


if __name__ == '__main__':
    main()  
