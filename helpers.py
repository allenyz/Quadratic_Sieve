import numpy as np
import math

#creates factoring interval
#INPUT: n, interval size
#OUTPUT: array of numbers using (sqrt(n) + i)^2 - n

def interval(n, i, start = 0):
    a = []
    for x in range(start, i):
        a.append(x_to_number(n, x))

    return a

def x_to_number(n, x):
    root = math.ceil(n**.5)
    return (root + x)**2 - n

def gcd(a,b): # Euclid's algorithm
    if b == 0:
        return a
    elif a >= b:
        return gcd(b,a % b)
    else:
        return gcd(b,a)

def isqrt(n): # Newton's method, returns exact int for large squares
    x = n
    y = (x + 1) // 2
    while y < x:
        x = y
        y = (x + n // x) // 2
    return x

