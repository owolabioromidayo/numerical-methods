#!/usr/bin/env python3

from math import exp 

def trapezoidal(f,a,b,n):
    "implements the composite trapezoidal rule for numerical integration"
    h = (b-a)/n
    result  = 0.5*(f(a) + f(b))
    for i in range(1,n):
        result += f(a + i*h)
    result *= h
    return result

def midpoint(f,a,b,n):
    "implements the composite midpoint method for numerical integration"
    result = 0
    w = (b-a)/n #width of each bin
    for i in range(n):
        result += f(a + w/2 + i*w)
    result *= w
    return result


if __name__ == "__main__":
    v = lambda t: 3*(t**2)*exp(t**3)
    n = 4
    numerical_approx = trapezoidal(v,0,1,n)
    print(numerical_approx)
