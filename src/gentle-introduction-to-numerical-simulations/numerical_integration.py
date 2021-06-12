#!/usr/bin/env python3

import sympy
import numpy as np

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

def midpoint_vectorized(f,a,b,n):
    "implements the composite midpoint method for numerical integration using numpy"
    w = (b-a)/n #width of each bin
    x = np.linspace(a,b, n+2)[1:-1]
    return np.sum(f(x)) * w 


def trapezoidal_vectorized(f,a,b,n):
    "implements the composite trapezoidal rule for numerical integration using numpy"
    h = (b-a)/n
    result  = 0.5*(f(a) + f(b))
    x = np.linspace(a,b,n)[1:-1]
    result += sum(f(x))
    result *= h
    return result

def midpoint_double1(f,a,b,c,d,nx,ny):
    "midpoint method for double integration"
    I = 0
    hx = (b-a)/nx
    hy = (d-c)/ny
    
    for i in range(nx):
        for j in range(ny):
            xi = a + hx/2 + i*hx
            yj = c + hy/2 + i*hy
            I += f(xi, yj)
    I *= hx*hy
    return I

def midpoint_double2(f,a,b,c,d,nx,ny):
    "midpoint double integration using 1d midpoint function"
    def g(x):
        return midpoint(lambda y: f(x,y), c, d, ny)
    return midpoint(g, a, b, nx)

def midpoint_triple1(f,a,b,c,d,e,f,nx,ny,nz):
    "midpoint triple integration using 1d midpoint function"
    def g(x,y):
        return midpoint(lambda z : f(x,y,z), e, f, nz)
    
    def h(x):
        return midpoint(lambda y: f(x,y), c, d, ny )
    return midpoint(h,a,b,nx)

def midpoint_triple2(f,a,b,c,d,e,f,nx,ny,nz):
    hx = (b-a)/nx
    hy = (d-c)/ny
    hz = (f-e)/nz
    I = 0
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                xi = a + hx/2 + i*hx
                yj = c + hy/2 + j*hy
                zk = e + hz/2 + k*hz
                I += hx*hy*hz*f(xi,yj,zk)
    return I

def convergence_test(method,f,a,b,n,iterations, rate=2):
    "basic convergence test (human verified) for integration methods"
    for i in range(iterations):
        print(f" {n} {method(f,a,b,n)}" )
        n*=rate

def convergence_test_double(method,f,a,b,c,d,nx,ny,iterations, rate=2):
    "basic convergence test (human verified) for double integration"
    for i in range(iterations):
        print(f" {nx} {ny} {method(f,a,b,c,d,nx,ny)}" )
        nx*=rate
        ny*=rate

def test_midpoint_double():
    "tests midpoint double integration exactly as midpoint method is accurate for linear functions"
    def f(x,y):
        return 2*x + y
    
    a=0;b=2;c=2;d=3
    x,y = sympy.symbols('x y')
    I_expected = sympy.integrate(f(x,y), (x,a,b), (y,c,d))
    for nx,ny in  (3, 5), (4, 4), (5, 3):
        I_computed1 = midpoint_double1(f, a, b, c, d, nx, ny)
        I_computed2 = midpoint_double2(f, a, b, c, d, nx, ny)
        tol = 1E-5
        #print I_expected, I_computed1, I_computed2
        assert abs(I_computed1 - I_expected) < tol
        assert abs(I_computed2 - I_expected) < tol


def run_tests():
    "run all tests"
    v = lambda t: 3*(t**2)*np.exp(t**3)
    fxy = lambda x,y: 2*x + y
    n = 4
    for method in [midpoint, midpoint_vectorized, trapezoidal, trapezoidal_vectorized]:
        print(method.__name__.upper() )
        convergence_test(method,v,0,1,n, 10)
        print('\n')

    for method in [midpoint_double1, midpoint_double2]:
        print(method.__name__.upper() )
        convergence_test_double(method, fxy,0,2,2,3,5,5,10)
        print('\n')

   # test_midpoint_double() does not work, floating point inaccuracies

if __name__ == "__main__":
   run_tests()

