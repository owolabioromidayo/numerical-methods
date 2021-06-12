#!/usr/bin/env python3

import numpy as np

def gaussElim(a,b):
    n = len(b)
    #elimination phase
    for k in range(0,n-1):
        for i in range(k+1, n):
            if a[i,k] != 0.0:
                _lambda = a[i,k] / a[k,k]
                a[i,k+1:n] = a[i,k+1:n] - _lambda*a[k,k+1:n]
                b[i] = b[i] - _lambda*b[k]

    #back-substitution phase
    for k in range(n-1, -1, -1):
        b[k] = (b[k] - np.dot(a[k,k+1:n], b[k+1:n]))/ a[k,k]
    return b

if __name__ == "__main__":
    a = np.array([[4,-2,1],[-2,4,-2],[1,-2,4]])
    b = np.array([11,-16,17])
    print(gaussElim(a,b))
