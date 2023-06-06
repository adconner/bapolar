from sage.all import PolynomialRing
from itertools import product

def matrixmult(m,n,l):
    dims = (m*n,n*l,l*m)
    R = PolynomialRing(QQ,['u%d%d' % (i,j) for i,j in product(range(m),range(n))]
                      +['v%d%d' % (i,j) for i,j in product(range(n),range(l))]
                      +['w%d%d' % (i,j) for i,j in product(range(l),range(m))])
    u = R.gens()[:m*n]
    v = R.gens()[m*n:m*n+n*l]
    w = R.gens()[m*n+n*l:]
    T = sum(u[i*n+j]*v[j*l+k]*w[k*m+i] for i,j,k in product(range(m),range(n),range(l)))
    return (T,dims)
    
