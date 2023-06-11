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
    
def skewcw2():
    dims=(3,3,3)
    R = PolynomialRing(QQ,['%s%d'%(x,i) for x in 'abc' for i in range(3)])
    T = sum([sigma.sign()*prod(R.gen(i*3+sigma(i+1)-1) for i in range(3))  for sigma in SymmetricGroup(3)])
    return T,dims

def tensor_kronecker_product(Sdat,Tdat):
    S,Sdim = Sdat
    T,Tdim = Tdat
    F = S.base_ring()
    assert F == T.base_ring()
    assert len(Sdim) == len(Tdim)
    R = PolynomialRing(F,['%s%d' % (x,i) 
                      for x,a,b in zip('abcdefghijklmnopqrstuvwxyz',Sdim,Tdim) 
                      for i in range(a*b)])
    ST = R.zero()
    for (e1,c1),(e2,c2) in product(S.dict().items(),T.dict().items()):
        term = c1*c2
        yi = 0
        for i,((x1,k1),d1,(x2,k2),d2) in enumerate(zip(e1.sparse_iter(),Sdim,e2.sparse_iter(),Tdim)):
            x1 -= sum(Sdim[:i])
            x2 -= sum(Tdim[:i])
            assert k1 == 1
            assert k2 == 1
            assert x1 <= d1
            assert x2 <= d2
            term *= R.gen(yi+x1*d2+x2)
            yi += d1*d2
        ST += term
    return ST,tuple(d1*d2 for d1,d2 in zip(Sdim,Tdim))

def tensor_sum(Sdat,Tdat):
    S,Sdim = Sdat
    T,Tdim = Tdat
    F = S.base_ring()
    assert F == T.base_ring()
    assert len(Sdim) == len(Tdim)
    R = PolynomialRing(F,['%s%d' % (x,i) 
                      for x,a,b in zip('abcdefghijklmnopqrstuvwxyz',Sdim,Tdim) 
                      for i in range(a+b)])
    SplusT = R.zero()
    Sxs = []
    Txs = []
    curd = 0
    for d1,d2 in zip(Sdim,Tdim):
        Sxs.extend(R.gens()[curd:curd+d1])
        Txs.extend(R.gens()[curd+d1:curd+d1+d2])
        curd += d1 + d2
    return S(Sxs) + T(Txs),tuple(d1+d2 for d1,d2 in zip(Sdim,Tdim))


