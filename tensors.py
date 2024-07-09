from sage.all import PolynomialRing
from itertools import product

def matrixmult(*idims):
    dims = tuple(a*b for a,b in zip(idims,idims[1:]+(idims[0],)))
    R = PolynomialRing(QQ,['%s%d%d' % (c,i,j) for c,a,b in zip('uvwxyzabcdefghijklmnopqrst', idims,idims[1:]+(idims[0],)) 
                           for i,j in product(range(a),range(b))])
    xs = []
    c = 0
    for a,b in zip(idims,idims[1:]+(idims[0],)):
        xs.append(R.gens()[c:c+a*b])
        c += a*b
    T = sum(prod(u[i*a+j] for u,i,j,a in zip(xs,ixs,ixs[1:]+(ixs[0],),idims[1:]+(idims[0],)))
            for ixs in product(*map(range,idims)))
    return (T,dims)

def matrixmult_ti(m,n,l):
    T,dims = matrixmult(m,n,l)
    vwts = matrix(ZZ,n+l+m-1,T.parent().ngens())
    for i in range(n):
        for j in range(m):
            vwts[i,j*n+i] = 1
        for k in range(l):
            vwts[i,m*n+i*l+k] = -1
    for i in range(l):
        for j in range(n):
            vwts[n+i,m*n+j*l+i] = 1
        for k in range(m):
            vwts[n+i,m*n+n*l+i*m+k] = -1
    for i in range(m-1):
        for j in range(l):
            vwts[n+l+i,m*n+n*l+j*m+i] = 1
        for k in range(n):
            vwts[n+l+i,i*n+k] = -1
    xs = []
    for i in range(m):
        for j in range(i+1,m):
            x = matrix(QQ,m*n+n*l+l*m,m*n+n*l+l*m,sparse=True)
            for k in range(l):
                x[m*n+n*l+k*m+i, m*n+n*l+k*m+j] = 1
            for k in range(n):
                x[j*n+k, i*n+k] = -1
            xs.append(x)
    for i in range(n):
        for j in range(i+1,n):
            x = matrix(QQ,m*n+n*l+l*m,m*n+n*l+l*m,sparse=True)
            for k in range(m):
                x[k*n+i, k*n+j] = 1
            for k in range(l):
                x[m*n+j*l+k, m*n+i*l+k] = -1
            xs.append(x)
    for i in range(l):
        for j in range(i+1,l):
            x = matrix(QQ,m*n+n*l+l*m,m*n+n*l+l*m,sparse=True)
            for k in range(n):
                x[m*n+k*l+i, m*n+k*l+j] = 1
            for k in range(m):
                x[m*n+n*l+j*m+k, m*n+n*l+i*m+k] = -1
            xs.append(x)
    return (T,dims,vwts,xs)

def tensor_W():
    dims = (2,2,2)
    R = PolynomialRing(QQ,['%s%d'%(x,i) for x,d in zip('abc',dims) for i in range(d)])
    T = R.gen(0)*R.gen(3)*R.gen(5) + R.gen(1)*R.gen(2)*R.gen(5) + R.gen(1)*R.gen(3)*R.gen(4)
    return T,dims

def skewcw(q=2):
    assert q % 2 == 0
    dims=(q+1,q+1,q+1)
    R = PolynomialRing(QQ,['%s%d'%(x,i) for x,d in zip('abc',dims) for i in range(d)])
    T = sum([sigma.sign()*prod(R.gen(i*(q+1)+
        (0 if sigma(i+1) == 1 else 1+(q//2)*(sigma(i+1)-2)+rho )) for i in range(3))
             for sigma in SymmetricGroup(3) for rho in range(q//2)])
    return T,dims

def cw(q=2):
    assert q % 2 == 0
    dims=(q+1,q+1,q+1)
    R = PolynomialRing(QQ,['%s%d'%(x,i) for x,d in zip('abc',dims) for i in range(d)])
    T = sum([prod(R.gen(i*(q+1)+
        (0 if sigma(i+1) == 1 else 1+(q//2)*(sigma(i+1)-2)+rho )) for i in range(3))
             for sigma in SymmetricGroup(3) for rho in range(q//2)])
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

def tensor_random_permute_coordinates(T):
    dims = T[1]
    sigmas = [SymmetricGroup(d).random_element() for d in dims]
    R = PolynomialRing(QQ,['%s%d'%(x,i) for x,d in zip('abcdefghijklmnopqrstuvwxyz',dims) for i in range(d)])
    subs = [R.gen(sum(dims[:i]) + sigmas[i](j+1)-1) for i in range(len(dims)) for j in range(dims[i])]
    return (T[0](subs), T[1])
    

unextendible_supports_333_dat = [
    [(1,1,3), (1,2,2), (2,1,2), (3,3,1)],
    [(1,1,3), (1,3,2), (2,3,1), (3,2,2)],
    [(1,1,3), (1,2,2), (1,3,1), (2,1,2), (3,2,1)],
    [(1,1,3), (1,2,2), (2,1,2), (2,3,1), (3,2,1)],
    [(1,1,3), (1,2,2), (2,3,1), (3,1,2), (3,2,1)],
    [(1,1,3), (1,3,2), (2,2,2), (3,1,2), (3,3,1)],
    [(1,1,3), (1,2,2), (1,3,1), (2,1,2), (2,2,1), (3,1,1)],
    [(1,1,3), (1,3,2), (2,2,2), (2,3,1), (3,1,2), (3,2,1)],
    [(1,2,3), (1,3,2), (2,1,3), (2,2,2), (2,3,1), (3,1,2), (3,2,1)]]
def unextendible_supports_333():
    dims = (3,3,3)
    R = PolynomialRing(QQ,'x',9)
    Ts = []
    for supp in unextendible_supports_333_dat:
        supp = [[i-1 for i in ix] for ix in supp]
        Ts.append((sum(R.gen(i)*R.gen(3+j)*R.gen(6+k) for i,j,k in supp),dims))
    T,_ = Ts.pop()
    Ts.append((T - 2*R.gen(2)*R.gen(3+1)*R.gen(6+0),dims))
    Ts.append((T,dims))
    return Ts

def fromT(T):
    dims = (len(T),) + T[0].dimensions()
    R = PolynomialRing(T[0].base_ring(),['%s%d' % (c,i+1) 
                 for c,d in zip('abc',dims) for i in range(d)])
    return (R.sum([e*R.gen(i)*R.gen(dims[0]+j)*R.gen(dims[0]+dims[1]+k) for i,m in enumerate(T) for (j,k),e in m.dict().items()]),dims)
