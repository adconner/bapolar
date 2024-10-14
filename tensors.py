from sage.all import PolynomialRing
from itertools import product

def matrixmult(*idims):
    dims = tuple(a*b for a,b in zip(idims,idims[1:]+idims[:1]))
    R = PolynomialRing(QQ,['%s%d%d' % (c,i,j) for c,a,b in zip('uvwxyzabcdefghijklmnopqrst', idims,idims[1:]+idims[:1]) 
                           for i,j in product(range(a),range(b))])
    xs = []
    c = 0
    for a,b in zip(idims,idims[1:]+idims[:1]):
        xs.append(R.gens()[c:c+a*b])
        c += a*b
    T = sum(prod(u[i*a+j] for u,i,j,a in zip(xs,ixs,ixs[1:]+(ixs[0],),idims[1:]+idims[:1]))
            for ixs in product(*map(range,idims)))
    return (T,dims)

def slpart110(p,n):
    xs = p.parent().gens()
    def mon(m):
        (i,_),(j,_) = list(m.exponents()[0].sparse_iter())
        j -= n*n
        if i % n == j // n:
            return m - sum(xs[(i//n)*n+k]*xs[n*n+k*n+j%n] for k in range(n))/n
        return m
    return sum(mon(m)*p.monomial_coefficient(m) for m in p.monomials())

def slpart111(p,n):
    xs = p.parent().gens()
    def mon(m):
        (i,_),(j,_),(k,_) = list(m.exponents()[0].sparse_iter())
        ixs = (i,j-n*n,k-2*n*n)
        
        nn = sum(
            prod(xs[e*n*n+(ixs[e]//n
                           if jxs[(e+2) % 3] is None else
                           jxs[(e+2)%3])*n+ (ixs[e]%n
                                 if jxs[e] is None
                                 else jxs[e])]
            for e in range(3))
                for jxs in product(*[range(n) if ixs[e] % n == ixs[(e+1)%3] // n else [None] 
                      for e in range(3)]))
        return m - nn/len(nn.monomials())
    return sum(mon(m)*p.monomial_coefficient(m) for m in p.monomials())

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

def det(n):
    m = matrix(n,n,PolynomialRing(QQ,['x%d%d' % (i+1,j+1) for i,j in product(range(n),range(n))]).gens())
    dims = (n**2,)  
    T = m.det()
    vwts = matrix(ZZ,2*(n-1),n**2)
    for i in range(n-1):
        for j in range(n):
            vwts[i, i*n+j] = 1
            vwts[i, (i+1)*n+j] = -1
            vwts[n-1+i, j*n+i] = 1
            vwts[n-1+i, j*n+i+1] = -1
    xs = []
    for i in range(n):
        for j in range(i+1,n):
            x = matrix(QQ,n*n,n*n,sparse=True)
            for k in range(n):
                x[i*n+k, j*n+k] = 1
            xs.append(x)
            x = matrix(QQ,n*n,n*n,sparse=True)
            for k in range(n):
                x[k*n+i, k*n+j] = 1
            xs.append(x)
    return (T,dims,vwts,xs)

def yao(dims):
    F = GF(32003)
    xss = get_defining_variables(F,dims)
    T = sum(F.random_element()*prod(xs[k] for xs,k in zip(xss,ks) ) 
            for ks in IntegerVectors((sum(dims)-3) // 2, len(dims), outer = [e-1 for e in dims]))
    return T,dims

def get_defining_variables(F,dims):
    R = PolynomialRing(F,['%s%d' % (x,i) 
                      for x,b in zip('abcdefghijklmnopqrstuvwxyz',dims)
                      for i in range(b)])
    xss = []
    c = 0
    for d in dims:
        xss.append(R.gens()[c:c+d])
        c += d
    return xss

def hermitian_symmetric_space(a,b):
    X,Y,Z,U = [matrix(a,b,xs) for xs in get_defining_variables(QQ,(a*b,)*4)]
    T = (X*Y.T).trace()* (Z*U.T).trace() - (X*U.T).trace()* (Y*Z.T).trace()
    return T, (a*b,)*4

def octonion():
    triads = [(1,2,4), (1,3,7), (1,5,6), (2,3,5), (2,6,7), (3,4,6), (4,5,7)]
    xs,ys,zs = get_defining_variables(QQ,(7,7,7))
    T = 0
    for ixs in triads:
        ixs = [i-1 for i in ixs]
        for a,b,c in zip(ixs,ixs[1:]+ixs[:1],ixs[2:]+ixs[:2]):
            T -= xs[a]*ys[b]*zs[c]
        ixs = [ixs[0],ixs[2],ixs[1]]
        for a,b,c in zip(ixs,ixs[1:]+ixs[:1],ixs[2:]+ixs[:2]):
            T += xs[a]*ys[b]*zs[c]
    return T,(7,7,7)
    
    
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
