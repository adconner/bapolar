# from sage.all import *

def diagonal_torus(T):
    eqs = matrix(ZZ,[e for e in T.exponents()])
    tights = eqs.right_kernel_matrix(basis='LLL')
    return tights.change_ring(T.base_ring())
    # # equivalently,
    # return liealgebraeqs(T, (1,)*T.parent().ngens()).right_kernel_matrix()

def liealgebra(T,dims = None):
    if dims == None:
        dims = (T.parent().ngens(),)
    eqs = liealgebraeqs(T,dims)
    K = eqs.right_kernel_matrix()
    L = []
    for r in K:
        cur = []
        lie_start = 0
        for d in dims:
            cur.append(matrix(T.base_ring(),d,d,r[lie_start:lie_start+d*d]))
            lie_start += d*d
        L.append(tuple(cur) if len(dims) > 1 else cur[0])
    return L

# gives equations for the lie algebra for T acting on the left. dims =
# (d_1,...,d_i) lets one optionally restrict to the d_i x d_i diagonal blocks
# of the lie algebra, which is useful, for instance, if T is multihomogeneous.
def liealgebraeqs(T,dims = None):
    R = T.parent()
    if dims == None:
        dims = (R.ngens(),)
    ms = []
    msix = {}
    eqs = {}
    lie_start = 0
    xs_start = 0
    for dim in dims:
        xs = R.gens()[xs_start:xs_start+dim]
        for (i,x),(j,y) in product(enumerate(xs),enumerate(xs)):
            Tcur = x*T.derivative(y)
            for m,c in Tcur.dict().items():
                if m not in msix:
                    msix[m] = len(ms)
                    ms.append(m)
                eqs[(msix[m] ,lie_start + i*dim + j)] = c
        xs_start += dim
        lie_start += dim*dim
    eqs = matrix(R.base_ring(),len(ms), lie_start, eqs)
    return eqs


# finds a large solvable lie algebra of T which can be used in bapolar. The
# found subalgebra respects given coordinates of T and uses the diagonal torus
# in this basis (so better may choices may exist after a change of basis).
def solvable_lie_algebra(T,stop_better_than=None,uppertri=False,search_bound=100):
    T,dims = T
    F = T.base_ring()
    vwts = diagonal_torus(T)
    eqs = liealgebraeqs(T,dims)
    
    start = 0
    lie_start = 0
    lie_weights = matrix(F,vwts.nrows(),sum(d*d for d in dims))
    for d in dims:
        for i,j in product(range(d),range(d)):
            if uppertri and i >= j:
                continue
            lie_weights[:,lie_start+i*d+j] = vwts[:, start+i] - vwts[:,start+j]
        start += d
        lie_start += d*d
        
    wts = {}
    for i,r in enumerate(lie_weights.columns()):
        wts.setdefault(r,[]).append(i)
    Mcols = []
    for wt,jxs in wts.items():
        if wt != 0:
            for _ in range(len(jxs) - eqs[:,jxs].rank()):
                Mcols.append(wt)
    if len(Mcols) == 0:
        return T,dims,vwts,[]
    M = matrix(Mcols).T
    sol = vector(best_hyperplane(M,stop_better_than,search_bound)[0])
    raising_wts = sorted(set([v for v in M.columns() if v.dot_product(sol) > 0]))
    xs = []
    for wt in raising_wts:
        jxs = wts[wt]
        K = eqs[:,jxs].right_kernel_matrix()
        for r in K:
            x = matrix(F,sum(dims),sum(dims))
            for i,e in r.dict().items():
                i = jxs[i]
                j = 0
                while i >= dims[j]*dims[j]:
                    i -= dims[j]*dims[j]
                    j += 1
                x[sum(dims[:j])+i//dims[j],sum(dims[:j])+i%dims[j]] = e
            xs.append(x)
    return T,dims,vwts,xs

