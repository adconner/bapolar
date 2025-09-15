# For koszul flattening compuations, three place tensors in A ot B ot C are
# represented by their contractions with a basis of A^*, ie lists of matrices
# representing elements of B ot C. To distinguish three place tensors
# represented this way with other tensors considered by the border apolarity
# code, the letter L is used for such tensors (meant to suggest viewing T as
# the linear space L = T(A^*))

# To use koszul flattenings to obtain lower bounds for a tensor T considered in
# this code, possibly of more than three places, it must be flattenend and
# converted to linear subspace of matrices form. The 
def flatten_to_matrix_subspace(T,row_factor_indices,col_factor_indices):
    p,dims = T
    row_factor_indices = sorted(row_factor_indices)
    col_factor_indices = sorted(col_factor_indices)
    assert all(i < len(dims) for i in row_factor_indices+col_factor_indices)
    assert len(set(row_factor_indices).intersection(set(col_factor_indices))) == 0
    Afactor_indices = sorted(set(range(len(dims))) - set(row_factor_indices) - set(col_factor_indices))
    deg = tensor_deg(T)
    def multi_subset_rank(s):
        return sum(binomial(e+i,i+1) for i,e in enumerate(sorted(s)))
    L = [{} for _ in range(prod(binomial(dims[i]+deg[i]-1,deg[i]) for i in Afactor_indices))]
    for ks,e in Tdict(T):
        Ai = 0
        for i in Afactor_indices:
            Ai *= binomial(dims[i]+deg[i]-1,deg[i])
            Ai += multi_subset_rank(ks[i])
        rowi = 0
        for i in row_factor_indices:
            rowi *= binomial(dims[i]+deg[i]-1,deg[i])
            rowi += multi_subset_rank(ks[i])
        coli = 0
        for i in col_factor_indices:
            coli *= binomial(dims[i]+deg[i]-1,deg[i])
            coli += multi_subset_rank(ks[i])
        L[Ai][(rowi,coli)] = e
    return [ matrix(p.base_ring(),
                    prod(binomial(dims[i]+deg[i]-1,deg[i]) for i in row_factor_indices),
                    prod(binomial(dims[i]+deg[i]-1,deg[i]) for i in col_factor_indices),m) 
            for m in L ]

def TAp(L,p):
    a=len(L)
    b,c = L[0].dimensions()
    F = L[0].base_ring()
    That = {}
    from itertools import combinations
    def subset_rank(s):
        return sum(binomial(e,i+1) for i,e in enumerate(sorted(s)))
    for ii,m in enumerate(L):
        for S in combinations(chain(range(ii),range(ii+1,a)),p):
            P = S+(ii,)
            sg = (-1)**len([_ for i in S if i < ii])
            br = subset_rank(P)
            bc = subset_rank(S)
            for (i,j),e in m.dict().items():
                That[(br*c + j,bc*b + i)] = sg * e
    return matrix(F,binomial(a,p+1)*c, binomial(a,p)*b,That)

# M is newa x a matrix
def tensor_restrict(L,M):
    assert M.ncols() == len(L)
    Tr = [{} for i in range(M.nrows())]
    for (i,j),e in M.dict().items():
        for (k,l),f in L[j].dict().items():
            Tr[i][(k,l)] = Tr[i].get((k,l),0) + e*f
    Tr = [ matrix(L[0].base_ring(),
        L[0].nrows(),L[0].ncols(),m) for m in Tr]
    return Tr

def generic_restrict(L,newa):
    F = L[0].base_ring()
    if F.is_finite():
        M = random_matrix(F,newa,len(L))
    else:
        M = random_matrix(ZZ,newa,len(L),x=-1000,y=1000)
    return tensor_restrict(L,M)

def koszul_lower_bound(T,p,known_bound = 0):
    Tp,dims = T
    # change to modular ring to improve rank computation speed
    T = (Tp.change_ring(GF(32003)),dims) 
    deg = tensor_deg(T)
    factordims = [binomial(dmax+dcnt-1, dcnt) for dmax,dcnt in zip(dims, deg)]
    def bound_upper_bound(Bixs,Cixs,newa):
        TApnrows = prod(factordims[i] for i in Bixs)*binomial(newa,p)
        TApncols = prod(factordims[i] for i in Cixs)*binomial(newa,p+1)
        return min(TApnrows, TApncols) / binomial(newa-1,p)
    ordered_set_partitions = [ [[j for j,e in enumerate(facmapping) if e==i ] for i in range(3)]
           for facmapping in product(*[range(3) if p > 0 else range(1,3)]*len(dims)) ]
    ordered_set_partitions = [Ss for Ss in ordered_set_partitions if all(len(s) > 0 for s in Ss)] if p > 0 else [Ss for Ss in ordered_set_partitions if all(len(s) > 0 for s in Ss[1:])]
    tests = [(Bixs,Cixs,newa,a) for Aixs,Bixs,Cixs in ordered_set_partitions
             for a in [prod(factordims[i] for i in Aixs)]
             for newa in range(p+1 if p > 0 else a, a+1) ]
    tests.sort(key = lambda rec: bound_upper_bound(rec[0],rec[1],rec[2]), reverse = True)
    for Bixs,Cixs,newa,a in tests:
        if bound_upper_bound(Bixs,Cixs,newa) <= known_bound:
            break
        L = flatten_to_matrix_subspace(T, Bixs, Cixs)
        print(f'p={p} ({len(L)}, {L[0].nrows()}, {L[0].ncols()}) ',end='')
        if newa < a:
            print(f'-> ({newa}, {L[0].nrows()}, {L[0].ncols()}) ', end='')
            L = generic_restrict(L, newa)
        sys.stdout.flush()
        M = TAp(L,p)
        print(f'TAp ({M.nrows()},{M.ncols()}) ',end='')
        sys.stdout.flush()
        r = M.rank()
        bound = ceil(r / binomial(newa-1,p))
        print(f'rank {r} bound {bound}')
        known_bound = max(bound,known_bound)
    return known_bound

def summarize(T,pmin=0,pmax=10,known_bound=0):
    for p in range(pmin,pmax+1):
        known_bound = koszul_lower_bound(T,p,known_bound)
    return known_bound
