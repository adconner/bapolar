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
    ix_sets = [ Subsets([i for i in range(dmax) for _ in range(dcnt)],dcnt,submultiset=True) 
               for dcnt,dmax in zip(deg,dims) ]
    L = [{} for _ in range(prod(ix_sets[i].cardinality() for i in Afactor_indices))]
    for deg,e in Tdict(T):
        Ai = 0
        for i in Afactor_indices:
            Ai *= ix_sets[i].cardinality()
            Ai += ix_sets[i].rank(deg[i])
        rowi = 0
        for i in row_factor_indices:
            rowi *= ix_sets[i].cardinality()
            rowi += ix_sets[i].rank(deg[i])
        coli = 0
        for i in col_factor_indices:
            coli *= ix_sets[i].cardinality()
            coli += ix_sets[i].rank(deg[i])
        L[Ai][(rowi,coli)] = e
    return [ matrix(p.base_ring(),
                    prod(ix_sets[i].cardinality() for i in row_factor_indices),
                    prod(ix_sets[i].cardinality() for i in col_factor_indices),m) 
            for m in L ]

def TAp(L,p):
    a=len(L)
    b,c = L[0].dimensions()
    F = L[0].base_ring()
    That = {}
    try:
        left_ixs = Subsets(range(a),p)
    except:
        embed()
        raise
        
    right_ixs = Subsets(range(a),p+1)
    for ii,m in enumerate(L):
        for S in Subsets(chain(range(ii),range(ii+1,a)),p):
            P = S+Set([ii])
            sg = (-1)^len([_ for i in S if i < ii])
            br = right_ixs.rank(P)
            bc = left_ixs.rank(S)
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

def summary1(L,p,restrict=True):
    if p > len(L)-1:
        raise ValueError("koszul flattening: dim A=%d too small for p=%d" % (len(L),p))

    if not restrict or 2*p+1 >= len(L):
        if restrict:
            print ('warning: TAp: a = %d <= %d = 2*p+1' % (len(L),2*p+1))
        r = TAp(L,p).rank()
        d = binomial(len(L)-1,p)
        return [p,r,d,ceil(r/d)]

    m = TAp(generic_restrict(L,2*p+1),p)
    r = m.rank()
    d = binomial(2*p,p)
    return [p,r,d,ceil(r/d)]

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
    tests = [(Bixs,Cixs,newa,a) for Aixs,Bixs,Cixs in OrderedSetPartitions(range(len(dims)),3) 
             for a in [prod(factordims[i] for i in Aixs)]
             for newa in range(p+1, a+1) ]
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

def summarize(T,pmin=1,pmax=3,known_bound=0):
    for p in range(pmin,pmax+1):
        known_bound = koszul_lower_bound(T,p,known_bound)
    print(f'best bound found {known_bound}')
