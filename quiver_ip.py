
def get_flag_inequalities2(m,g,interval_length_bound=oo):
    g = g.transitive_reduction()
    below_cnts = [0 for _ in range(m.ncols())]
    for i,j,_ in g.edges():
        below_cnts[j] += 1
    minelts = set([j for j,cnt in enumerate(below_cnts) if cnt == 0])
    def ideal_push(j):
        minelts.remove(j)
        try:
            for _,k,_ in g.outgoing_edges(j):
                below_cnts[k] -= 1
                if below_cnts[k] == 0:
                    minelts.add(k)
        except:
            embed()
            raise
    def ideal_pop(j):
        for _,k,_ in g.outgoing_edges(j):
            if below_cnts[k] == 0:
                minelts.remove(k)
            below_cnts[k] += 1
        minelts.add(j)
    for rowcutoff in range(max(m.nrows() - interval_length_bound + 1,1)):
        mcur = m[rowcutoff:]
        rowcols = []
        rowcols_sets_seen = set()
        def dfs():
            nonlocal mcur
            if len(rowcols) >= interval_length_bound:
                return
            for j in list(minelts):
                col = mcur.column(j)
                i = next(i for i in range(mcur.nrows()-1,-1,-1) if col[i] != 0)
                rowcols.append((rowcutoff+i,j))
                frozen_rowcols = frozenset(rowcols)
                if frozen_rowcols in rowcols_sets_seen:
                    rowcols.pop()
                    continue
                rowcols_sets_seen.add(frozen_rowcols)
                yield sorted(rowcols)
                row = mcur.row(i)
                col /= col[i]
                ideal_push(j)
                mcur -= col.column() * row.row()
                zeroed = []
                while any(mcur[:,kx := k].is_zero() for k in minelts):
                    ideal_push(kx)
                    zeroed.append(kx)
                for res in dfs():
                    yield res
                while zeroed:
                    ideal_pop(zeroed.pop())
                mcur += col.column() * row.row()
                ideal_pop(j)
                rowcols.pop()
        zeroed = []
        while any(mcur[:,kx := k].is_zero() for k in minelts):
            ideal_push(kx)
            zeroed.append(kx)
        for res in dfs():
            yield res
        while zeroed:
            ideal_pop(zeroed.pop())
    
                
def get_flag_inequalities(ms,interval_length_bound=oo):
    if len(ms) == 0:
        return
    tosize = ms[0].nrows()
    assert all(m.nrows() == tosize for m in ms)
    sol_ss = {}
    ranks = {}
    # for l in range(1,tosize+1):
    # for l in list(range(1,min(interval_length_bound,tosize-1)+1))+[tosize]:
    for l in range(1,min(interval_length_bound,tosize)+1):
        for a in range(tosize-l+1):
            b = a+l
            sol_ss[(a,b)] = SetSystem()
            ss = sol_ss[(a,b)]
            Ms = []
            for M in ms:
                # want largest fb so that the subspace indexed up to fb maps into 
                # the subspace indexed up to b, ie, so that M[b:,:fb].is_zero()
                fb = next(fb for fb in range(M.ncols(),-1,-1) if M[b:,:fb].is_zero())
                Ms.append(M[a:b,:fb])
            r = block_matrix([Ms]).rank()
            ranks[(a,b)] = r
            if r == 0:
                continue
            cols_sets_seen = set()
            cols = []
            def dfs(Ms):
                if len(cols) == r:
                    ss.add_set(cols,None)
                    have = False
                    for l1 in range(1,l):
                        if (a,a+l1) in sol_ss and (a+l1,b) in sol_ss:
                            assert r >= ranks[(a,a+l1)] + ranks[(a+l1,b)]
                            if r == ranks[(a,a+l1)] + ranks[(a+l1,b)]:
                                ss1 = sol_ss[(a,a+l1)]
                                ss2 = sol_ss[(a+l1,b)]
                                if any(True for psol,_ in ss1.iter_sets(Shi=cols) 
                                       for _ in ss2.iter_sets(Shi=[col for col in cols if col not in psol])):
                                    have = True
                                    break
                    if not have:
                        yield deepcopy(cols)
                    return
                for ii,M in enumerate(Ms):
                    if M.is_zero():
                        continue
                    jx = next(j for j,c in enumerate(M.columns()) if not c.is_zero())
                    cols.append((ii,jx))
                    frozen_cols = frozenset(cols)
                    if frozen_cols in cols_sets_seen:
                        cols.pop()
                        continue
                    cols_sets_seen.add(frozen_cols)
                    col = M.column(jx)
                    pivot_row = next(i for i,e in enumerate(col) if e != 0)
                    col /= col[pivot_row]
                    Ms_next = [copy(M) for M in Ms]
                    for mi,M in enumerate(Ms_next):
                        Ms_next[mi] -= col.column() * M.row(pivot_row).row()
                    for full_cols in dfs(Ms_next):
                        yield full_cols
                    cols.pop()
            for full_cols in dfs(Ms):
                yield (list(range(a,b)), full_cols)
                
def ip_to_z3(ip):
    import z3
    s = z3.Solver()
    s.set("sat.pb.solver", "circuit")
    for lo,(xis,alphas),hi in ip.constraints():
        if all(alpha == 1 for alpha in alphas) or all(alpha == -1 for alpha in alphas):
            xs = [z3.Bool('p%d'%xi) for xi in xis]
            if alphas[0] == -1:
                lo,hi = None if hi is None else -hi,None if lo is None else -lo
            if lo is not None:
                s.add( z3.AtLeast(*xs+ [int(lo)]) )
            if hi is not None:
                s.add( z3.AtMost(*xs+ [int(hi)]) )
        else:
            xs = [(z3.Bool('p%d'%xi),int(alpha)) for xi,alpha in zip(xis,alphas)]
            if lo is not None and hi is not None and lo == hi:
                s.add( z3.PbEq(xs, int(hi)) )
            else:
                if lo is not None:
                    s.add( z3.PbGe(xs, int(lo)) )
                if hi is not None:
                    s.add( z3.PbLe(xs, int(hi)) )
    return s

def ip_to_scipy(ip):
    import scipy
    import numpy as np
    I = []
    J = []
    V = []
    lb = []
    ub = []
    for ci,(lo,(xis,alphas),hi) in enumerate(ip.constraints()):
        I.extend([ci]*len(xis))
        J.extend(xis)
        V.extend(alphas)
        lb.append(-np.inf if lo is None else lo)
        ub.append(np.inf if hi is None else hi)
    A = scipy.sparse.coo_matrix((V,(I,J)), (ip.number_of_constraints(), ip.number_of_variables()) )
    constraints = scipy.optimize.LinearConstraint(A,lb=lb,ub=ub)
    lb = np.zeros(ip.number_of_variables())
    ub = np.zeros(ip.number_of_variables())
    for x in ip.default_variable().keys():
        lb[ next(iter(ip[x].dict().keys())) ] = ip.get_min(ip[x])
        ub[ next(iter(ip[x].dict().keys())) ] = ip.get_max(ip[x])
    bounds = scipy.optimize.Bounds(lb=lb,ub=ub)
    integrality = [1]*A.shape[1]
    return scipy.optimize.milp([0]*A.shape[1],
        integrality = integrality, bounds = bounds, constraints = constraints, options = {'disp' : True})
