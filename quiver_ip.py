from functools import cache

def echelonize_graph(m,g):
    m = copy(m)
    below = {}
    pivots = {}
    for j in g.topological_sort():
        if g.in_degree(j) > 0:
            M = block_matrix([[below[j2]] for j2,_,_ in g.incoming_edges(j)] )
            M.echelonize()
            M = M[:M.rank()]
            pivots[j] = M.pivots()
        else:
            M = matrix(m.base_ring(),0,m.nrows())
            pivots[j] = ()
        col = m.column(j)
        col -= vector([col[i] for i in pivots[j]])*M
        m[:,j] = col.column()
        M = M.T.augment(col.column()).T
        M.echelonize()
        below[j] = M
    pivots = [pivots[j] for j in range(m.ncols())]
    return m,pivots

def get_flag_inequalities2(m,g,relation_size_bound=oo):
    # m,pivots = echelonize_graph(m[::-1],g)
    # m = m[::-1]
    # pivots = [[m.nrows()-1-j for j in jxs] for jxs in pivots]
    
    m2 = block_matrix([[m,identity_matrix(m.base_ring(),m.nrows())]])
    g2 = copy(g)
    g2.add_vertex(m.ncols())
    for j in range(m.ncols(),m2.ncols()-1):
        g2.add_edge(j,j+1)

    lastnzs = [None for _ in range(m2.ncols())]
    for j in g2.topological_sort():
        col = m2.column(j)
        lastnz = 0 if col.is_zero() else next(len(col)-i for i,e in enumerate(col[::-1]) if e != 0)
        lastnzs[j] = max(lastnz,max((lastnzs[j2] for j2,_,_ in g2.incoming_edges(j)),default=0))
        
    # lp = MixedIntegerLinearProgram(solver="SCIP")
    # lp.set_binary(lp.default_variable())
    # vs = {}
    # for j,r in product(range(m2.ncols()),range(m2.nrows())):
    #     if r <= lastnzs[j]-1:
    #         vs[(j,r)] = lp[(j,r)]
    #     else:
    #         vs[(j,r)] = lp.linear_functions_parent().zero()
    # for j in range(m2.ncols()):
    #     lp.add_constraint(lp.sum(vs[(j,r)] for r in range(lastnzs[j])) == lp[j])
    # for (_,j1),(_,j2) in combinations(sorted([(nz,i) for i,nz in enumerate(lastnzs)]),2):
    #     for r2 in range(lastnzs[j2]):
    #         for r1 in range(r2+1,lastnzs[j1]):
    #             lp.add_constraint( vs[(j1,r1)] + vs[(j2,r2)] <= 1 )
    # for r in range(m2.nrows()):
    #     lp.add_constraint(lp.sum(vs[(j,r)] for j in range(m2.ncols())) <= 1)
    # for r1,r2 in combinations(range(m2.nrows()),2):
    #     lp.add_constraint(lp.sum(vs[(j,r1)] for j in range(m2.ncols())) >= 
    #                       lp.sum(vs[(j,r2)] for j in range(m2.ncols())))
    # lp.add_constraint(lp.sum(lp[j] for j in range(m2.ncols())) <= relation_size_bound)
    # lp.add_constraint(lp.sum(vs[(j,lastnzs[j]-1)] for j in range(m.ncols()) if lastnzs[j] > 0) == 1)
    # import pyscipopt
    # lp.get_backend()._get_model().setEmphasis(pyscipopt.SCIP_PARAMEMPHASIS.FEASIBILITY)
    # from sage.numerical.mip import MIPSolverException
    # def prunef(jxs,jxsout):
    #     lp.get_backend()._get_model().freeTransform()
    #     lp.set_min(lp.default_variable(),0)
    #     lp.set_max(lp.default_variable(),1)
    #     for j in jxs:
    #         lp.set_min(lp[j],1)
    #     for j in jxsout:
    #         lp.set_max(lp[j],0)
    #     # lp.get_backend()._get_model().hideOutput(False)
    #     try:
    #         lp.solve()
    #         return True
    #     except MIPSolverException:
    #         return False
        
    # from pyscipopt import quicksum,Model
    # M = Model()
    # M.enableReoptimization()
    # vs = {}
    # for j,r in product(range(m2.ncols()),range(m2.nrows())):
    #     if r <= lastnzs[j]-1:
    #         vs[(j,r)] = M.addVar(str((j,r)),'B')
    #     else:
    #         vs[(j,r)] = 0
    # lpcols = []
    # for j in range(m2.ncols()):
    #     lpcols.append(quicksum(vs[(j,r)] for r in range(lastnzs[j])))
    #     M.addCons(lpcols[-1] <= 1)
    # for (_,j1),(_,j2) in combinations(sorted([(nz,i) for i,nz in enumerate(lastnzs)]),2):
    #     for r2 in range(lastnzs[j2]):
    #         for r1 in range(r2+1,lastnzs[j1]):
    #             M.addCons( vs[(j1,r1)] + vs[(j2,r2)] <= 1 )
    # for r in range(m2.nrows()):
    #     M.addCons(quicksum(vs[(j,r)] for j in range(m2.ncols())) <= 1)
    # for r1,r2 in combinations(range(m2.nrows()),2):
    #     M.addCons(quicksum(vs[(j,r1)] for j in range(m2.ncols())) >= 
    #               quicksum(vs[(j,r2)] for j in range(m2.ncols())))
    # M.addCons(quicksum(lpcols[j] for j in range(m2.ncols())) <= relation_size_bound)
    # M.addCons(quicksum(vs[(j,lastnzs[j]-1)] for j in range(m.ncols()) if lastnzs[j] > 0) == 1)
    # def prunef(jxs,jxsout):
    #     M.freeReoptSolve()
    #     M.chgReoptObjective(quicksum(lpcols[j] for j in jxs) - quicksum(lpcols[j] for j in jxsout),"maximize")
    #     M.hideOutput(True)
    #     M.optimize()
    #     return M.getStatus() == 'optimal' and int(M.getObjVal()) == len(jxs)
    
    def prunef(jxs,jxsout):
        nzs = [(lastnzs[j],j) for j in jxs]
        nzs.sort()
        complete = False
        for i,(nz,j) in enumerate(nzs):
            assert nz >= i+1
            if complete and j < m.ncols():
                return False
            if nz == i+1 and j < m.ncols():
                complete = True
        return True
        
    for jxs in get_fillings(m2,g2,[1]*m.ncols()+[0]*m.nrows(),relation_size_bound,prunef,lambda j: lastnzs[j]):
        rows = set(range(m.nrows()))
        cols = []
        for j in jxs:
            if j < m.ncols():
                cols.append(j)
            else:
                rows.remove(j - m.ncols())
        yield (list(rows),cols)

def get_fillings(m,g,costs=None,cost_bound=0,prunef=lambda cols,skip: True,sortkey = None):
    assert(len(g) == m.ncols())
    if costs is None:
        costs = [0]*m.ncols()
    g = g.transitive_reduction()
    below_cnts = [0 for _ in range(m.ncols())]
    for i,j,_ in g.edges(sort=False):
        below_cnts[j] += 1
    minelts = set([j for j,cnt in enumerate(below_cnts) if cnt == 0])
    def ideal_push(j):
        minelts.remove(j)
        for _,k,_ in g.outgoing_edges(j):
            below_cnts[k] -= 1
            if below_cnts[k] == 0:
                minelts.add(k)
    def ideal_pop(j):
        for _,k,_ in g.outgoing_edges(j):
            if below_cnts[k] == 0:
                minelts.remove(k)
            below_cnts[k] += 1
        minelts.add(j)
    while any(m[:,kx := k].is_zero() for k in minelts):
        ideal_push(kx)
    cols = []
    skip = set()
    def dfs(c):
        nonlocal m
        if len(minelts) == 0:
            yield list(cols)
        skiphere = []
        for j in sorted(minelts-skip,key=sortkey):
            try:
                cols.append(j)
                if c + costs[j] > cost_bound or not prunef(cols,skip):
                    cols.pop()
                    continue
                col = m.column(j)
                i = next(i for i,e in enumerate(col) if e != 0)
                row = m.row(i)
                col /= col[i]
                ideal_push(j)
                m -= col.column() * row.row()
                zeroed = []
                while any(m[:,kx := k].is_zero() for k in minelts):
                    ideal_push(kx)
                    zeroed.append(kx)
                for res in dfs(c+costs[j]):
                    yield res
                while zeroed:
                    ideal_pop(zeroed.pop())
                m += col.column() * row.row()
                ideal_pop(j)
                cols.pop()
            finally:
                skip.add(j)
                skiphere.append(j)
        for j in skiphere:
            skip.remove(j)
    return dfs(0)
    
                
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
