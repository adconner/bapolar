
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
                    Ms_next = deepcopy(Ms)
                    for mi,M in enumerate(Ms_next):
                        Ms_next[mi] -= col.column() * M.row(pivot_row).row()
                    for full_cols in dfs(Ms_next):
                        yield full_cols
                    cols.pop()
            for full_cols in dfs(Ms):
                yield (list(range(a,b)), full_cols)
                
