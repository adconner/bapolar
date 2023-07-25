
def lproduct(*its):
    its = [iter(it) for it in its]
    pools = [[] for e in its]
    for k in count():
        added = False
        for i,(pool,it) in enumerate(zip(pools,its)):
            try:
                pool.append(next(it))
                added = True
                for res in product(*pools[:i]+[pool[-1:]]+pools[i+1:]):
                    yield res
            except StopIteration:
                if len(pool) == 0:
                    return
        if not added:
            return
        
def get_set_system():
    sets = {}
    set_universe = []
    set_universe_ix = {}
    def add_set(S,label):
        for s in S:
            if s not in set_universe_ix:
                set_universe_ix[s] = len(set_universe)
                set_universe.append(s)
        S = sorted(S, key = lambda s: set_universe_ix[s])
        cur_sets = sets
        for s in S:
            cur_sets = cur_sets.setdefault(s,{})
        cur_sets[None] = label
    from collections import Counter
    from copy import copy
    def iter_subsets(S):
        S = Counter(S)
        cur = []
        def dfs(cur_sets):
            if None in cur_sets:
                yield (copy(cur),cur_sets[None])
            for e,sets1 in cur_sets.items():
                if e is None:
                    continue
                if S.get(e,0) > 0:
                    cur.append(e)
                    S[e] -= 1
                    for r in dfs(sets1):
                        yield r
                    S[e] += 1
                    cur.pop()
        return dfs(sets)
    return add_set,iter_subsets

def lp_integer_points(lp,xs=None):
    from copy import deepcopy
    from sage.numerical.mip import MIPSolverException
    lp = deepcopy(lp)
    if xs is None:
        xs = list(lp.default_variable().keys())
    def dfs():
        try:
            lp.solve()
        except MIPSolverException:
            return
        remxs = [x for x in xs if lp.get_min(lp[x]) < lp.get_max(lp[x])]
        sol = lp.get_values(lp.default_variable())
        if len(remxs) == 0:
            yield sol
            return
        x = max(remxs, key=lambda x: abs(sol[x]-round(sol[x])))
        print('%s%s %d %.2f %d' % (' '*(len(xs)-len(remxs)),str(x),lp.get_min(lp[x]),
                                 sol[x],lp.get_max(lp[x])))
        v = floor(sol[x]+1e-10)
        omax = lp.get_max(lp[x])
        if v == omax:
            v -= 1
        hi_first = sol[x] - v >= 0.5
        if hi_first:
            omin = lp.get_min(lp[x])
            lp.set_min(lp[x],v+1)
            for sol in dfs():
                yield sol
            lp.set_min(lp[x],omin)
        lp.set_max(lp[x],v)
        for sol in dfs():
            yield sol
        lp.set_max(lp[x],omax)
        if not hi_first:
            omin = lp.get_min(lp[x])
            lp.set_min(lp[x],v+1)
            for sol in dfs():
                yield sol
            lp.set_min(lp[x],omin)
    return dfs()
        
