    
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
    from functools import cache
    sets = {}
    ei = -1
    @cache
    def key(s):
        nonlocal ei
        ei += 1
        return ei
    def add_set(S,label):
        S = sorted(S, key=key)
        cur_sets = sets
        for s in S:
            cur_sets = cur_sets.setdefault(s,{})
        cur_sets[None] = label
    from collections import Counter
    from copy import copy
    from heapq import heapify,heappush,heappop
    def iter_sets(Slo,Shi):
        Slo = [(key(s),s) for s in Slo]
        heapify(Slo)
        Shi = Counter(Shi)
        cur = []
        def dfs(cur_sets):
            if None in cur_sets and len(Slo) == 0:
                yield (copy(cur),cur_sets[None])
            for e,sets1 in cur_sets.items():
                if e is None:
                    continue
                if Shi.get(e,0) > 0 and (
                        len(Slo) == 0 or key(e) <= Slo[0][0]):
                    cur.append(e)
                    Shi[e] -= 1
                    if len(Slo) > 0 and key(e) == Slo[0][0]:
                        topush = heappop(Slo)
                    else:
                        topush = None
                    for r in dfs(sets1):
                        yield r
                    if topush is not None:
                        heappush(Slo,topush)
                    Shi[e] += 1
                    cur.pop()
        return dfs(sets)
    return add_set,iter_sets

def lp_integer_points(lp,xs=None,fullsol=True,prunef=lambda psol: True):
    from copy import deepcopy
    from sage.numerical.mip import MIPSolverException
    lp = deepcopy(lp)
    if xs is None:
        xs = list(lp.default_variable().keys())
    sol = {x : lp.get_min(lp[x]) for x in xs 
           if lp.get_min(lp[x]) == lp.get_max(lp[x])}
    st = []
    def dfs(check_solvable=True):
        try:
            if check_solvable:
                lp.solve()
        except MIPSolverException:
            return
        csol = lp.get_values(lp.default_variable())
        remxs = [x for x in xs if lp.get_min(lp[x]) < lp.get_max(lp[x])]
        if len(remxs) == 0:
            if fullsol:
                yield csol
            else:
                yield copy(sol)
            return
        x = max(remxs, key=lambda x: abs(csol[x]-round(csol[x])))
        v = floor(csol[x]+1e-10)
        omin = lp.get_min(lp[x])
        omax = lp.get_max(lp[x])
        if v == omax:
            v -= 1
        assert v >= lp.get_min(lp[x])
        assert v+1 <= lp.get_max(lp[x])
        print('%s%s %d %.2f %d' % (''.join(st),str(x),lp.get_min(lp[x]),
                                 csol[x],lp.get_max(lp[x])))
        hi_first = csol[x] - v >= 0.5
        if hi_first:
            lp.set_min(lp[x],v+1)
            st.append(' ')
            if v+1 == lp.get_max(lp[x]):
                sol[x] = v+1
                if prunef(sol):
                    for res in dfs(abs(csol[x]-(v+1)) > 1e-10):
                        yield res
                del sol[x]
            else:
                for res in dfs(abs(csol[x]-(v+1)) > 1e-10):
                    yield res
            lp.set_min(lp[x],omin)
            st[-1] = '.'
        else:
            st.append(' ')
        lp.set_max(lp[x],v)
        if v == lp.get_min(lp[x]):
            sol[x] = v
            if prunef(sol):
                for res in dfs(abs(csol[x]-v) > 1e-10):
                    yield res
            del sol[x]
        else:
            for res in dfs(abs(csol[x]-v) > 1e-10):
                yield res
        lp.set_max(lp[x],omax)
        if not hi_first:
            st[-1] = '.'
            lp.set_min(lp[x],v+1)
            if v+1 == lp.get_max(lp[x]):
                sol[x] = v+1
                if prunef(sol):
                    for res in dfs(abs(csol[x]-(v+1)) > 1e-10):
                        yield res
                del sol[x]
            else:
                for res in dfs(abs(csol[x]-(v+1)) > 1e-10):
                    yield res
            lp.set_min(lp[x],omin)
        st.pop()
    return dfs()
        
