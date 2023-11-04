from copy import deepcopy

# reimplementation of product from itertools which reorders its output to
# achieve maximum laziness 
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
        
# creates a data structure representing families of submultisets of some
# universe set with two functions add_set and iter_sets. add_set adds a
# multiset to the collection and iter_sets returns an iterator over all
# multisets in the collection containing Slo and contained in Shi. Elements of
# sets must be hashable
class SetSystem:
    def __init__(self):
        self.sets = {}
        self.ei = -1
        from functools import cache
        self.key = cache(self.key)
    def key(self,s):
        self.ei += 1
        return self.ei
    def add_set(self,S,label=None):
        S = sorted(S, key=self.key)
        cur_sets = self.sets
        for s in S:
            cur_sets = cur_sets.setdefault(s,{})
        cur_sets[None] = label
    def remove_set(self,S):
        S = sorted(S, key=self.key)
        def rs(cur_sets,i):
            if i == len(S):
                del cur_sets[None]
            else:
                rs(cur_sets[S[i]], i+1)
                if len(cur_sets[S[i]]) == 0:
                    del cur_sets[S[i]]
        rs(self.sets,0)
    def iter_sets(self,Slo=[],Shi=None,Slolex=[]):
        from collections import Counter
        from copy import copy
        from heapq import heapify,heappush,heappop
        Slo = [self.key(s) for s in Slo]
        heapify(Slo)
        Slolex = [s for s in Slolex]
        heapify(Slolex)
        if Shi is not None:
            Shi = Counter(Shi)
        cur = []
        def dfs(cur_sets,Slolex):
            if None in cur_sets and len(Slo) == 0 and len(Slolex) == 0:
                yield (copy(cur),cur_sets[None])
            for e,sets1 in cur_sets.items():
                if e is None:
                    continue
                if (Shi is None or Shi.get(e,0) > 0) and \
                        (len(Slo) == 0 or self.key(e) <= Slo[0]) and \
                        (len(Slolex) == 0 or Slolex[0] <= e):
                    cur.append(e)
                    if Shi is not None:
                        Shi[e] -= 1
                    topushSlo = heappop(Slo) if len(Slo) > 0 and self.key(e) == Slo[0] else None
                    if len(Slolex) > 0 and e == Slolex[0]:
                        topushSlolex = heappop(Slolex)
                        for r in dfs(sets1,Slolex):
                            yield r
                        heappush(Slolex,topushSlolex)
                    else:
                        for r in dfs(sets1,[]):
                            yield r
                    if topushSlo is not None:
                        heappush(Slo,topushSlo)
                    if Shi is not None:
                        Shi[e] += 1
                    cur.pop()
        return dfs(self.sets,Slolex)
    def iter_unions(self,Shi,replacement=False):
        Shi = Counter(Shi)
        prev = []
        for s,k in self.iter_sets(Shi=Shi):
            prev.append(((s,),(k,)))
            yield (s,),(k,)
        for cnt in range(2,Shi.total()+1):
            nextprev = []
            for ss,ks in prev:
                for s in ss:
                    for e in s:
                        Shi[e] -= 1
                for s,k in self.iter_sets(Shi=Shi,Slolex=ss[-1]):
                    if replacement or s != ss[-1]:
                        nextprev.append( (ss+(s,),ks+(k,)) )
                        yield (ss+(s,),ks+(k,))
                for s in ss:
                    for e in s:
                        Shi[e] += 1
            prev = nextprev
                

# function to enumerate integer points of a polytope defined by input linear program
# (sage MixedIntegerLinearProgram). If xs, a subset of variables, is provided,
# this problem is solved for the projection of the polytope away from the
# variables not in xs (solutions with integer values for the xs and possibly
# noninteger values for the rest, achieving all the possible xs values).
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
        v = math.floor(csol[x]+1e-10)
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

class IntegerProgram:
    def __init__(self,lp,xs=None):
        self.lp = deepcopy(lp)
        self.bounds = [(x, lp.get_min(lp[x]), lp.get_max(lp[x]))
                    for x in lp.default_variable().keys()]
        if xs is None:
            self.xs = list(self.lp.default_variable().keys())
        else:
            self.xs = xs
        self.infeas = SetSystem()
        self.sols = SetSystem()
    def add_infeas(self, bans):
        rem = list(self.infeas.iter_sets(Slo=bans))
        for r,_ in rem:
            self.infeas.remove_set(r)
        self.infeas.add_set(bans)
    def reset_bounds(self):
        for x,lo,hi in self.bounds:
            self.lp.set_min(self.lp[x],lo)
            self.lp.set_max(self.lp[x],hi)
    # if unsafe_no_copy_lp, this generator must either be exhasted before
    # calling any other methods, or reset_bounds() must be called (and no
    # more solutions emitted)
    def extend_psol(self,psol=[],unsafe_no_copy_lp=False):
        from sage.numerical.mip import MIPSolverException
        bans = [(x,i) for x,v in psol for i in range(int(self.lp.get_min(self.lp[x])), v) ]
        bans.extend([(x,i) for x,v in psol for i in range(v+1, int(self.lp.get_max(self.lp[x])+1))])
        if any(True for _ in self.infeas.iter_sets(Shi=bans)):
            return
        if unsafe_no_copy_lp:
            lp = deepcopy(self.lp)
        else:
            lp = self.lp
        for x,v in psol:
            lp.set_min(lp[x],v)
            lp.set_max(lp[x],v)
        st = []
        def dfs(check_solvable=True):
            try:
                if check_solvable:
                    lp.solve()
            except MIPSolverException:
                return
            csol = lp.get_values(lp.default_variable())
            remxs = [x for x in self.xs if lp.get_min(lp[x]) < lp.get_max(lp[x])]
            if len(remxs) == 0:
                sol = [(x,csol[x]) for x in self.xs]
                self.sols.add_set(sol)
                yield sol
                return
            x = max(remxs, key=lambda x: abs(csol[x]-round(csol[x])))
            v = math.floor(csol[x]+1e-10)
            omin = int(lp.get_min(lp[x]))
            omax = int(lp.get_max(lp[x]))
            if v == omax:
                v -= 1
            assert v >= omin and v+1 <= omax
            print('%s%s %d %.2f %d' % (''.join(st),str(x),lp.get_min(lp[x]),
                                     csol[x],lp.get_max(lp[x])))
            def tryhi():
                for i in range(omin,v+1):
                    bans.append((x,i))
                if not any(True for _ in self.infeas.iter_sets(Shi=bans)):
                    lp.set_min(lp[x],v+1)
                    havesol = False
                    for res in dfs(v+1-csol[x] >= 1e-10):
                        havesol = True
                        yield res
                    if not havesol:
                        self.add_infeas(bans)
                    lp.set_min(lp[x],omin)
                for i in range(omin,v+1):
                    bans.pop()
            def trylo():
                for i in range(v+1,omax+1):
                    bans.append((x,i))
                if not any(True for _ in self.infeas.iter_sets(Shi=bans)):
                    lp.set_max(lp[x],v)
                    havesol = False
                    for res in dfs(csol[x]-v >= 1e-10):
                        havesol = True
                        yield res
                    if not havesol:
                        self.add_infeas(bans)
                    lp.set_max(lp[x],omax)
                for i in range(v+1,omax):
                    bans.pop()
            hi_first = csol[x] - v >= 0.5
            st.append(' ')
            if hi_first:
                for sol in tryhi():
                    yield sol
                st[-1] = '.'
            for sol in trylo():
                yield sol
            if not hi_first:
                st[-1] = '.'
                for sol in tryhi():
                    yield sol
            st.pop()
        havesol = False
        for sol in dfs():
            havesol = True
            yield sol
        if not havesol:
            self.add_infeas(bans)
    def can_extend_psol(self,psol=[]):
        if any(True for _ in self.sols.iter_sets(Slo=psol)):
            return True
        res = any(True for _ in self.extend_psol(psol,unsafe_no_copy_lp=True))
        if res:
            self.reset_bounds()
        return res
        

        
