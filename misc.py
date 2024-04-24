
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

def get_positive_constraints(points,xs):
    xsix = {x:i for i,x in enumerate(xs)}
    n = len(xsix)
    R = PolynomialRing(QQ,'x',n)
    return frobby.alexander_dual(R.ideal([prod(R.gen(xsix[m]) for m in pt if m in xsix) for pt in points]))
    # return R.ideal(1).intersection(*[R.ideal([R.gen(xsix[m]) for m in pt if m in xsix]) for pt in points])

def add_positive_constraints(lp,xs,constraints):
    assert len(xs) == constraints.ring().ngens()
    for not_all_missing in constraints.gens():
        if not not_all_missing.is_zero():
            lp.add_constraint(lp.sum(lp[xs[i]] for i,_ in not_all_missing.exponents()[0].sparse_iter()) <= not_all_missing.degree()-1)
            
    # Rfull = PolynomialRing(QQ,'x',2*n)
    # bad_condition_full = Rfull.ideal(1)
    # triv = Rfull.ideal([Rfull.gen(2*i)*Rfull.gen(2*i+1) for i in range(n)])
    # for pt in points:
    #     bad_condition = bad_condition.intersection(R.ideal([
    #         R.gen(xsix[m]) for m in sol if m in xsix]))
    #     bad_condition = R.ideal([m for m in bad_condition.gens() if m.degree() <= 20])
    #     print (Counter( map(lambda p: p.degree(), bad_condition.gens()) ))
    #     # bad_condition_full = bad_condition_full.intersection(Rfull.ideal([
    #     #     Rfull.gen(2*xsix[m]) if v==0 else Rfull.gen(2*xsix[m]+1) for m,v in solall.items() if m in xsix]))
    #     # bad_condition_full = Rfull.ideal([triv.reduce(m) for m in bad_condition_full.gens()])
    #     # print (Counter( map(lambda p: p.degree(), bad_condition_full.gens()) ))
    # for bad in bad_condition_full.gens():
    #     lp.add_constraint(lp.sum(lp[xs[i//2]] if i % 2 == 0 else 1-lp[xs[i//2]] 
    #                              for i,_ in bad.exponents()[0].sparse_iter()) <= bad.degree()-1)

class BinaryProgram:
    def __init__(self,lp,xs=None,use_infeas=True,prunef=lambda psol: True):
        self.lp = lp
        self.xs = list(self.lp.default_variable().keys()) if xs is None else xs
        self.use_infeas = use_infeas
        self.infeas = SetSystem()
        self.sols = SetSystem()
        self.psol = set()
        self.prunef = prunef
    def set_min(self, x, v):
        v = int(v)
        self.lp.set_min(self.lp[x], v)
        if v == 1:
            assert (x,1) not in self.psol
            self.psol.add((x,1))
        else:
            self.psol.remove((x,1))
    def set_max(self, x, v):
        v = int(v)
        self.lp.set_max(self.lp[x], v)
        if v == 0:
            assert (x,0) not in self.psol
            self.psol.add((x,0))
        else:
            self.psol.remove((x,0))
    def get_min(self, x):
        return int(self.lp.get_min(self.lp[x]))
    def get_max(self, x):
        return int(self.lp.get_max(self.lp[x]))
    def extend(self):
        sol = list(islice(self.sols.iter_sets(Slo=self.psol),1))
        if len(sol) == 1:
            return sol[0][0]
        def prunef(psol):
            return not any(True for _ in self.infeas.iter_sets(Shi=psol.items())) and self.prunef(psol)
        if not self.use_infeas:
            prunef = self.prunef
        res = ip_to_scipy(self.lp)
        res = [{x: res.x[next(iter(self.lp[x].dict().keys()))] for x in self.lp.default_variable().keys()}] if res.success else []
        # res = list(islice(lp_integer_points(self.lp,self.xs,fullsol=True,prunef=prunef),1))
        if len(res) == 0:
            if self.use_infeas:
                rem = list(self.infeas.iter_sets(Slo=self.psol))
                for r,_ in rem:
                    self.infeas.remove_set(r)
                self.infeas.add_set(self.psol)
        else:
            sol = [(x,int(round(v))) for x,v in res[0].items() if abs(v-round(v)) < 1e-10]
            self.sols.add_set(sol)
            return sol
        
# class BinaryProgram:
#     def __init__(self,lp,xs=None,use_infeas=True,prunef=lambda psol: True):
#         self.lp = deepcopy(lp)
#         self.lp.set_binary(self.lp.default_variable())
#         self.xs = list(self.lp.default_variable().keys()) if xs is None else xs
#         self.use_infeas = use_infeas
#         self.infeas = SetSystem()
#         self.sols = SetSystem()
#         self.psol = set()
#         self.prunef = prunef
#     def set_min(self, x, v):
#         v = int(v)
#         self.lp.set_min(self.lp[x], v)
#         if v == 1:
#             assert (x,1) not in self.psol
#             self.psol.add((x,1))
#         else:
#             self.psol.remove((x,1))
#     def set_max(self, x, v):
#         v = int(v)
#         self.lp.set_max(self.lp[x], v)
#         if v == 0:
#             assert (x,0) not in self.psol
#             self.psol.add((x,0))
#         else:
#             self.psol.remove((x,0))
#     def get_min(self, x):
#         return int(self.lp.get_min(self.lp[x]))
#     def get_max(self, x):
#         return int(self.lp.get_max(self.lp[x]))
#     def extend(self):
#         sol = list(islice(self.sols.iter_sets(Slo=self.psol),1))
#         if len(sol) == 1:
#             return sol[0][0]
#         def prunef(psol):
#             return not any(True for _ in self.infeas.iter_sets(Shi=psol.items())) and self.prunef(psol)
#         if not self.use_infeas:
#             prunef = self.prunef
#         res = list(islice(lp_integer_points(self.lp,self.xs,fullsol=True,prunef=prunef),1))
#         from sage.numerical.mip import MIPSolverException
#         try:
#             self.lp.solve()
#             sol = self.lp.get_values(self.lp.default_variable())
#             sol = [(x,int(round(v))) for x,v in sol.items() if abs(v-round(v)) < 1e-10]
#             self.sols.add_set(sol)
#             return sol
#         except MIPSolverException:
#             if self.use_infeas:
#                 rem = list(self.infeas.iter_sets(Slo=self.psol))
#                 for r,_ in rem:
#                     self.infeas.remove_set(r)
#                 self.infeas.add_set(self.psol)
                   
class BinaryPrograms:
    def __init__(self,lp,xss,xsshi,use_infeas=True):
        ms = dict([(next(iter(v.dict().keys())),m) for m,v in lp.default_variable().items()])
        msix = dict([(m,next(iter(v.dict().keys()))) for m,v in lp.default_variable().items()])
        xssi = {}
        for i,xs in enumerate(xss):
            for x in xs:
                xssi.setdefault(x,[]).append(i)
        xsbpis = [set(xs) for xs in xss]
        xsbp = set(x for xs in xss for x in xs)
        for _,(xis,_),_ in lp.constraints():
            curxs = set([ms[xi] for xi in xis])
            curbpis = reduce(lambda a,b : a|b, 
                             (set(xssi.get(x,[])) for x in curxs), set())
            print(curxs,curbpis)
            for bpi in curbpis:
                xsbpis[bpi] |= curxs
            if len(curbpis) >= 2:
                xsbp |= curxs
        # xssalli = deepcopy(xssi)
        # for i,xs in enumerate(xsshi):
        #     for x in xs:
        #         xssalli.setdefault(x,[]).append(len(xs)+i)
        # for bpi,xs in enumerate(xsbpis):
        #     xsbpis[bpi] = reduce(lambda a,b: a|b, (set(xss[i])|set([x]) if i < len(xss) else set(xsshi[i-len(xss)])|set([x]) 
        #          for x in xs for xsi in xssalli[x]), set())
        # xsbp = reduce(lambda a,b: a|b, (set(xss[i])|set([x]) if i < len(xss) else set(xsshi[i-len(xss)])|set([x]) 
        #      for x in xsbp for xsi in xssalli[x]), set())
        
        self.bps = []
        for xs,ys in zip(xss,xsbpis):
            curxis = set(msix[x] for x in ys)
            lpcur = MixedIntegerLinearProgram(solver='GLPK')
            lpcur.set_min(lpcur.default_variable(),0)
            lpcur.set_max(lpcur.default_variable(),1)
            for lo,(xis,cs),hi in lp.constraints():
                if all(xi in curxis for xi in xis):
                    lpcur.add_constraint(lpcur.sum(c*lpcur[ms[xi]] for xi,c in zip(xis,cs)),min=lo,max=hi)
            self.bps.append(BinaryProgram(lpcur,xs,use_infeas))
        def prunef(psol):
            revert = []
            for x,v in psol.items():
                for i in xssi.get(x,[]):
                    if v == 0:
                        assert self.bps[i].get_min(x) == 0
                        if self.bps[i].get_max(x) == 1:
                            revert.append((i,x,v))
                            self.bps[i].set_max(x,0)
                    else:
                        assert self.bps[i].get_max(x) == 1
                        if self.bps[i].get_min(x) == 0:
                            revert.append((i,x,v))
                            self.bps[i].set_min(x,1)
            res = all(bp.extend() is not None for bp in self.bps)
            for i,x,v in revert:
                if v == 0:
                    self.bps[i].set_max(x,1)
                else:
                    self.bps[i].set_min(x,0)
            return res
        
        lpcur = MixedIntegerLinearProgram(solver='GLPK')
        lpcur.set_min(lpcur.default_variable(),0)
        lpcur.set_max(lpcur.default_variable(),1)
        curxis = set(msix[x] for x in xsbp)
        for lo,(xis,cs),hi in lp.constraints():
            if all(xi in curxis for xi in xis):
                lpcur.add_constraint(lpcur.sum(c*lpcur[ms[xi]] for xi,c in zip(xis,cs)),min=lo,max=hi)
        from operator import concat
        self.bp = BinaryProgram(lpcur,reduce(concat,xss),use_infeas=False,prunef=prunef)
        self.ysbis = {}
        for y in self.bp.lp.default_variable().keys():
            self.ysbis.setdefault(y,[]).append(-1)
        for bpi, bp in enumerate(self.bps):
            for y in bp.lp.default_variable().keys():
                self.ysbis.setdefault(y,[]).append(bpi)
        self.psol = set()
    def set_min(self, x, v):
        v = int(v)
        if v == 1:
            assert (x,1) not in self.psol
            self.psol.add((x,1))
        else:
            self.psol.remove((x,1))
        for bi in self.ysbis[x]:
            if bi == -1:
                self.bp.set_min(x,v)
            else:
                self.bps[bi].set_min(x,v)
    def set_max(self, x, v):
        v = int(v)
        if v == 0:
            assert (x,0) not in self.psol
            self.psol.add((x,0))
        else:
            self.psol.remove((x,0))
        for bi in self.ysbis[x]:
            if bi == -1:
                self.bp.set_max(x,v)
            else:
                self.bps[bi].set_max(x,v)
    def get_min(self, x):
        bi = self.ysbis[x][0]
        return self.bp.get_min(x) if bi == -1 else self.bps[bi].get_min(x)
    def get_max(self, x):
        bi = self.ysbis[x][0]
        return self.bp.get_max(x) if bi == -1 else self.bps[bi].get_max(x)
    def extend(self):
        return self.bp.extend()
