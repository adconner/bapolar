
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
