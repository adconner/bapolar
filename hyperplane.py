# in this file we provide routines to approximately and exactly solve the problem:
# given vectors v in R^d, find vector w so that (v,w) > 0 for as many as possible

# For all routines here, the vectors v are given as the columns of a matrix M
    
def best_hyperplanes(M):
    from sage.numerical.mip import MIPSolverException

    def get_lp(include, exclude, integer=False):
        lp = MixedIntegerLinearProgram(maximization=True)
        if integer:
            lp.set_integer(lp.default_variable())
        for vi in include:
            lp.add_constraint(lp.sum(M[i,vi]*lp[i] for i in range(M.nrows())) >= 1)
        for vi in exclude:
            lp.add_constraint(lp.sum(M[i,vi]*lp[i] for i in range(M.nrows())) <= 0)
        return lp
    
    include = []
    exclude = []
    rest = set(range(M.ncols()))
    best_so_far = 0
    sols = []
    def dfs(i):
        nonlocal best_so_far, sols
        print(' '*(len(include)+len(exclude)), len(include), best_so_far)
        if len(include) + len(rest) < best_so_far:
            return
        if len(include) > best_so_far:
            best_so_far = len(include)
            sols = [sorted(include)]
        elif len(rest) == 0:
            sols.append(sorted(include))
        include.append(i)
        rest.remove(i)
        lp = get_lp(include,exclude)
        try:
            lp.solve()
            tcur = lp.get_values(lp.default_variable())
            tcur = vector([tcur[k] for k in range(M.nrows())])
            j = max(rest,key=lambda j: tcur.dot_product(M.column(j)))
            dfs(j)
        except MIPSolverException:
            pass
        include.pop()
        exclude.append(i)
        try:
            lp.solve()
            tcur = lp.get_values(lp.default_variable())
            tcur = vector([tcur[k] for k in range(M.nrows())])
            j = max(rest,key=lambda j: tcur.dot_product(M.column(j)))
            dfs(j)
        except MIPSolverException:
            pass
        exclude.pop()
        rest.add(i)

    t0 = best_hyperplane_approx(M)
    sol = [i for i in range(M.ncols()) if M.column(i).dot_product(t0) > 0]
    sols.append(sol)
    best_so_far = len(sol)
    
    istart = max(rest,key=lambda j: t0.dot_product(M.column(j)))
        
    try:
        dfs(istart) # could choose smarter start
    except KeyboardInterrupt:
        pass
    t0s = []
    for sol in sols:
        lp = get_lp(sol,[],integer=True)
        lp.solve()
        sol = lp.get_values(lp.default_variable())
        sol = [int(sol[i]) for i in range(M.nrows())]
        t0s.append(sol)
    return t0s
    
# This solves the problem approximately but fast
def best_hyperplane_approx(M,iters = 30,lambda_final=0.4):
    vs = M.columns()
    vs = [v for v in vs if not v.is_zero()]
    d = len(vs[0])
    vsn = [vv/vv.norm() for v in vs for vv in [v.change_ring(RDF)]]
    wcur = sum(vs).change_ring(RDF)
    if wcur.is_zero():
        print('no initial guess')
        wcur = random_vector(RDF,d)
    wcur /= wcur.norm() 
    lam = lambda it: ((it+1)/iters)*lambda_final+(1-(it+1)/iters)*1.0
    wbest = wcur
    costf = lambda w: len([v for v in vs if v.dot_product(w) <= 0])
    for it in range(iters):
        curlam = lam(it)
        # maximize sum_i (wnext,v_i)^lam ~ sum_i (wcur,v_i)^(lam-1) * (wnext,v_i)
        posvs = [v for v in vsn if abs(v.dot_product(wcur)) > 1e-5]
        if len(posvs) == 0:
            continue
        wnext = sum(v * abs(v.dot_product(wcur))**(curlam-1) for v in posvs)
        if wnext.norm() < 1e-10:
            continue
        wnext /= wnext.norm()
        print (it,curlam,costf(wcur))
        if costf(wcur) < costf(wbest):
            wbest = wcur
        wcur = wnext
    print ('using',costf(wbest))
    return wbest

