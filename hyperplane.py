# in this file we provide routines to approximately and exactly solve the problem:
# given vectors v in R^d, find vector w so that (v,w) > 0 for as many as possible

# For all routines here, the vectors v are given as the columns of a matrix M
    
def best_hyperplanes(M):
    from sage.numerical.mip import MIPSolverException

    def get_lp(include, exclude):
        lp = MixedIntegerLinearProgram(maximization=True)
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
    def dfs():
        nonlocal best_so_far
        print(' '*(len(include)+len(exclude)), len(include), best_so_far)
        if len(include) + len(rest) < best_so_far:
            return
        lp = get_lp(include,exclude)
        try:
            lp.solve()
        except MIPSolverException:
            return
        if len(include) > best_so_far:
            best_so_far = len(include)
            sols.clear()
            sols.append(sorted(include))
        elif len(rest) == 0:
            sols.append(sorted(include))
        if len(rest) == 0:
            return
        tcur = lp.get_values(lp.default_variable())
        tcur = vector([tcur[k] for k in range(M.nrows())])
        i = max(rest,key=lambda j: tcur.dot_product(M.column(j)))
        rest.remove(i)
        include.append(i)
        dfs()
        include.pop()
        exclude.append(i)
        dfs()
        exclude.pop()
        rest.add(i)

    t0 = best_hyperplane_approx(M)
    sol = [i for i in range(M.ncols()) if M.column(i).dot_product(t0) > 0]
    sols.append(sol)
    best_so_far = len(sol)
    
    istart = max(rest,key=lambda j: t0.dot_product(M.column(j)))
        
    try:
        rest.remove(istart)
        include.append(istart)
        dfs()
        include.pop()
        exclude.append(istart)
        dfs()
        exclude.pop()
    except KeyboardInterrupt:
        pass
    t0s = []
    for sol in sols:
        lp = get_lp(sol,[])
        lp.solve()
        sol = lp.get_values(lp.default_variable())
        sol = [int(round(sol[i])) for i in range(M.nrows())]
        t0s.append(sol)
    print ('bsf',best_so_far)
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

