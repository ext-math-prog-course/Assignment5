"""
Solves the program

min         0.5 x'Qx + x'q
x = [w,z]

s.t.        l1 ≤ Ax ≤ u1
            x2 solves LMCP(M, N, q, l2, u2, x1)
           (z)                              (w)

Pseudo code:

x₁ = x_init
(zdim, wdim) = size(N)

for k = 1,2,...
    wₖ = xₖ[1:wdim]
    zₖ = solve_lmcp(M, N, o, l2, u2, wₖ) NOTE THE VARIABLE IS NAMED "o" TO DISAMBIGUIATE WITH "q" in QP
    local_sols = find_local_solution_set(M,N,o,l2,u2,wₖ,zₖ)
    for each sol in local_sols
        form and solve approximate QP:
            x* = argmin  0.5 x'Qx + x'q
                    x   
                    s.t. l1 ≤ Ax ≤ u1
                         x ∈ sol
    if all( x* ≈ [wₖ; zₖ]), then x* is solution (can assume each approximate QP has unique solution)
    else, pick some x* != [wₖ; zₖ] and set xₖ₊₁ = x* and repeat
end
"""
function solve_bilevel(Q, q, A, l1, u1, M, N, o, l2, u2, x_init; max_iters=50)
    x = x_init
    zdim, wdim = size(N)

    for iter = 1:max_iters
        w = x[1:wdim]
        z = solve_lmcp(M, N, o, l2, u2, w)
        local_sols = find_local_solution_set(M, N, o, l2, u2, w, z)
        all_same = true
        for poly in local_sols
            l = [l1; poly.l]
            u = [u1; poly.u]
            A2 = poly.A
            AA = [A; A2[:,end-wdim+1:end] A2[:, 1:zdim]]
            
            m = OSQP.Model()
            OSQP.setup!(m; P=sparse(Q), q, A=sparse(AA), l, u, verbose=false, polish=true, eps_abs=1e-8, eps_rel=1e-8)

            res = OSQP.solve!(m)
            xopt = res.x
            if !isapprox(xopt, x; atol=1e-5)
                all_same = false
                x .= xopt
                break
            end
        end
        if all_same
            return x
        end
    end
    error("Can't find solution")
end

