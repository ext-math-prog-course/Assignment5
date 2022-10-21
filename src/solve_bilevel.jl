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
function solve_bilevel(Q, q, A, l1, u1, M, N, o, l2, u2, x_init)
    return x_init # FIX ME
end

