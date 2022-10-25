function ensure_correct(z, M, N, q, l, u, w; tol=1e-4)
    r = M*z + N*w + q

    I₋ = r .< -tol
    I₊ = r .> tol
    I₀ = isapprox.(r, 0; atol=tol)
    
    @test all(isapprox.(z[I₊] - l[I₊], 0; atol=tol))
    @test all(isapprox.(z[I₋] - u[I₋], 0; atol=tol))
    @test all(l.-tol .≤ z .≤ u.+tol)
end

@testset "lmcp" begin
    rng = MersenneTwister(420)
    for i = 1:10
        M = randn(rng, 5,10)
        M = M'*M
        N = randn(rng, 10, 3)
        w = randn(rng, 3)
        q = randn(rng, 10)
        l = randn(rng, 10)
        u = randn(rng, 10)
        inconsistent = u.<l
        u[inconsistent] = l[inconsistent]
        z = solve_lmcp(M, N, q, l, u, w)
        ensure_correct(z, M, N, q, l, u, w)
    end
end

