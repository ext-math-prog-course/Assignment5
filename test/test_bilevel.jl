@testset "bilevel" begin
    Q = [1 0;
         0 1.0]
    A = zeros(0,2)
    l1 = u1 = zeros(0)

    M = ones(1,1)
    N = -ones(1,1)
    o = zeros(1)
    l2 = zeros(1)
    u2 = fill(Inf, 1)

    x_init = [9.0, 8]

    q = -[1.0; -2]
    x = solve_bilevel(Q, q, A, l1, u1, M, N, o, l2, u2, x_init)
    @test isapprox(x, [0.0, 0]; atol=1e-5)
    
    q = -[3.0; 0]
    x = solve_bilevel(Q, q, A, l1, u1, M, N, o, l2, u2, x_init)
    @test isapprox(x, [1.5, 1.5]; atol=1e-5)
    
    q = -[-5.0; -2]
    x = solve_bilevel(Q, q, A, l1, u1, M, N, o, l2, u2, x_init)
    @test isapprox(x, [-5.0, 0]; atol=1e-5)
    
    q = -[0.0; 0]
    x = solve_bilevel(Q, q, A, l1, u1, M, N, o, l2, u2, x_init)
    @test isapprox(x, [0.0, 0]; atol=1e-5)
end
