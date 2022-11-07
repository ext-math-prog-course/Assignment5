

"""
Finds z which solves the parametric LMCP:
Let n = length(q)

r = Mz + Nw + q
for all i ∈ 1:n, one of the three holds:
    zᵢ = lᵢ      && rᵢ ≥ 0
    lᵢ ≤ zᵢ ≤ uᵢ && rᵢ = 0
    zᵢ = uᵢ      && rᵢ ≤ 0
Use the PATH solver to find z.
Instructions for setting up and using can be
found here:
https://github.com/chkwon/PATHSolver.jl

In particular, take a look at this function in the PATHSolver.jl API:
https://github.com/chkwon/PATHSolver.jl/blob/2087cc0669fa9e1a6faf994bd6b942fadb324a40/src/C_API.jl#L702
"""
function solve_lmcp(M, N, q, l, u, w)
    z0 = zeros(length(q))
    Msp = SparseMatrixCSC{Float64, Int32}(M)
    (path_status, z, info) =  PATHSolver.solve_mcp(Msp, N*w+q, l, u, z0, silent=true)
    return z
end


"""
Identify the local polyhedral regions representing the set of w,z
so that z solves LMCP(M,N,q,l,u,w).

For a particular solution z and parameter w,

Define the index sets J:

J[1] = {i : lᵢ = zᵢ     , (Mz+Nw+q)ᵢ > 0 }
J[2] = {i : lᵢ = zᵢ     , (Mz+Nw+q)ᵢ = 0 }
J[3] = {i : lᵢ < zᵢ < uᵢ, (Mz+Nw+q)ᵢ = 0 }
J[4] = {i :      zᵢ = uᵢ, (Mz+Nw+q)ᵢ = 0 }
J[5] = {i :      zᵢ = uᵢ, (Mz+Nw+q)ᵢ < 0 }
J[6] = {i : lᵢ = zᵢ = uᵢ, (Mz+Nw+q)ᵢ = 0 }


all_local_polys = []
for each S2 in powerset(J[2])
    C2 = all elements in J[2] not in S2
    for each S4 in powerset(J[4])
        C4 = all elements in J[4] not in S4
    
        Now define new sets:
        K[1] = J[1] ∪ C2
        K[2] = J[3] ∪ S2 ∪ S4
        K[3] = J[5] ∪ C4
        K[4] = J[6]

        append Set(K) to all_local_polys
    end
end

Set(K) is defined as the polyhedron:
∀ i ∈ K[1] : (Mz+Nw+q)ᵢ ≥ 0, zᵢ = lᵢ
∀ i ∈ K[2] : (Mz+Nw+q)ᵢ = 0, lᵢ ≤ zᵢ ≤ uᵢ
∀ i ∈ K[3] : (Mz+Nw+q)ᵢ ≤ 0, zᵢ = uᵢ
∀ i ∈ K[4] : -Inf ≤ (Mz+Nw+q)ᵢ ≤ Inf, lᵢ = zᵢ = uᵢ
"""
function find_local_solution_set(M,N,q,l,u,w,z)
    all_local_polys = Vector{Poly}()
     
    J = comp_indices(M, N, q, l, u, z, w)
    J2 = Set(J[2])
    J4 = Set(J[4])
     
    Ks = map(Iterators.product(powerset(J[2]), powerset(J[4]))) do (S2, S4)
            C2 = setdiff(J2, Set(S2)) |> collect
            C4 = setdiff(J4, Set(S4)) |> collect
            K = Dict(1=>Set([J[1];C2]), 2=>Set([J[3];S2;S4]),3=>Set([J[5];C4]),4=>Set(J[6]))
            #local_piece(avi,n,m,K)
    end
    pieces = (local_piece(M, N, q, l, u, K) for K in Ks)
end


function comp_indices(M, N, q, l, u, z, w; tol=1e-4)
    J = Dict{Int, Vector{Int}}()
    r = M*z+N*w+q
    equal_bounds = isapprox.(l, u; atol=tol)
    riszero = isapprox.(r, 0; atol=tol)
    J[1] = findall( isapprox.(z, l; atol=tol) .&& r .> tol )
    J[2] = findall( isapprox.(z, l; atol=tol) .&& riszero .&& .!equal_bounds)
    J[3] = findall( (l.+tol .< z .< u.-tol) .&& riszero )
    J[4] = findall( isapprox.(z, u; atol=tol) .&& riszero .&& .!equal_bounds)
    J[5] = findall( isapprox.(z, u; atol=tol) .&& r .< -tol )
    J[6] = findall( equal_bounds .&& riszero )
    return J
end

function local_piece(M, N, q, l, u, K)
    n, m = size(N)
    A = [M N;
         I(n) spzeros(n,m)]
    bounds = mapreduce(vcat, 1:n) do i
        if i ∈ K[1]
            [-q[i] Inf l[i] l[i]]
        elseif i ∈ K[2]
            [-q[i] -q[i] l[i] u[i]] 
        elseif i ∈ K[3]
            [-Inf -q[i] u[i] u[i]]
        else
            [-Inf Inf l[i] u[i]]
        end
    end
    ll = [bounds[:,1]; bounds[:,3]]
    uu = [bounds[:,2]; bounds[:,4]]
    Poly(A, ll, uu)
end


"""
A polyhedron object:
l ≤ Ax ≤ u
"""
struct Poly
    A::Matrix{Float64}
    l::Vector{Float64}
    u::Vector{Float64}
end

"""
Can be used like 

if x in P
    # do something
end

or like

if x ∈ P
    # do something
end
"""
function Base.in(x, P::Poly; tol=1e-5)
    r = P.A*x
    all((P.l .- tol) .≤ r .≤ (P.u .+ tol))
end
