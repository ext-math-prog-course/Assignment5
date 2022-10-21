

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
    z = similar(q) # TODO fix me 
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
    # TODO Fill these in
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
