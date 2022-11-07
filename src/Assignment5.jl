module Assignment5

using PATHSolver
# TODO add this license string
PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")

using LinearAlgebra
using OSQP
using SparseArrays
using Combinatorics

include("solve_lmcp.jl")
include("solve_bilevel.jl")

export solve_lmcp, solve_bilevel

end # module Assignment5
