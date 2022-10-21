# Assignment5
# Instructions
1. A linear mixed complementarity problem is a generalization of the LCP. This
   formulation is parameterized by a square matrix M, and vectors q, l, and u.
   A vector z is a solution to the LMCP if for all i in 1...n, one of the three
   cases is satisfied: 
    z_i = l_i and (Mz+q)i >= 0.
    l_i <= z_i <= u_i and (Mz+q)i = 0
    z_i = u_i and (Mz+q)i <= 0
   Complete src/solve_lmcp.jl, using the PATH solver to solve LMCPs so that
   tests in test/test_lmcp.jl pass.
2. Given a solution z to a parametric LMCP, write a function to compute the
   polyhedral regions representing the set of all local solutions. See
   src/solve_lmcp.jl for details and definitions.
3. Implement an algorithm for solving bilevel optimization problems. This
   algorithm does not need to handle EPECs, or compute the solution SET of the
   blievel problem. See src/solve_bilevel.jl for details. 
4. Tests passing DOES NOT necessariliy imply full credit for this assignment. It is not
   trivial to write tests for bilevel programs in a way that doesn't provide
   a means of solving them.

# Julia resources

- [Julia Manual](https://docs.julialang.org/en/v1/manual/getting-started/)
- [Julia Package/Environment Guide](https://pkgdocs.julialang.org/v1/)
- [JuMP (JuliaMathematicalProgramming) Documentation](https://jump.dev/JuMP.jl/stable/)
- [Julia workflow tips](https://m3g.github.io/JuliaNotes.jl/stable/workflow/)
- [Julia Academy](https://juliaacademy.com/courses)
- [Algorithms for Optimization book with Julia Examples](https://algorithmsbook.com/optimization/)
