using Groebner, Nemo

function eliminate(sys, vars)
    @assert !isempty(sys) && !isempty(vars) && allunique(vars)
    @assert all(x -> parent(x) == parent(sys[1]), vars)
    all_vars = gens(parent(sys[1]))
    vars_block_1 = vars
    vars_block_2 = setdiff(all_vars, vars)
    ordering = DegRevLex(vars_block_1) * DegRevLex(vars_block_2)
    gb = groebner(sys, ordering=ordering)
    gen = filter(poly -> all(x -> degree(poly, x) == 0, vars_block_1), gb)
    gen
end

_r, (y_0, y_1, y_2) = polynomial_ring(Nemo.QQ, ["y_0", "y_1", "y_2"])
R, (x_0, x_1, x_2, mu_0) = polynomial_ring(fraction_field(_r), ["x_0", "x_1", "x_2", "mu_0"])
# R, (x_0, x_1, mu_0) = polynomial_ring(fraction_field(_r), ["x_0", "x_1", "mu_0"])

sys = [
  -x_0^2 - x_0 + y_0,
 mu_0*x_0 + x_1,
 -2*x_0*x_1 - x_1 + y_1,
 mu_0*x_1 + x_2,
 # -2*x_1^2 - 2*x_0*x_2 - x_2 + y_2,
]
el = eliminate(sys, [x_1, x_2])
# two solutions

include("case_by_case/crauste/sys_crauste_AAA_abstract_algebra.txt")
main_vars = filter(x -> !endswith(string(x), r"_[1-9]"), union(map(vars,sys)...))
diff_vars = setdiff(union(map(vars,sys)...), main_vars)

el = eliminate(sys, diff_vars)

function solve_one_by_one(sys)
           SOL = Dict()
           sys_small = sys
           seq = []
           while true
               sys_small = filter(!iszero, sys_small)
               isempty(sys_small) && (@info "All solved!"; break)
               trivial = findall(f -> total_degree(f) == 1 && length(vars(f)) == 1, sys_small)
               length(trivial) == 0 && (@info "No more trivial equations!"; break)
               _vars = map(only ∘ vars, sys_small[trivial])
               sols = map(f -> - trailing_coefficient(f) / leading_coefficient(f), sys_small[trivial])
               varmap = Dict(_vars .=> sols)
               _vars = [v for v in collect(keys(varmap))]
               sols = [varmap[v] for v in _vars]
               @info "Substitute $_vars => $sols"
               sys_small = map(f -> evaluate(f, _vars, sols), sys_small)
               for (var,sol) in varmap
                   @assert !haskey(SOL, var)
                   SOL[var] = sol
               end
                append!(seq, _vars)
           end
           sys_small, seq, SOL
       end

