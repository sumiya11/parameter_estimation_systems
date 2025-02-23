using RationalUnivariateRepresentation, RS, AbstractAlgebra

# Load system
sys = include("case_by_case/nfkb_lilliput/sys_nfkb_lilliput_AAA_abstract_algebra.txt")

# Round coefficients to 4 "digits"
sys_orig = sys
sys = map(f -> map_coefficients(c -> rationalize(BigInt, round(BigFloat(c), digits=4)), f), sys)

# Compute RUR
rur, sep = zdim_parameterization(sys, get_separating_element=true);

# Isolate
sol_prec = Int32(100)
sol = RS.rs_isolate(rur, sep, output_precision=sol_prec)

# True values:
# i1 = 0.200, i1a = 0.300, k_prod = 0.400,  t1 = 0.500,  t2 = 0.600
p_true = [0.2, 0.3, 0.4, 0.5, 0.6];

# Filter based on p[1] = i1 = 0.2.
sol_best = filter(x -> abs(x[1][1] - p_true[1]) < 0.01, sol)[1]

# Relative error in %
@info "Error = $(round(Float64(maximum((abs.(RS.mid_point.(sol_best[1:5]) .- p_true)) ./ p_true) * 100), digits=2)) %"

# Newton refinement
sol_refined, newton_status = RS.rs_newton(
	sys_orig, symbols(R), RS.mid_point.(sol_best),
	system_precision=Int32(62),
	point_precision=sol_prec,
	output_precision=Int32(1000))

# Relative error in %
@info "Error = $(round(Float64(maximum((abs.(RS.mid_point.(sol_refined[1][1:5]) .- p_true)) ./ p_true) * 100), digits=2)) %"

