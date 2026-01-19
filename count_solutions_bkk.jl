using Pkg
# Pkg.activate("env")

using Revise
using Nemo
using HomotopyContinuation
using StructuralIdentifiability
using ParameterEstimation
using StatsBase

examples = [
  # "case_by_case/akt_pathway_small/sys_akt_pathway_small_AAA_abstract_algebra.txt",
  # "case_by_case/paper-ex/sys_paper-ex_AAA_abstract_algebra.txt",
  # "case_by_case/goodwin/sys_goodwin_AAA_abstract_algebra.txt",
  # "case_by_case/treatment/sys_treatment_AAA_abstract_algebra.txt",
  # "case_by_case/pk1/sys_PK1_AAA_abstract_algebra.txt",
  "case_by_case/crn/sys_crn_AAA_abstract_algebra.txt",
  # "case_by_case/seir_36/sys_seir_36_AAA_abstract_algebra.txt",
  # "case_by_case/nfkb_lilliput/sys_nfkb_lilliput_AAA_abstract_algebra.txt",
  # "case_by_case/akt_pathway_small/sys_akt_pathway_small_AAA_abstract_algebra.txt",
  # # "case_by_case/lotka-volterra/sys_lotka-volterra_AAA_abstract_algebra.txt",
  # "case_by_case/crauste/sys_crauste_AAA_abstract_algebra.txt",
  # "case_by_case/crauste_two_squared/sys_crauste_two_squared_AAA_abstract_algebra.txt",
  # "case_by_case/crauste_four_squared/sys_crauste_four_squared_AAA_abstract_algebra.txt",
]

function nemo2hc(expr_tree::Union{Expr, Symbol})
	#traverse expr_tree
	if typeof(expr_tree) == Symbol
		return HomotopyContinuation.Expression(HomotopyContinuation.variables(expr_tree)[1])
	end
	if typeof(expr_tree) == Expr
		if expr_tree.head == :macrocall
			return eval(expr_tree)
		end
		if expr_tree.head == :call
			if expr_tree.args[1] in [:+, :-, :*, :/, :^, ://]
				if length(expr_tree.args) == 2
					return eval(expr_tree.args[1])(nemo2hc(expr_tree.args[2]))
				else
					# println(expr_tree.args[2:end])
					return reduce(eval(expr_tree.args[1]),
						map(nemo2hc, expr_tree.args[2:end]))
				end
			end
		end
	end
end

function nemo2hc(expr_tree::QQMPolyRingElem)
	# println(expr_tree)
	return nemo2hc(Meta.parse(string(expr_tree)))
end

function nemo2hc(expr_tree::Number)
	return expr_tree
end

function nemo2hc(expr_tree::Nemo.Generic.FracFieldElem)
	numer, denom = Nemo.numerator(expr_tree), Nemo.denominator(expr_tree)
	return nemo2hc(numer) / nemo2hc(denom)
end

hc_sys = nothing
if true
for example in examples
  println()
  @info "Example: $example"
  include(example)
  hc_vars = [HomotopyContinuation.ModelKit.variables(string(s))[1] for s in gens(R)]
  global hc_sys = nemo2hc.(sys)
  main_vars = filter(x -> !endswith(string(x), r"_[1-9]"), variables(hc_sys))
  diff_vars = setdiff(variables(hc_sys), main_vars)
  println("main_vars : $(length(main_vars))")
  println("degrees in main vars = ", countmap(map(f -> HomotopyContinuation.ModelKit.degree(f, main_vars), hc_sys)))
  println("degrees in vars = ", countmap(map(f -> HomotopyContinuation.ModelKit.degree(f), hc_sys)))
  println(hc_sys)
  mv = HomotopyContinuation.mixed_volume(System(hc_sys))
  println("MV = $mv")
end
end

#=
[ Info: Example: case_by_case/paper-ex/sys_paper-ex_AAA_abstract_algebra.txt
main_vars = Variable[mu_0, x_0]
degrees in main vars = [2, 2, 1, 1, 1]
MV = 9
[ Info: Example: case_by_case/goodwin/sys_goodwin_AAA_abstract_algebra.txt
main_vars = Variable[Ki_0, k1_0, k2_0, k4_0, k5_0, k6_0, x1_0, x2_0, x3_0]
degrees in main vars = [1, 12, 1, 2, 0, 11, 0, 1, 2, 0, 11, 0, 1, 1, 0, 11, 0, 1, 1, 0, 1, 1]
MV = 40
[ Info: Example: case_by_case/treatment/sys_treatment_AAA_abstract_algebra.txt
main_vars = Variable[In_0, N_0, S_0, Tr_0, b_0, d_0, g_0, nu_0]
degrees in main vars = [1, 2, 1, 0, 0, 1, 4, 0, 1, 3, 4, 0, 1, 3, 0, 3, 0, 1, 3, 3, 0, 0, 1, 3, 3, 0, 0, 1, 3, 0, 3]
MV = 0
[ Info: Example: case_by_case/pk1/sys_PK1_AAA_abstract_algebra.txt
main_vars = Variable[k2_0, k3_0, k4_0, k5_0, k6_0, k7_0, s3_0, u1_0, x1_0, x2_0, x3_0, x4_0]
degrees in main vars = [1, 2, 2, 2, 0, 1, 2, 2, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1]
MV = 12

[ Info: Example: case_by_case/crn/sys_crn_AAA_abstract_algebra.txt
main_vars = Variable[k1_0, k2_0, k3_0, k4_0, k5_0, k6_0, x1_0, x2_0, x3_0, x4_0, x5_0, x6_0]
degrees in main vars = [1, 3, 1, 3, 0, 2, 3, 3, 0, 2, 3, 3, 0, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2, 0, 2, 2, 2]
Mixed volume:  1625    Time: 0:00:06
MV = 1625

[ Info: Example: case_by_case/seir_36/sys_seir_36_AAA_abstract_algebra.txt
main_vars = Variable[De_0, Di_0, E_0, F_0, I_0, N_0, S_0, beta_0, beta_d_0, gamma_0, gamma_d_0, mu_0_0, mu_d_0, mu_i_0, nu_0, phi_0, phi_e_0, q_0, s_0, s_d_0]
degrees in main vars = [1, 0, 1, 2, 1, 2, 1, 2, 1, 0, 1, 0, 0, 1, 2, 0, 1, 4, 0, 1, 0, 1, 1, 0, 1, 3, 4, 0, 1, 0, 1, 1, 0, 1, 3, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 3, 0, 0, 3, 0, 0, 1, 1, 1, 0, 1, 3, 3, 0, 0, 0, 0, 1, 3, 0, 3, 0, 0]
MV = 0


[ Info: Example: case_by_case/akt_pathway_small/sys_akt_pathway_small_AAA_abstract_algebra.txt
main_vars = Variable[Akt_0, EGF_EGFR_0, S6_0, pAkt_0, pAkt_S6_0, pEGFR_0, pEGFR_Akt_0, pS6_0, reaction_2_k2_0, reaction_3_k1_0, reaction_4_k1_0, reaction_5_k1_0, reaction_6_k1_0, reaction_7_k1_0, reaction_8_k1_0, reaction_9_k1_0]
degrees in main vars = [1, 2, 2, 1, 3, 3, 1, 2, 0, 1, 1, 2, 2, 0, 2, 2, 3, 0, 1, 0, 1, 1, 1, 1, 0, 2, 2, 2, 0, 1, 0, 1, 1, 1, 1, 0, 2, 2, 2, 0, 1, 0, 1, 1, 1, 1, 0, 2, 2, 2, 0, 1, 0, 1, 1, 1, 1]
Mixed volume:  8095    Time: 0:02:42
MV = 8095
Sols: 110


[ Info: Example: case_by_case/crauste/sys_crauste_AAA_abstract_algebra.txt
main_vars = Variable[E_0, M_0, N_0, P_0, S_0, delta_EL_0, delta_LM_0, delta_NE_0, mu_EE_0, mu_LE_0, mu_LL_0, mu_M_0, mu_N_0, mu_PE_0, mu_PL_0, mu_P_0, rho_E_0, rho_P_0]
degrees in main vars = [1, 3, 1, 3, 1, 3, 1, 2, 3, 0, 2, 0, 2, 0, 2, 0, 2, 1, 0, 2, 0, 2, 0, 2, 0, 2, 1, 0, 2, 0, 2, 0, 2, 1, 0, 2, 0, 1, 2, 0, 2, 2, 2]
Mixed volume:  53    Time: 0:00:00
MV = 53


[ Info: Example: case_by_case/crauste_two_squared/sys_crauste_two_squared_AAA_abstract_algebra.txt
main_vars = Variable[E_0, M_0, N_0, P_0, S_0, delta_EL_0, delta_LM_0, delta_NE_0, mu_EE_0, mu_LE_0, mu_LL_0, mu_M_0, mu_N_0, mu_PE_0, mu_PL_0, mu_P_0, rho_E_0, rho_P_0]
degrees in main vars = [1, 3, 1, 4, 1, 3, 1, 2, 3, 0, 2, 0, 3, 0, 2, 0, 2, 1, 0, 2, 0, 3, 0, 2, 0, 2, 1, 0, 2, 0, 2, 0, 2, 1, 0, 2, 0, 1, 2, 0, 2, 2, 3]
Mixed volume:  212    Time: 0:00:00
MV = 212
[ Info: Example: case_by_case/crauste_four_squared/sys_crauste_four_squared_AAA_abstract_algebra.txt
main_vars = Variable[E_0, M_0, N_0, P_0, S_0, delta_EL_0, delta_LM_0, delta_NE_0, mu_EE_0, mu_LE_0, mu_LL_0, mu_M_0, mu_N_0, mu_PE_0, mu_PL_0, mu_P_0, rho_E_0, rho_P_0]
degrees in main vars = [1, 4, 1, 4, 1, 3, 1, 2, 3, 0, 3, 0, 3, 0, 2, 0, 2, 1, 0, 3, 0, 3, 0, 2, 0, 2, 1, 0, 3, 0, 2, 0, 2, 1, 0, 2, 0, 1, 2, 0, 2, 3, 3]
Mixed volume:  848    Time: 0:00:00
MV = 848
=#

