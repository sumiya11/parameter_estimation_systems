using Pkg
Pkg.activate("BB")

using Revise
using ParameterEstimation

name = "sum_of_species_august_21"
ParameterEstimation.MODEL[] = name
ParameterEstimation.OUTPUT_DIR[] = "case_by_case/$name"

include("case_by_case/$(name)/$(name).jl")

using RationalUnivariateRepresentation, Groebner, Nemo, HomotopyContinuation, StatsBase

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

function nemo2hc(expr_tree)
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

##### AA ######

# Load system
sys = include("case_by_case/$name/sys_$(name)_AAA_abstract_algebra.txt")

# Round coefficients to 4 "digits"
sys_orig = sys
sys = map(f -> map_coefficients(c -> rationalize(BigInt, round(BigFloat(c), digits=4)), f), sys)

# Compute RUR
rur, sep = zdim_parameterization(sys, get_separating_element=true);

@info "" length(rur[1])

##### HC ######

hc_sys = nemo2hc.(sys)
mv_ = HomotopyContinuation.mixed_volume(HomotopyContinuation.System(hc_sys))
@info "MV" mv_

println("degrees in vars = ", countmap(map(f -> HomotopyContinuation.ModelKit.degree(f), hc_sys)))
