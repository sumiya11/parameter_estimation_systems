using ModelingToolkit, Catalyst
using SIAN, Nemo, StructuralIdentifiability, Groebner
using LinearAlgebra, OffsetArrays

include("../make_poly_from_ode.jl")
include("../plug_numbers_in_poly.jl")
include("../quotient_basis.jl")

n = 5
t = default_t()

# States

str = join(["S"*string(i)*"(t)" for i in 0:n], " ")
expr = "@species $str"
S = OffsetArray(eval(Meta.parse(expr)), 0:n)

str = join(["SK"*string(i)*"(t)" for i in 0:n], " ")
expr = "@species $str"
SK = OffsetArray(eval(Meta.parse(expr)), 0:n)

str = join(["SF"*string(i)*"(t)" for i in 0:n], " ")
expr = "@species $str"
SF = OffsetArray(eval(Meta.parse(expr)), 0:n)

@species K(t) F(t)

# Parameters 

str = join(["a"*string(i)*"" for i in 0:n+1], " ")
expr = "@parameters $str"
A = OffsetArray(eval(Meta.parse(expr)), 0:n+1)

str = join(["b"*string(i)*"" for i in 0:n+1], " ")
expr = "@parameters $str"
B = OffsetArray(eval(Meta.parse(expr)), 0:n+1)

rxs = vcat(
  [Reaction(231, [S[0], K], [SK[0]])],
  [Reaction(129, [SK[i]], [SK[i+1]]) for i in 0:n-2],
  [Reaction(111, [SK[n-1]], [K, S[n]])],
  [Reaction(229, [SK[0]], [S[0], K])],
  # [Reaction(A[n+1], [S[0], K], [SK[0]])],
  # [Reaction(A[i], [SK[i]], [SK[i+1]]) for i in 0:n-2],
  # [Reaction(A[n-1], [SK[n-1]], [K, S[n]])],
  # [Reaction(A[n], [SK[0]], [S[0], K])],
  
  [Reaction(1 // 3, [SF[n]], [S[n], F])], # [Reaction(B[0], [SF[n]], [S[n], F])],
  [Reaction(1 // (7*i), [SF[i]], [SF[i-1]]) for i in n:-1:2], # [Reaction(B[i], [SF[i]], [SF[i-1]]) for i in n:-1:2],
  [Reaction(1 // 11, [SF[1]], [F, S[0]])], # [Reaction(B[1], [SF[1]], [F, S[0]])],
  [Reaction(1//13, [S[n], F], [SF[n]])] # [Reaction(B[n+1], [S[n], F], [SF[n]])]
)

@named crn = ReactionSystem(rxs, t)
crn = complete(crn)

ode = convert(ODESystem, crn)

@variables Y1(t) Y2(t)
measured_quantities = [
  Y1 ~ S[0],
  Y2 ~ F
]

# str = join(["y"*string(i)*"(t)" for i in 1:length(ModelingToolkit.unknowns(ode))], " ")
# expr = "@variables $str"
# Y = eval(Meta.parse(expr))
# measured_quantities = [
#   Y[i] ~ ModelingToolkit.unknowns(ode)[i]
#   for i in 1:length(ModelingToolkit.unknowns(ode))
# ]

res = ode_to_poly(ode, measured_quantities)
poly = res["polynomial_system"]
states = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(collect(S), collect(SK), collect(SF), F, K)), gens(parent(poly[1])))
parameters = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(collect(A),collect(B))), gens(parent(poly[1])))
states = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in states]
parameters = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in parameters]
parameters = Vector{eltype(states)}(parameters)
poly_num, subs = plug_numbers_in_poly(poly, states, parameters)

@info "" length(poly_num)
@info "Linear equations" length(filter(x -> x in gens(parent(poly_num[1])), map(leading_monomial, poly_num)))

r_drl, _ = polynomial_ring(base_ring(parent(poly_num[1])), symbols(parent(poly_num[1])), internal_ordering=:degrevlex)
poly_num = map(f -> StructuralIdentifiability.parent_ring_change(f, r_drl), poly_num)
poly_num_zp = map(f -> map_coefficients(c -> GF(2^30+3)(c), f), poly_num)
gb = groebner(poly_num_zp);
@info "" length(quotient_basis(gb))
