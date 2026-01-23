using ModelingToolkit, Catalyst
using SIAN, Nemo, StructuralIdentifiability, Groebner
using LinearAlgebra, OffsetArrays

include("../make_poly_from_ode.jl")
include("../plug_numbers_in_poly.jl")
include("../quotient_basis.jl")

n = 2
t = default_t()

# states

str = join(["S"*string(i)*"_0(t)" for i in 1:n], " ")
expr = "@species $str"
S0 = eval(Meta.parse(expr))

str = join(["S"*string(i)*"_1(t)" for i in 0:n], " ")
expr = "@species $str"
S1 = OffsetArrays.Origin(0)(eval(Meta.parse(expr)))

str = join(["F"*string(i)*"(t)" for i in 1:n], " ")
expr = "@species $str"
F = eval(Meta.parse(expr))

str = join(["Y"*string(i)*"_0(t)" for i in 1:n], " ")
expr = "@species $str"
Y0 = eval(Meta.parse(expr))

str = join(["Y"*string(i)*"_1(t)" for i in 1:n], " ")
expr = "@species $str"
Y1 = eval(Meta.parse(expr))

# params 

str = join(["a"*string(i)*"_0" for i in 1:n], " ")
expr = "@parameters $str"
A0 = eval(Meta.parse(expr))

str = join(["b"*string(i)*"_0" for i in 1:n], " ")
expr = "@parameters $str"
B0 = eval(Meta.parse(expr))

str = join(["c"*string(i)*"_0" for i in 1:n], " ")
expr = "@parameters $str"
C0 = eval(Meta.parse(expr))

str = join(["a"*string(i)*"_1" for i in 1:n], " ")
expr = "@parameters $str"
A1 = eval(Meta.parse(expr))

str = join(["b"*string(i)*"_1" for i in 1:n], " ")
expr = "@parameters $str"
B1 = eval(Meta.parse(expr))

str = join(["c"*string(i)*"_1" for i in 1:n], " ")
expr = "@parameters $str"
C1 = eval(Meta.parse(expr))

rxs = vcat(
  [Reaction(A0[i], [S1[i-1], S0[i]], [Y0[i]]) for i in 1:n],
  [Reaction(B0[i], [Y0[i]], [S1[i-1], S0[i]]) for i in 1:n],
  [Reaction(1 // (7i), [Y0[i]], [S1[i-1], S1[i]]) for i in 1:n],

  [Reaction(1 // (4i), [F[i], S1[i]], [Y1[i]])    for i in 1:n],
  [Reaction(1 // (2i), [Y1[i]], [F[i], S1[i]])    for i in 1:n],
  [Reaction(1 // (3i), [Y1[i]], [F[i], S0[i]])    for i in 1:n],
)

@named crn = ReactionSystem(rxs, t)
crn = complete(crn)
ode = convert(ODESystem, crn)

str = join(["y"*string(i)*"(t)" for i in 1:length(ModelingToolkit.unknowns(ode))], " ")
expr = "@variables $str"
Y = eval(Meta.parse(expr))
measured_quantities = [
  Y[i] ~ ModelingToolkit.unknowns(ode)[i]
  for i in 1:length(ModelingToolkit.unknowns(ode))
]

res = ode_to_poly(ode, measured_quantities)
poly = res["polynomial_system"]
states = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(collect(S0),collect(S1),collect(F),collect(Y0),collect(Y1))), res["vars"])
parameters = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(collect(A0),collect(B0),collect(C0),collect(A1),collect(B1),collect(C1))), res["vars"])
states = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in states]
parameters = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in parameters]
poly_num, subs = plug_numbers_in_poly(poly, states, parameters)

@info "" length(poly_num)
@info "Linear equations" length(filter(x -> x in gens(parent(poly_num[1])), map(leading_monomial, poly_num)))

# gb = groebner(poly_num);
# @assert length(quotient_basis(gb)) == 1

