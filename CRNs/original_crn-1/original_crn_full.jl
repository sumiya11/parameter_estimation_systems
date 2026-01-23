# Eq. 3.4 in https://arxiv.org/pdf/1705.10913

using ModelingToolkit, Catalyst
using SIAN, Nemo, StructuralIdentifiability, Groebner
using LinearAlgebra, OffsetArrays

include("../make_poly_from_ode.jl")
include("../plug_numbers_in_poly.jl")
include("../quotient_basis.jl")

t = default_t()

# crn = @reaction_network begin
#     (k1, k2), S0 + K <--> S0K
#     k3, S0K --> Sn + K
#     (k4, k5), Sn + F <--> SnF
#     k6, SnF --> S0 + F
# end

crn = @reaction_network begin
    (1//2, 1//3), S0 + K <--> S0K
    1//17, S0K --> Sn + K
    (1//11, 1//7), Sn + F <--> SnF
    1//5, SnF --> S0 + F
end

ode = convert(ODESystem, crn)

@variables Y1(t) Y2(t)
measured_quantities = [
  Y1 ~ S0,
  Y2 ~ F
]

res = ode_to_poly(ode, measured_quantities)
poly = res["polynomial_system"]
states = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), ModelingToolkit.unknowns(ode)), gens(parent(poly[1])))
parameters = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), ModelingToolkit.parameters(ode)), gens(parent(poly[1])))
states = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in states]
parameters = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in parameters]
parameters = Vector{eltype(states)}(parameters)
poly_num, subs = plug_numbers_in_poly(poly, states, parameters)

@info "" length(poly_num)
@info "Linear equations" length(filter(x -> x in gens(parent(poly_num[1])), map(leading_monomial, poly_num)))

r_drl, _ = polynomial_ring(base_ring(parent(poly_num[1])), symbols(parent(poly_num[1])), internal_ordering=:degrevlex)
poly_num = map(f -> StructuralIdentifiability.parent_ring_change(f, r_drl), poly_num)

gb = groebner(poly_num);
@info "" length(quotient_basis(gb))

# @assert length(quotient_basis(gb)) == 7
