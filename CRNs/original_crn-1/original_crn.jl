using ModelingToolkit, Catalyst
using SIAN, Nemo, StructuralIdentifiability, Groebner
using LinearAlgebra, OffsetArrays

include("../make_poly_from_ode.jl")
include("../plug_numbers_in_poly.jl")
include("../quotient_basis.jl")

@parameters k1 k2 k3 k4 k5 k6
@variables t x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) y1(t) y2(t)
D = Differential(t)
states = [x1, x2, x3, x4, x5, x6]
parameters = [k1, k2, k3, k4, k5, k6]

@info "CRN: $(length(states)) states, $(length(parameters)) parameters"

@named ode = ODESystem([
                             D(x1) ~ -k1 * x1 * x2 + k2 * x4 + k4 * x6,
                             D(x2) ~ -k1 * x1 * x2 + k2 * x4 + k3 * x4,
                             D(x3) ~ k3 * x4 + k5 * x6 - k6 * x3 * x5,
                             D(x4) ~ k1 * x1 * x2 - k2 * x4 - k3 * x4,
                             D(x5) ~ k4 * x6 + k5 * x6 - k6 * x3 * x5,
                             D(x6) ~ -k4 * x6 - k5 * x6 + k6 * x3 * x5,
                         ], t, states, parameters)
measured_quantities = [y1 ~ x3, y2 ~ x2]

res = ode_to_poly(ode, measured_quantities)
poly = res["polynomial_system"]
states = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), states), gens(parent(poly[1])))
parameters = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), parameters), gens(parent(poly[1])))
states = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in states]
parameters = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in parameters]
poly_num, subs = plug_numbers_in_poly(poly, states, parameters)
r_drl, _ = polynomial_ring(base_ring(parent(poly_num[1])), symbols(parent(poly_num[1])), internal_ordering=:degrevlex)
poly_num = map(f -> StructuralIdentifiability.parent_ring_change(f, r_drl), poly_num)

gb = groebner(poly_num);
@assert length(quotient_basis(gb)) == 7
