using ModelingToolkit, Catalyst, ReactionNetworkImporters
using SIAN, Nemo, StructuralIdentifiability, Groebner
using LinearAlgebra, OffsetArrays

include("../make_poly_from_ode.jl")
include("../plug_numbers_in_poly.jl")
include("../quotient_basis.jl")

n = 5
t = default_t()

# network from basic stoichiometry using ReactionNetworkImporters
@parameters k1 k2 k3 k4 k5
@species A(t) B(t) C(t)
species = [A, B, C]
pars = [k1, k2, k3, k4, k5]
substoich = [2 0 1 0 0;
             0 1 1 0 0;
             0 0 0 1 3]
prodstoich = [0 2 0 1 3;
              1 0 0 1 0;
              0 0 1 0 0]
mn = MatrixNetwork(ones(Int, length(pars)), substoich, prodstoich; species = species,
    params = pars) # a matrix network
prn = loadrxnetwork(mn; name = :testnetwork) # dense version

crn = complete(prn.rn)

ode = convert(ODESystem, crn)

n_observed = 2
str = join(["y"*string(i)*"(t)" for i in 1:length(ModelingToolkit.unknowns(ode))], " ")
expr = "@variables $str"
Y = eval(Meta.parse(expr))
measured_quantities = [
  Y[i] ~ ModelingToolkit.unknowns(ode)[i]
  for i in 1:n_observed
]

res = ode_to_poly(ode, measured_quantities)
poly = res["polynomial_system"]
states = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(species)), gens(parent(poly[1])))
parameters = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(pars)), gens(parent(poly[1])))
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

