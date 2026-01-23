using ModelingToolkit, Catalyst
using SIAN, Nemo, StructuralIdentifiability, Groebner
using LinearAlgebra, OffsetArrays

include("../make_poly_from_ode.jl")
include("../plug_numbers_in_poly.jl")
include("../quotient_basis.jl")

t = default_t()

function run()
for n in 2:6
# n = 5
@info "n = $n"

# States

str = join(["X"*string(i)*"(t)" for i in 1:n], " ")
expr = "@species $str"
X = eval(Meta.parse(expr))

# Parameters 

str = join(["a"*string(i)*"" for i in 1:n], " ")
expr = "@parameters $str"
A = OffsetArray(eval(Meta.parse(expr)), 1:n)

str = join(["b"*string(i)*"" for i in 1:n], " ")
expr = "@parameters $str"
B = OffsetArray(eval(Meta.parse(expr)), 1:n)

rxs = vcat(
  [Reaction(A[i], [X[i-1]], [X[i]]) for i in 2:n],
  [Reaction(B[i], [X[i]], [X[i-1]]) for i in 2:n],
  [Reaction(A[1], nothing, [X[1]])],
  [Reaction(B[1], [X[n]], nothing)],
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
states = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), collect(X)), res["vars"])
parameters = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(collect(A),collect(B))), res["vars"])
states = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in states]
parameters = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in parameters]
poly_num, subs = plug_numbers_in_poly(poly, states, parameters)

@info "" length(poly_num)
@info "Linear equations" length(filter(x -> x in gens(parent(poly_num[1])), map(leading_monomial, poly_num)))

gb = groebner(poly_num);
@assert length(quotient_basis(gb)) == 1
end
end

run()
