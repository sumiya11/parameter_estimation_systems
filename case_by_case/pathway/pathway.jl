using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

using Catalyst
using OffsetArrays

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
  [Reaction(C0[i], [Y0[i]], [S1[i-1], S1[i]]) for i in 1:n],

  [Reaction(1 / (4i), [F[i], S1[i]], [Y1[i]])    for i in 1:n],
  [Reaction(1 / (2i), [Y1[i]], [F[i], S1[i]])    for i in 1:n],
  [Reaction(1 / (3i), [Y1[i]], [F[i], S0[i]])    for i in 1:n],
)

@named crn = ReactionSystem(rxs, t)
crn = complete(crn)
model = convert(ODESystem, crn)

str = join(["y"*string(i)*"(t)" for i in 1:length(ModelingToolkit.unknowns(model))], " ")
expr = "@variables $str"
Y = eval(Meta.parse(expr))
measured_quantities = [
  Y[i] ~ ModelingToolkit.unknowns(model)[i]
  for i in 1:length(ModelingToolkit.unknowns(model))
]

ic = [1.0 for _ in ModelingToolkit.unknowns(model)]
time_interval = [0.0, 1.0]
datasize = 50
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.3 for _ in ModelingToolkit.parameters(model)] # True Parameters

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
	solver = solver)


