using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

using Catalyst
using OffsetArrays

n = 3
t = default_t()

# states

str = join(["X"*string(i)*"(t)" for i in 1:n], " ")
expr = "@species $str"
X = eval(Meta.parse(expr))

# params 

str = join(["a"*string(i)*"" for i in 1:n], " ")
expr = "@parameters $str"
A = eval(Meta.parse(expr))

str = join(["b"*string(i)*"" for i in 1:n], " ")
expr = "@parameters $str"
B = eval(Meta.parse(expr))

rxs = vcat(
  [Reaction(A[i], [X[i-1]], [X[i]]) for i in 2:n],
  [Reaction(B[i], [X[i]], [X[i-1]]) for i in 2:n],
  [Reaction(1., nothing, [X[1]])],
  [Reaction(1., [X[n]], nothing)],
)

@named crn = ReactionSystem(rxs, t)
crn = complete(crn)

crn = @reaction_network begin
    (a, 1), X1 <--> X2
    (1, d), X2 <--> X3
    (e, 1), X3 <--> X1
end

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


