using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

using Catalyst

t = default_t()

crn = @reaction_network begin
    (a, b), X1 <--> X3 + X4
    (c, d), X1 + X2 <--> X4
end

model = convert(ODESystem, crn)

str = join(["y"*string(i)*"(t)" for i in 1:length(ModelingToolkit.unknowns(model))], " ")
expr = "@variables $str"
Y = eval(Meta.parse(expr))
measured_quantities = [
  Y[i] ~ ModelingToolkit.unknowns(model)[i]
  for i in 1:length(ModelingToolkit.unknowns(model))
]
measured_quantities = [
  Y[1] ~ ModelingToolkit.unknowns(model)[1]
]

ic = [1/i for (i,) in enumerate(ModelingToolkit.unknowns(model))]
time_interval = [0.0, 1.0]
datasize = 50
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.3 for _ in ModelingToolkit.parameters(model)] # True Parameters

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
	solver = solver)


