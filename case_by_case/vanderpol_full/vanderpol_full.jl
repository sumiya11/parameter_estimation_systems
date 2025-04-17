using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters a b
@variables t x1(t) x2(t) y1(t) y2(t)
D = Differential(t)
# TODO
states = [x1, x2]
parameters = [a, b]
@named model = ODESystem([
                             D(x1) ~ a * x2,
                             D(x2) ~ -(x1) - b * (x1^2 - 1) * (x2),
                         ], t, states, parameters)
measured_quantities = [
        y1 ~ x1,
        y2 ~ x2,
]

ic = [0.582, 0.536]
p_true = [0.539, 0.672]
time_interval = [-0.5, 0.5]
datasize = 21

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
	solver = solver)

