using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters g a b
@variables t V(t) R(t) y1(t) y2(t)
D = Differential(t)
states = [V, R]
parameters = [g, a, b]

ic = [1.0, -1.0]
time_interval = [0.0, 1.0]
datasize = 50
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [2, 2 / 10, 2 / 10] # True Parameters
measured_quantities = [y1 ~ V]

@named model = ODESystem([
                             D(V) ~ g * (V - V^3 / 3 + R),
                             D(R) ~ 1 / g * (V - a + b * R),
                         ], t, states, parameters)

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
	solver = solver)

