using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters k1 k2 k3
@variables t r(t) w(t) y1(t)
D = Differential(t)

ic = [0.3, 0.4]
time_interval = [0.0, 1.0]
datasize = 21
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.02, 0.03, 0.05] # True Parameters
measured_quantities = [y1 ~ r]
states = [r, w]
parameters = [k1, k2, k3]

@named model = ODESystem([D(r) ~ k1 * r - k2 * r * w,
		D(w) ~ k2 * r * w - k3 * w],
	t, states, parameters)

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
	solver = solver)

