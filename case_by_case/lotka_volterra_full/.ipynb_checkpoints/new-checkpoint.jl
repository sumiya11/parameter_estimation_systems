using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Tsit5()

@parameters k01 k12 k13 k14 k21 k31 k41
@variables t x1(t) x2(t) x3(t) x4(t) y1(t) y2(t) y3(t) y4(t)
D = Differential(t)

ic = [1.0, 2.0, 1.0, -1.0]
time_interval = [0.0, 10.0]
datasize = 20
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.2, 0.3, 0.5, 0.6, -0.2, 1.1, 0.02] # True Parameters

states = [x1, x2, x3, x4]
parameters = [k01, k12, k13, k14, k21, k31, k41]
@named model = ODESystem([
                             D(x1) ~ -k01 * x1 + k12 * x2 + k13 * x3 + k14 * x4 - k21 * x1 -
                                     k31 * x1 - k41 * x1,
                             D(x2) ~ -k12 * x2 + k21 * x1,
                             D(x3) ~ -k13 * x3 + k31 * x1,
                             D(x4) ~ -k14 * x4 + k41 * x1],
                         t, states, parameters)
measured_quantities = [y1 ~ x1, y2 ~ x2 + x4, y3 ~ x3]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
                                              p_true, ic, datasize; solver = solver)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample)

#=
using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters a12 a22
@variables t x1(t) x2(t) y1(t)
D = Differential(t)

ic = [100.0, 100.0]
time_interval = [0.0, 1.0]
datasize = 21
sampling_times = range(time_interval[1], time_interval[2], length = datasize)
p_true = [0.02, 0.03] # True Parameters
measured_quantities = [y1 ~ x1^2+x1]
states = [x1, x2]
parameters = [a12, a22]

@named model = ODESystem(
    [D(x1) ~ 0 * x1 + a12 * x2,
	D(x2) ~ 1 * x1 + a22 * x2],
	t, states, parameters)

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
	p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
	solver = solver)
=#
