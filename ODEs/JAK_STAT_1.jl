# JAK_STAT_1
# https://github.com/SciML/StructuralIdentifiability.jl/blob/d34f4dd4c858f3ec7f816bfa749cfc4546793f01/benchmarking/IdentifiableFunctions/benchmarks.jl#L571C17-L592C10

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters u t1 t10 t11 t12 t13 t14 t15 t16 t17 t18 t19 t2 t20 t21 t22 t3 t4 t5 t6 t7 t8 t9
@variables t x1(t) x2(t) x3(t) x4(t) x5(t) x6(t) x7(t) x8(t) x9(t) x10(t) y1(t) y2(t) y3(t) y4(t) y5(t) y6(t) y7(t) y8(t) y9(t)
D = Differential(t)
states = [x1, x2, x3, x4, x5, x6, x7, x8, x9, x10]
parameters = [u, t1, t10, t11, t12, t13, t14, t15, t16, t17, t18, t19, t2, t20, t21, t22, t3, t4, t5, t6, t7, t8, t9]

@info "JAK_STAT_1: $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
                D(x1) ~ -t1 * x1 * 2 * u - t5 * x1 + t6 * x2,
                D(x2) ~ t5 * x1 - t6 * x2,
                D(x3) ~ t1 * 2 * u * x1 - t2 * x3 * (-x6 + 3),
                D(x4) ~ t2 * x3 * (-x6 + 3) - t3 * x4,
                D(x5) ~ t3 * x4 - t4 * x5,
                D(x6) ~
                    -t7 * x3 * x6 / (1 + t13 * x1) -
                    t7 * x4 * x6 / (1 + t13 * x10) + t8 * (-x6 + 3) * 92,
                D(x7) ~ -t9 * x7 * (-x6 + 3) + t10 * (-x7 + 165) * 92,
                D(x8) ~ t11 * (-x7 + 165),
                D(x9) ~ -t12 * 2 * u * x9,
		D(x10) ~ x8 * t14 / (t15 + x8) - t16 * x10 
	], t, states, parameters)
measured_quantities = [
    y1 ~ x1 + x3 + x4,
    y2 ~ t18 * (x3 + x4 + x5 + (1 / 3 - x9)),
    y3 ~ t19 * (x4 + x5),
    y4 ~ t20 * (-x6 + 3),
    y5 ~ t21 * x8,
    y6 ~ t22 * x8 * t17 / t11,
    y7 ~ x10,
    y8 ~ -x7 + 165,
]

ic = [1.0 for i in 1:length(states)]
time_interval = [0.0, 1.0]
datasize = 21
p_true = [0.1 + i*(1 / (2length(parameters))) for i in 1:length(parameters)]
data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)

@info "" data_sample
println(parameters)
println(p_true)

res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
        solver = solver)

