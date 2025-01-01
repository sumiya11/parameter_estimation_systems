# TumorHu2019
# https://github.com/SciML/StructuralIdentifiability.jl/blob/85c2f525d4e1eaaebae761319e455c538cf0c843/benchmarking/IdentifiableFunctions/benchmarks.jl#L699

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters d1 k2 k3 k4 k5 m1 m3 m4 mu2 mu3 mu4 mu5 r1 r2 r3 r4
@variables t x(t) y(t) z(t) w(t) v(t) y1(t)
D = Differential(t)
states = [x, y, z, w, v]
parameters = [d1, k2, k3, k4, k5, m1, m3, m4, mu2, mu3, mu4, mu5, r1, r2, r3, r4]

@info "$(@__FILE__): $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
            D(x) ~
                (r1 + 0.1 * y) * x * (1 - 0.1 * x) - d1 * x * z / (m1 + w), #pancreatic cancer cell population
            D(y) ~
                (r2 + 0.1 * w / (k2 + w)) * y * (1 - 0.1 * y) - mu2 * y, #pancreatic stellate cell population
            D(z) ~ 0.1 * z * v / ((k3 + v) * (m3 + w)) - mu3 * z + r3, #effector cells, including CD8+T cells and NK cells
            D(w) ~
                0.1 * x * z / (k4 + x) - mu4 * w +
                r4 * x * y / (m4 + v), #concentration of tumor promoting cytokines, including TGF-beta and IL-6
            D(v) ~ 0.1 * x * z / (k5 + x) - mu5 * v, #concentration of tumor suppressing cytokines, including INF-gamma and IL-2
         ], t, states, parameters)
measured_quantities = [
    y1 ~ z
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

