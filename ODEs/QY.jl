# QY
# https://github.com/SciML/StructuralIdentifiability.jl/blob/d34f4dd4c858f3ec7f816bfa749cfc4546793f01/benchmarking/IdentifiableFunctions/benchmarks.jl#L347C17-L425C10

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

@parameters Ks M Mar alpa beta beta_SA beta_SI phi siga1 siga2
@variables t P0(t) P1(t) P2(t) P3(t) P4(t) P5(t) y1(t)
D = Differential(t)
states = [P0, P1, P2, P3, P4, P5]
parameters = [Ks, M, Mar, alpa, beta, beta_SA, beta_SI, phi, siga1, siga2]

@info "QY: $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
            D(P0) ~ P1,
            D(P1) ~ P2,
            D(P2) ~ P3,
            D(P3) ~ P4,
            D(P4) ~
                -(
                    Ks * M * siga1 * siga2 * P1 +
                    (
                        Ks * M * siga1 +
                        Ks * M * siga2 +
                        Ks * siga1 * siga2 +
                        siga1 * siga2 * M
                    ) * P2 +
                    (
                        Ks * M +
                        Ks * siga1 +
                        Ks * siga2 +
                        M * siga1 +
                        M * siga2 +
                        siga1 * siga2
                    ) * P3 +
                    (Ks + M + siga1 + siga2) * P4
                ) -
                (
                    Mar * P5 +
                    beta +
                    beta_SA / (siga2 * M) * (
                        P3 +
                        P2 * (Ks + M + Mar) +
                        P1 * (Ks * M + Ks * Mar + M * Mar) +
                        P0 * Ks * M * Mar
                    ) +
                    beta_SI / M * (P2 + P1 * (Ks + Mar) + P0 * Ks * Mar) +
                    beta_SA * phi / ((1 - phi) * siga2 * M) * (
                        P3 +
                        P2 * (Ks + M + siga2) +
                        P1 * (Ks * M + Ks * siga2 + M * siga2) +
                        P0 * Ks * M * siga2
                    )
                ) * (
                    alpa +
                    Ks * M * siga1 * siga2 * P0 +
                    (
                        Ks * M * siga1 +
                        Ks * M * siga2 +
                        Ks * siga1 * siga2 +
                        siga1 * siga2 * M
                    ) * P1 +
                    (
                        Ks * M +
                        Ks * siga1 +
                        Ks * siga2 +
                        M * siga1 +
                        M * siga2 +
                        siga1 * siga2
                    ) * P2 +
                    (Ks + M + siga1 + siga2) * P3 +
                    P4
                ),
            D(P5) ~
                -Mar * P5 - (
                    beta +
                    beta_SA / (siga2 * M) * (
                        P3 +
                        P2 * (Ks + M + Mar) +
                        P1 * (Ks * M + Ks * Mar + M * Mar) +
                        P0 * Ks * M * Mar
                    ) +
                    beta_SI / M * (P2 + P1 * (Ks + Mar) + P0 * Ks * Mar) +
                    beta_SA * phi / ((1 - phi) * siga2 * M) * (
                        P3 +
                        P2 * (Ks + M + siga2) +
                        P1 * (Ks * M + Ks * siga2 + M * siga2) +
                        P0 * Ks * M * siga2
                    )
                ),
         ], t, states, parameters)
measured_quantities = [y1 ~ P0]

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


