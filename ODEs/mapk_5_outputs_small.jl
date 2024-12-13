# MAPK (5 outputs) (with 6 parameters set to constants)
# https://github.com/SciML/StructuralIdentifiability.jl/blob/7faf054d535e84165e0bdda8c0699dc098183573/benchmarking/IdentifiableFunctions/benchmarks.jl#L64

using ParameterEstimation
using ModelingToolkit, DifferentialEquations
solver = Vern9()

# Set the following parameters to 1.0 : a00, a01, a10, alpha01, alpha10, alpha11.
@parameters a00 a01 a10 alpha01 alpha10 alpha11 b00 b01 b10 beta01 beta10 beta11 c0001 c0010 c0011 c0111 c1011 gamma0100 gamma1000 gamma1100 gamma1101 gamma1110
@variables t KS00(t) KS01(t) KS10(t) FS01(t) FS10(t) FS11(t) K(t) F(t) S00(t) S01(t) S10(t) S11(t) y1(t) y2(t) y3(t) y4(t) y5(t)
D = Differential(t)
states = [KS00, KS01, KS10, FS01, FS10, FS11, K, F, S00, S01, S10, S11]
parameters = [b00, b01, b10, beta01, beta10, beta11, c0001, c0010, c0011, c0111, c1011, gamma0100, gamma1000, gamma1100, gamma1101, gamma1110]

@info "MAPK (5 outputs): $(length(states)) states, $(length(parameters)) parameters"

@named model = ODESystem(
        [
         D(KS00) ~
                -0.1 * K * S00 +
                b00 * KS00 +
                gamma0100 * FS01 +
                gamma1000 * FS10 +
                gamma1100 * FS11,
         D(KS01) ~
                -0.1 * K * S01 + b01 * KS01 + c0001 * KS00 -
                0.1 * F * S01 +
                beta01 * FS01 +
                gamma1101 * FS11,
         D(KS10) ~
                -0.1 * K * S10 + b10 * KS10 + c0010 * KS00 -
                0.1 * F * S10 +
                beta10 * FS10 +
                gamma1110 * FS11,
         D(FS01) ~
                -0.1 * F * S11 +
                beta11 * FS11 +
                c0111 * KS01 +
                c1011 * KS10 +
                c0011 * KS00,
                D(FS10) ~ 0.1 * K * S00 - (b00 + c0001 + c0010 + c0011) * KS00,
                D(FS11) ~ 0.1 * K * S01 - (b01 + c0111) * KS01,
          D(K) ~ 0.1 * K * S10 - (b10 + c1011) * KS10,
          D(F) ~ 0.1 * F * S01 - (beta01 + gamma0100) * FS01,
          D(S00) ~ 0.1 * F * S10 - (beta10 + gamma1000) * FS10,
          D(S01) ~
                0.1 * F * S11 -
                (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11,
          D(S10) ~
                -0.1 * K * S00 + (b00 + c0001 + c0010 + c0011) * KS00 -
                0.1 * K * S01 + (b01 + c0111) * KS01 -
               0.1 * K * S10 + (b10 + c1011) * KS10,
          D(S11) ~
                -0.1 * F * S01 + (beta01 + gamma0100) * FS01 -
                0.1 * F * S10 + (beta10 + gamma1000) * FS10 -
                0.1 * F * S11 +
                (beta11 + gamma1101 + gamma1110 + gamma1100) * FS11,
         ], t, states, parameters)
measured_quantities = [y1 ~ F,
            y2 ~ S00,
            y3 ~ S01,
            y4 ~ S10,
            y5 ~ S11]

ic = [1.0 for i in 1:length(states)]
time_interval = [0.0, 1.0]
datasize = 21
p_true = [0.1 + i*(1 / (2length(parameters))) for i in 1:length(parameters)]

@info "Truth:" ic p_true

data_sample = ParameterEstimation.sample_data(model, measured_quantities, time_interval,
        p_true, ic, datasize; solver = solver, abstol = 1.0e-14, reltol = 1.0e-14)
res = ParameterEstimation.estimate(model, measured_quantities, data_sample;
        solver = solver)

