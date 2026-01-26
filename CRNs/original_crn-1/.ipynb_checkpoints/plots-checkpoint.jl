using Catalyst, GraphMakie, NetworkLayout
using CairoMakie

crn = @reaction_network begin
    (k1, k2), S0 + K <--> S0K
    k3, S0K --> Sn + K
    (k4, k5), Sn + F <--> SnF
    k6, SnF --> S0 + F
end

f, ax, p = plot_complexes(crn; show_rate_labels = true, elabels_fontsize = 8)
hidedecorations!(ax); hidespines!(ax)
colsize!(f.layout, 1, Aspect(1, 0.6))
save(joinpath(@__DIR__, "graph.png"), f)
