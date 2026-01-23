using Catalyst, GraphMakie, NetworkLayout
using CairoMakie

include("pathway-5.jl")

f, ax, p = plot_complexes(crn)
hidedecorations!(ax); hidespines!(ax)
colsize!(f.layout, 1, Aspect(1, 0.6))
save(joinpath(@__DIR__, "graph.png"), f)
