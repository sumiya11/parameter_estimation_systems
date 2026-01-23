using Catalyst, GraphMakie, NetworkLayout
using CairoMakie

include("original-crn-2.jl")

f, ax, p = plot_complexes(crn; )
hidedecorations!(ax); hidespines!(ax)
colsize!(f.layout, 1, Aspect(1, 1.0))
save(joinpath(@__DIR__, "graph.png"), f)
