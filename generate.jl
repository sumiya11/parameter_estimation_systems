using Pkg
Pkg.activate("new")

using Revise
using ParameterEstimation

ParameterEstimation.OUTPUT_DIR[] = "examples/new"
PREFIX = joinpath((@__DIR__), "ODEs")


# Only globally identifiable examples
for (dir, _, files) in walkdir(PREFIX)
	for file in files
		!endswith(file, ".jl") && continue
		if file in ["mapk_5_outputs_lilliput.jl", "JAK_STAT_1.jl", "PK1.jl", "QY.jl", "SEIR_36_ref.jl", "akt_pathway.jl", "akt_pathway_small.jl", "mapk_5_outputs.jl", "mapk_5_outputs_small.jl", "mapk_5_outputs_tiny.jl", "treatment.jl", "nfkb.jl", "akt_pathway_small.jl", ] || occursin("nfkb", file)
			@info "Skipping $file"
			continue
		end
		@error "Processing $file"
		ParameterEstimation.MODEL[] = chop(file, tail=3)
		try
			include(PREFIX * "/" * file)
		catch e
			@warn "ParameterEstimation.jl errored from $file, skipping" e
		end
	end
end

