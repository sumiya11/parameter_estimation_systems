using Pkg

Pkg.activate("env")

using ParameterEstimation

ParameterEstimation.OUTPUT_DIR[] = "examples/new"
PREFIX = joinpath((@__DIR__), "ODEs")


# Only globally identifiable examples
for (dir, _, files) in walkdir(PREFIX)
	for file in files
		!endswith(file, ".jl") && continue
	 	@info "Processing $file"
		ParameterEstimation.MODEL[] = chop(file, tail=3)
		try
			include(PREFIX * "/" * file)
		catch e
			@warn "ParameterEstimation.jl errored from $file, skipping" e
		end
	end
end

