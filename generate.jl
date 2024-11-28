using ParameterEstimation, Pkg

Pkg.activate("env")

ParameterEstimation.OUTPUT_DIR[] = "/home/ademin/parameter_estimation_systems/examples/global"
PREFIX = "/home/ademin/ParameterEstimation.jl/examples/all-global"


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

