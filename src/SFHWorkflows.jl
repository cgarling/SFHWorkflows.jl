module SFHWorkflows

export fit_sfh

include("parse_config.jl") # Code to parse configuration file
include("script.jl")       # Main code to run SFH fitting from provided configuration

end # module