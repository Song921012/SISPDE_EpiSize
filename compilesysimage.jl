using PackageCompiler
 create_sysimage([:DifferentialEquations,:Plots,:ModelingToolkit,:MethodOfLines], sysimage_path="JuliaSysimage.dll", precompile_execution_file="./src/functions.jl")