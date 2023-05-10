using PackageCompiler
create_sysimage([:DifferentialEquations,:Plots,:ModelingToolkit,:MethodOfLines,:Integrals], sysimage_path="JuliaSysimage.so", precompile_execution_file="./src/functions.jl")q