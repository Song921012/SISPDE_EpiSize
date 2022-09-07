##
include("./src/functions.jl")

## 
# Case one
function γ(x)
    y = 1.0
    return y
end
function initS(x)
    y = 0.9 + 0.1 * sin(2 * pi * x)
    return y
end
function initI(x)
    y = 0.1 + 0.1 * cos(2 * pi * x)
    return y
end
dx = 0.05
n = 20

function ratio(x, brn, ϵ)
    y = brn + ϵ * sin(2 * pi * x)
    return y
end
prob = probgeneration!(ratio, γ, initS, initI, dx)
probsinf = sinfprobgeneration!(ratio, γ, initI, dx)
probiinf = iinfprobgeneration!(ratio, γ, initS, dx)
brn = 3.0
ϵ = 2.0
p = [1.0, 1.0, brn, ϵ]
psinf = [1.0, brn, ϵ]
piinf = [1.0, brn, ϵ]

@time episize!(prob, p, n);
@time sinfepisize!(probsinf, psinf, n);
@time iinfepisize!(probiinf, piinf, n);
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
display(plotdI(Ilim, epiresultI, vartype))
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
display(plotdI(Ilim, epiresultI, vartype))
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
display(plotdI(Ilim, epiresultI, vartype))

##
# Test
# Level Set of (dS,dI)
Ilim = Dict("min" => -5.0, "max" => 10.0, "len" => 100)
Slim = Dict("min" => -5.0, "max" => 10.0, "len" => 100)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
#@save "./output/case2/levelsi2.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 100
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
#savefig("./output/case2/levelsi2.png")
##
# ds to infinity
brn = 3.0
ϵ = 2.0
p = [1.0, 1.0, brn, ϵ]
psinf = [1.0, brn, ϵ]
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case1/sinfepidi.png")

##
# ds to infinity case 2
brn = 3.0
ϵ = 3.0
p = [1.0, 1.0, brn, ϵ]
psinf = [1.0, brn, ϵ]
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case1/sinfepidi2.png")

##
# ds to infinity
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case1/iinfepids.png")

##
# Monontone of BRN
Ilim = Dict("min" => 2.0, "max" => 50.0, "len" => 100)
vartype = "R"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label=L"Epidemic size of $c$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"$c$ (basic reproduction number in homogenous environment)")
ylabel!("Epidemic size")
savefig("./output/case1/epibrn.png")

##
# Level Set of (dS,dI)
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
Slim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case1/levelsi1.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 20
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
savefig("./output/case1/levelsi1.png")

##
# test
# Level Set of (dS,dI)
Ilim = Dict("min" => -15.0, "max" => 2.0, "len" => 100)
Slim = Dict("min" => -15.0, "max" => 10.0, "len" => 100)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
#@save "./output/case1/levelsi1.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 600
contour(I_range, S_range, epiresultSI, levels=nlevels, lw=4, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
#savefig("./output/case1/levelsi1.png")


using PackageCompiler
create_sysimage([:DifferentialEquations, :Plots, :Symbolics, :ModelingToolkit, :Integrals, :SciMLSensitivity, :MethodOfLines, :NonlinearSolve], sysimage_path="JuliaSysimage.dll", precompile_execution_file="./src/functions.jl")


##
# Level Set of (dI,tau)
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
Slim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
leveltype = "IT"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case1/levelit1.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 20
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(\tau)")
ylabel!(L"\ln(d_{I})")
title!(L"Level set of $(d_{I},\tau)$")
savefig("./output/case1/levelit1.png")

##
# Level Set of (dS,tau)
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
Slim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
leveltype = "ST"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case1/levelst1.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 20
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(\tau)")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},\tau)$")
savefig("./output/case1/levelst1.png")


## 
# Level Set of (T, E)
p = [0.5, 0.05, brn, ϵ]
Ilim = Dict("min" => -2.0, "max" => 2.0, "len" => 20)
Slim = Dict("min" => -5.0, "max" => 2.0, "len" => 20)
leveltype = "TE"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case1/levelte1.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 100
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\epsilon")
ylabel!(L"\ln(\tau)")
title!(L"Level set of $(\tau,\epsilon)$")
savefig("./output/case1/levelte1.png")

##
p = [0.01, 0.05, brn, ϵ]
Ilim = Dict("min" => -3.0, "max" => 3.0, "len" => 100)
vartype = "E"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n)
τ = p[2] / p[1];
A1_value = A1(τ, brn)
A2_value = A2(τ, brn)
A3_value = A3(τ, brn)
println("Parameter values: A1: $A1_value; A2: $A2_value; A3: $A3_value")
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label="Epidemic size")
title!(L"Epidemic size of $\epsilon$")
xlabel!(L"\epsilon")
ylabel!("Epidemic size")
savefig("./output/case1/enve1.png")

## 
p = [10, 0.05, brn, ϵ]
Ilim = Dict("min" => -3.0, "max" => 3.0, "len" => 100)
vartype = "E"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n)
τ = p[2] / p[1];
A1_value = A1(τ, brn)
A2_value = A2(τ, brn)
A3_value = A3(τ, brn)
println("Parameter values: A1: $A1_value; A2: $A2_value; A3: $A3_value")
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label="Epidemic size")
title!(L"Epidemic size of $\epsilon$")
xlabel!(L"\epsilon")
ylabel!("Epidemic size")
savefig("./output/case1/enve2.png")
##
p = [0.5, 0.05, brn, ϵ]
Ilim = Dict("min" => -3.0, "max" => 3.0, "len" => 100)
vartype = "E"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n)
τ = p[2] / p[1];
A1_value = A1(τ, brn)
A2_value = A2(τ, brn)
A3_value = A3(τ, brn)
println("Parameter values: A1: $A1_value; A2: $A2_value; A3: $A3_value")
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label="Epidemic size")
title!(L"Epidemic size of $\epsilon$")
xlabel!(L"\epsilon")
ylabel!("Epidemic size")
savefig("./output/case1/enve3.png")




## 
## Limit of dS to 0 and monotone on ds
p = [exp(-10), 1.0, brn, ϵ]
Ilim = Dict("min" => 0.0, "max" => 15.0, "len" => 100)
vartype = "S"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label=L"Epidemic size of $d_{S}$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case1/dSmon.png")

p = [exp(-15), 0.05, brn, ϵ]
Ilim = Dict("min" => -3.0, "max" => 3.0, "len" => 100)
vartype = "E"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n)

τ = p[2] / p[1];
A1_value = A1(τ, brn)
A2_value = A2(τ, brn)
A3_value = A3(τ, brn)
println("Parameter values: A1: $A1_value; A2: $A2_value; A3: $A3_value")
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label=L"Epidemic size($d_{S}\rightarrow 0$)")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\epsilon")
ylabel!(L"Epidemic size($d_{S}\rightarrow 0$)")
savefig("./output/case1/dSlim.png")

##
p = [0.1, 1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 100)
vartype = "fixtS"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label="Epidemic size")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case1/dSmonchange.png")