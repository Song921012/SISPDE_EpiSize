##
include("./src/functions.jl")

## 
# Case three
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
    y = brn + ϵ * (sin(2 * pi * x) + sin(4 * pi * x))
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
Ilim = Dict("min" => -5.0, "max" => 5.0, "len" => 50)
vartype = "I"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
display(plotdI(Ilim, epiresultI, vartype))
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
display(plotdI(Ilim, epiresultI, vartype))
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
display(plotdI(Ilim, epiresultI, vartype))

##
# ds to infinity
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3,foreground_color_legend = nothing,label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case3/sinfepidi.png")

##
# Level Set of (dS,dI)
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
Slim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case3/levelsi3.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 20
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
savefig("./output/case3/levelsi3.png")

##
p = [0.1, 1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 100)
vartype = "fixtS"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI,lw=3,foreground_color_legend = nothing, label="Epidemic size")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case3/dSmonchange3.png")


## 
# Level Set of (T, E)
p = [0.5, 0.05, brn, ϵ]
Ilim = Dict("min" => -2.0, "max" => 2.0, "len" => 20)
Slim = Dict("min" => -5.0, "max" => 2.0, "len" => 20)
leveltype = "TE"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case3/levelte3.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 100
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\epsilon")
ylabel!(L"\ln(\tau)")
title!(L"Level set of $(\tau,\epsilon)$")
savefig("./output/case3/levelte3.png")














## 
# Case Four
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
    y = brn + ϵ * (sin(2 * pi * x) + sin(4 * pi * x) + sin(6 * pi * x))
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
Ilim = Dict("min" => -5.0, "max" => 5.0, "len" => 50)
vartype = "I"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
display(plotdI(Ilim, epiresultI, vartype))
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
display(plotdI(Ilim, epiresultI, vartype))
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
display(plotdI(Ilim, epiresultI, vartype))

##
# ds to infinity
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case4/sinfepidi.png")


##
# Level Set of (dS,dI)
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
Slim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case4/levelsi4.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 20
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
savefig("./output/case4/levelsi4.png")

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
savefig("./output/case4/dSmonchange4.png")


# Case 5
##

# ds to infinity
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
# 20221207 Eq2.25
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
    y = brn + ϵ * (10+x)
    return y
end
prob = probgeneration!(ratio, γ, initS, initI, dx)
brn = 3.0
ϵ = 0.1
p = [10.0, 1.0, brn, ϵ]
psinf = [1.0, brn, ϵ]
piinf = [1.0, brn, ϵ]

@time episize!(prob, p, n);
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
display(plotdI(Ilim, epiresultI, vartype))
