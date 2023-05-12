##
@time include("./src/functions.jl")



# di to infnty
##  
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


Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, legend=:bottomright, label=L"disease prevalence of $d_{S}$ as $d_{I} \rightarrow \infty$")
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("disease prevalence")
savefig("./output/iinfepids.png")
##  
function γ(x)
    y = exp(sin(2 * pi * x)) + (1 - x)
    return y
end
function ratio(x, η, ν)
    y = (η + ν) / γ(x) + 1
    return y
end
probiinf = iinfprobgeneration!(ratio, γ, initI, dx)
η = 0.3
ν = 1.5
piinf = [1.0, ν, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label=L"disease prevalence of $d_{S}$ as $d_{I} \rightarrow \infty$ $\nu=0.9$"))
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("disease prevalence")
savefig("./output/iinfepidi3.png")

##  
function γ(x)
    y = exp(sin(2 * pi * x)) + (1 - x)
    return y
end
function ratio(x, η, ν)
    y = 1 + (η + x)^ν
    return y
end
probiinf = iinfprobgeneration!(ratio, γ, initI, dx)
η = 0.3
ν = 1.5
piinf = [1.0, ν, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label=L"disease prevalence of $d_{S}$ as $d_{I} \rightarrow \infty$ $\nu=1.5$"))
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("disease prevalence")
savefig("./output/iinfepidi2.png")

##
function γ(x)
    y = 1 + x
    return y
end
function ratio(x, η, ν)
    y = 1 + (η + ν) / γ(x)
    return y
end
probiinf = iinfprobgeneration!(ratio, γ, initI, dx)
η = 0.3
ν = 1.5
piinf = [1.0, ν, η]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])

display(plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label=L"disease prevalence of $d_{S}$ as $d_{I} \rightarrow \infty$"))
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("disease prevalence")
savefig("./output/iinfepidi4.png")


## dsmon.png
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
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, legend=:bottomright, label=L"disease prevalence of $d_{S}$ as $d_{I} \rightarrow \infty$")
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("disease prevalence")
savefig("./output/iinfepids.png")
p = [exp(-10), 1.0, brn, ϵ]
Ilim = Dict("min" => 0.0, "max" => 15.0, "len" => 100)
vartype = "S"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label=L"disease prevalence of $d_{S}$")
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("disease prevalence")
savefig("./output/dSmon.png")



# ds to infinity
##
function γ(x)
    y = 1.0 + 0.1 * x + exp(1.0 * cos(2 * pi * x))
    return y
end
function ratio(x, brn, ϵ)
    y = 1 + (brn + ϵ) / γ(x)
    return y
end
probsinf = sinfprobgeneration!(ratio, γ, initS, dx)
brn = 1.0
ϵ = 0.5
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label=L"disease prevalence of $d_{I}$ as $d_{S} \rightarrow \infty$"))
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("disease prevalence")
savefig("./output/sinfepidi9.png")


##
function γ(x)
    y = 1.0 + x
    return y
end
function ratio(x, brn, ϵ)
    y = brn + ϵ * sin(2 * pi * x)
    return y
end
probsinf = sinfprobgeneration!(ratio, γ, initS, dx)
brn = 3.0
ϵ = 2.5
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label=L"disease prevalence of $d_{I}$ as $d_{S} \rightarrow \infty$"))
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("disease prevalence")
savefig("./output/sinfepidi7.png")

##
function γ(x)
    y = 1.0
    return y
end
function ratio(x, brn, ϵ)
    y = brn + ϵ * (1.0 + x)
    return y
end
probsinf = sinfprobgeneration!(ratio, γ, initS, dx)
brn = 3.0
ϵ = 2.0
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label=L"disease prevalence of $d_{I}$ as $d_{S} \rightarrow \infty$"))
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("disease prevalence")
savefig("./output/sinfepidi6.png")




# ds monchange

##
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
    y = brn + ϵ * (sin(2 * pi * x))
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
p = [0.1, 1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 100)
vartype = "fixtS"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label="disease prevalence")
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("disease prevalence")
savefig("./output/dSmonchange0.png")

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
p = [0.1, 1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 100)
vartype = "fixtS"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label="disease prevalence")
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("disease prevalence")
savefig("./output/dSmonchange3.png")


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
    y = brn + ϵ * (sin(2 * pi * x)+sin(4 * pi * x)+sin(6 * pi * x)+1.0-x^2)
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
p = [0.1, 1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 100)
vartype = "fixtS"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3, foreground_color_legend=nothing, label="disease prevalence")
#title!(L"disease prevalence of $d_{S}$")
xlabel!(L"\ln(d_{I})")
display(ylabel!("disease prevalence"))
savefig("./output/dSmonchange2.png")
# level set

##

function γ(x)
    y = exp(sin(2 * pi * x)) + (1 - x)
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

function ratio(x, η, ν)
    y = 1 + (η + x)^ν
    return y
end

prob = probgeneration!(ratio, γ, initS, initI, dx)
η = 0.1
ν = 4.0
p = [exp(10.0), exp(10.0), η, ν]

@time episize!(prob, p, n);
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
plotdI(Ilim, epiresultI, vartype)

vartype = "I"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
plotdI(Ilim, epiresultI, vartype)
η = 0.3
ν = 0.9
p = [exp(10.0), 1.0, η, ν]

@time episize!(prob, p, n);
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
plotdI(Ilim, epiresultI, vartype)
η = 0.1
ν = 4.0
p = [exp(10.0), 1.0, η, ν]
Ilim = Dict("min" => -10, "max" => 10, "len" => 100)
Slim = Dict("min" => -10, "max" => 10, "len" => 100)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
#@save "./output/case2/levelsi2.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = [0.1865614, 0.2060475, 0.2255336, 0.2308393, 0.2328650, 0.2429940, 0.1475892, 0.1378461, 0.1315759, 0.1295501, 0.1275243, 0.1254985, 0.1113180, 0.0667508]#,0.0647250] # η = 0.3 ν = 0.9 si3
contour(I_range, S_range, epiresultSI, lw=3, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!("Level set of disease prevalence")
savefig("./output/reviselevelsi4.png")
η = 0.3
ν = 0.9
p = [exp(10.0), 1.0, η, ν]
Ilim = Dict("min" => -10, "max" => 10.0, "len" => 100)
Slim = Dict("min" => -10, "max" => 10, "len" => 100)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = [0.412602, 0.414626, 0.418039, 0.419987, 0.428799, 0.429724, 0.430021, 0.430247, 0.439462, 0.449200, 0.447252, 0.452285, 0.451070] # η = 0.3 ν = 0.9 si3
contour(I_range, S_range, epiresultSI, lw=3, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!("Level set of disease prevalence")
savefig("./output/reviselevelsi3.png")