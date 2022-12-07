##
@time include("./src/functions.jl")
path = "./output/case2/"
case = 2

## 
# Case one
# 
function γ(x)
    y = exp(sin(2 * pi * x)) + (1 - x)
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
##
η = 0.3
ν = 0.9
p = [exp(10.0),1.0, η, ν]

@time episize!(prob, p, n);
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = episingle!(prob, vartype, Ilim, p, n);
plotdI(Ilim, epiresultI, vartype)
##
# Test
# Level Set of (dS,dI)
η = 0.1
ν = 4.0
p = [exp(10.0),1.0, η, ν]
Ilim = Dict("min" => -10, "max" => 10, "len" => 100)
Slim = Dict("min" => -10, "max" => 10, "len" => 100)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
#@save "./output/case2/levelsi2.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = [0.1865614,0.2060475,0.2255336,0.2308393,0.2328650,0.2429940,0.1475892,0.1378461,0.1315759,0.1295501,0.1275243,0.1254985, 0.1113180,0.0667508]#,0.0647250] # η = 0.3 ν = 0.9 si3
contour(I_range, S_range, epiresultSI, lw=3,levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
savefig("./output/case2/levelsi4.png")
##
# Test
# Level Set of (dS,dI)
η = 0.3
ν = 0.9
p = [exp(10.0),1.0, η, ν]
Ilim = Dict("min" => -10, "max" => 10.0, "len" => 100)
Slim = Dict("min" => -10, "max" => 10, "len" => 100)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
#@save "./output/case2/levelsi2.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = [0.412602,0.414626,0.418039,0.419987,0.428799,0.429724,0.430021,0.430247,0.439462,0.449200,0.447252,0.452285,0.451070] # η = 0.3 ν = 0.9 si3
contour(I_range, S_range, epiresultSI, lw=3,levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
savefig("./output/case2/levelsi3.png")
##
# Level Set of (dS,dI)
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
Slim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case2/levelsi2.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 20
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
savefig("./output/case2/levelsi2.png")


##
# test
# Level Set of (dS,dI)
Ilim = Dict("min" => -6.0, "max" => 3.0, "len" => 100)
Slim = Dict("min" => -15.0, "max" => -5.0, "len" => 100)
leveltype = "SI"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
#@save "./output/case1/levelsi1.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 600
contour(I_range, S_range, epiresultSI, levels=nlevels,lw=4, contour_labels=true)
xlabel!(L"\ln(d_{I})")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},d_{I})$")
#savefig("./output/case1/levelsi1.png")

##
# Level Set of (dI,tau)
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
Slim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
leveltype = "IT"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case2/levelit2.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 20
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(\tau)")
ylabel!(L"\ln(d_{I})")
title!(L"Level set of $(d_{I},\tau)$")
savefig("./output/case2/levelit2.png")

##
# Level Set of (dS,tau)
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
Slim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
leveltype = "ST"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case2/levelst2.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 20
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\ln(\tau)")
ylabel!(L"\ln(d_{S})")
title!(L"Level set of $(d_{S},\tau)$")
savefig("./output/case2/levelst2.png")


## 
# Level Set of (T, E)
p = [0.5, 0.05, brn, ϵ]
Ilim = Dict("min" => -2.0, "max" => 2.0, "len" => 20)
Slim = Dict("min" => -5.0, "max" => 2.0, "len" => 20)
leveltype = "TE"
@time epiresultSI = levelset(prob, leveltype, Slim, Ilim, p, n);
@save "./output/case2/levelte2.bson" epiresultSI
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
S_range = range(Slim["min"], Slim["max"], length=Slim["len"])
nlevels = 100
contour(I_range, S_range, epiresultSI, levels=nlevels, contour_labels=true)
xlabel!(L"\epsilon")
ylabel!(L"\ln(\tau)")
title!(L"Level set of $(\tau,\epsilon)$")
savefig("./output/case2/levelte2.png")