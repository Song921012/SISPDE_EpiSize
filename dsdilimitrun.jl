##
@time include("./src/functions.jl")


# Test ds to infinity
##  
function γ(x)
    y = exp(sin(2 * pi * x)) + (1 - x)
    return y
end
function ratio(x, η, ν)
    y = 1 + (η + x)^ν
    return y
end
probsinf = sinfprobgeneration!(ratio, γ, initI, dx)
η = 0.3
ν = 0.01
psinf = [1.0, ν, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI, label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$ $\nu=0.01$"))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
#savefig("./output/case2/sinfepidi2.png")

# Test di to infinity
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
display(plot(I_range, epiresultI,lw=3,foreground_color_legend = nothing, label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$ $\nu=1.5$" ))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case2/iinfepidi2.png")

# Test di to infinity
##  
function γ(x)
    y = exp(sin(2 * pi * x)) + (1 - x)
    return y
end
function ratio(x, η, ν)
    y = (η+ν)/γ(x) +1
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
display(plot(I_range, epiresultI, lw=3,foreground_color_legend = nothing,label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$ $\nu=0.9$" ))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case2/iinfepidi3.png")


# Test di to infinity
## 
function γ(x)
    y = 1.0
    return y
end
function ratio(x, brn, ϵ)
    y = brn + ϵ * (exp(sin(2 * pi * x))+1.0-x)
    return y
end
probiinf = iinfprobgeneration!(ratio, γ, initI, dx)
brn = 3.0
ϵ = 2.0
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -0.0, "max" => 5.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI))

#label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$ $\epsilon=2.0$" ))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case1/iinfepidi4.png")

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
ν = 1.0
piinf = [1.0, ν, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI, label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$"))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case2/iinfepidi2.png")


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
plot(I_range, epiresultI, lw=3,foreground_color_legend = nothing,label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case1/sinfepidi.png")
##
# di to infinity
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3,foreground_color_legend = nothing,label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case1/iinfepids.png")

##
# ds to infinity change epsilon
brn = 3.0
ϵ = 3.0
p = [1.0, 1.0, brn, ϵ]
psinf = [1.0, brn, ϵ]
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3,foreground_color_legend = nothing,label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$")
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case1/sinfepidi2.png")


##
# di to infinity
brn = 3.0
ϵ = 3.0
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
display(plot(I_range, epiresultI, label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$"))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case1/iinfepids2.png")



##
# ds to infinity change epsilon
brn = 3.0
ϵ = 3.0
p = [1.0, 1.0, brn, ϵ]
psinf = [1.0, brn, ϵ]
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "I"
@time epiresultI = sinfepisingle!(probsinf, vartype, Ilim, psinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])
plot(I_range, epiresultI, lw=3,foreground_color_legend = nothing,label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$")
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




# beta and gamma are constants
# di goes to infinity
##  
function γ(x)
    y = 1+x
    return y
end
function ratio(x, η, ν)
    y = 1+ (η+ν)/γ(x)
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

display(plot(I_range, epiresultI,lw=3,foreground_color_legend = nothing, label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$" ))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case4/iinfepidi4.png")



# ds to infinity
##
function γ(x)
    y = 1.0
    return y
end
function ratio(x, brn, ϵ)
    y = brn + ϵ * (1.0+x)
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
display(plot(I_range, epiresultI, lw=3,foreground_color_legend = nothing,label=L"Epidemic size of $d_{I}$ as $d_{S} \rightarrow \infty$"))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{I})")
ylabel!("Epidemic size")
savefig("./output/case1/sinfepidi6.png")



# di to infinity test monotonocity
##  
function γ(x)
    y = 1.0
    return y
end
function ratio(x, η, ν)
    y = brn+ ϵ*(sin(x)+exp(sin(x)))
    return y
end
probiinf = iinfprobgeneration!(ratio, γ, initI, dx)
η = brn
ν = ϵ
piinf = [1.0, brn, ϵ]
Ilim = Dict("min" => -10.0, "max" => 10.0, "len" => 50)
vartype = "S"
@time epiresultI = iinfepisingle!(probiinf, vartype, Ilim, piinf, n);
I_range = range(Ilim["min"], Ilim["max"], length=Ilim["len"])

display(plot(I_range, epiresultI,lw=3,foreground_color_legend = nothing, label=L"Epidemic size of $d_{S}$ as $d_{I} \rightarrow \infty$" ))
#title!(L"Epidemic size of $d_{S}$")
xlabel!(L"\ln(d_{S})")
ylabel!("Epidemic size")
savefig("./output/case1/iinfepidi6.png")