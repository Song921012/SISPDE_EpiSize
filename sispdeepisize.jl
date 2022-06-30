# Loading packages
using DifferentialEquations, ModelingToolkit, MethodOfLines, DomainSets, Plots, NLsolve, Sundials, ProgressBars, LoopVectorization

# Define functions
## recovery_rate function
function probgeneration(ratio::Function,γ::Function,initS::Function,initI::Function,dx)
    # Parameters, variables, and derivatives
    @parameters t x
    @parameters dS dI brn ϵ
    @variables S(..) I(..)
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    # 1D PDE and boundary conditions
    eq = [Dt(S(t, x)) ~ dS * Dxx(S(t, x)) - ratio(x, brn, ϵ) * γ(x) * S(t, x) * I(t, x) / (S(t, x) + I(t, x)) + γ(x) * I(t, x),
        Dt(I(t, x)) ~ dI * Dxx(I(t, x)) + ratio(x, brn, ϵ) * γ(x) * S(t, x) * I(t, x) / (S(t, x) + I(t, x)) - γ(x) * I(t, x)]
    bcs = [S(0, x) ~ initS(x),
        I(0, x) ~ initI(x),
        Dx(S(t, 0)) ~ 0.0,
        Dx(S(t, 1)) ~ 0.0,
        Dx(I(t, 0)) ~ 0.0,
        Dx(I(t, 1)) ~ 0.0]
    # Space and time domains
    domains = [t ∈ Interval(0.0, 5.0),
        x ∈ Interval(0.0, 1.0)]
    # PDE system
    @named pdesys = PDESystem(eq, bcs, domains, [t, x], [S(t, x), I(t, x)], [dS => 0.5, dI => 0.1, brn => 3.0, ϵ => 1.0])
    # Method of lines discretization
    # Need a small dx here for accuracy
    order = 2
    discretization = MOLFiniteDifference([x => dx], t)
    # Convert the PDE problem into an ODE problem
    prob = discretize(pdesys, discretization)
    return prob
end

function γ(x)
    y = 1.0
    return y
end
## initial value function
function initS(x)
    y = 0.9 + 0.1 * sin(2 * pi * x)
    return y
end
function initI(x)
    y = 0.1 + 0.1 * cos(2 * pi * x)
    return y
end
function ratio(x, brn, ϵ)
    y = brn + ϵ * sin(2 * pi * x)
    return y
end
prob = probgeneration(ratio, γ, initS, initI, 0.05)
steadystateprob = SteadyStateProblem(prob)
state = solve(steadystateprob, DynamicSS(Tsit5()))
sum(state)

function episize(dS, dI)
    newprob = remake(prob, p=[dS, dI, 3.0, 2.0])
    steadystateprob = SteadyStateProblem(newprob)
    state = solve(steadystateprob, DynamicSS(Tsit5()))
    y = sum(state[20:end]) / 19
    return y
end
episize(1.0, 1.0)

@time episize(exp(-5.0), 2.0)

