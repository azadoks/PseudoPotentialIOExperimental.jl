Pkg.activate("play")
begin
    using PseudoPotentialIOExperimental
    import PseudoPotentialIOExperimental: ArbitraryMesh, UniformMesh, LogMeshWithZero, LogMeshWithoutZero, deriv
    import PseudoPotentialIOExperimental: NumericLocalPotential, NumericProjector, NumericDensity
    import PseudoPotentialIOExperimental: hankel_transform
    import PseudoPotentialIOExperimental: Simpson, Trapezoid, QESimpson, integration_weights
    import PseudoPotentialIOExperimental: CoulombCorrection, ErfCoulombCorrection
    import PseudoPotentialIOExperimental: RealSpace, FourierSpace
    using CairoMakie
    using CubicSplines
    using LinearAlgebra
    using QuadGK
    using SpecialFunctions
end

function linear_grid_exp_sin(n; xn=Float64(π))
    x = UniformMesh(0.0, xn, n)
    f(x) = exp(-(x^2)) * sin(x)
    y = f.(x)

    I_true = real(-(√π * (-2erfi(1/2) + erfi(1/2 - im*xn) + erfi(1/2 + im*xn))) / (4 * exp(1/4)))

    return (; x, f, y, I_true)
end

function log_grid_exp_sin(n; xn=Float64(π))
    x = LogMeshWithoutZero(0.0, xn, n)

    # i = collect(1:n)
    # b = xn / n
    # a = log(xn / b) / (n - 1)
    # x = logarithmic_mesh1.(1:n, a, b)
    # dx = @. a * b * exp.(a * (i - 1))
    # Δx = diff(x)
    f(x) = exp(-(x^2)) * sin(x)
    y = f.(x)

    I_true = real(-(√π * (-2erfi(1/2) + erfi(1/2 - im*xn) + erfi(1/2 + im*xn))) / (4 * exp(1/4)))

    return (; x, f, y, I_true)
end
