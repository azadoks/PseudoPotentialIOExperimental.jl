abstract type NumericPsPQuantity{T<:Real,S<:EvaluationSpace} <: AbstractPsPQuantity end

minimum_radius(q::NumericPsPQuantity) = minimum(q.r; init=Inf)
maximum_radius(q::NumericPsPQuantity) = maximum(q.r; init=0)
extrema_radii(q::NumericPsPQuantity) = (minimum_radius(q), maximum_radius(q))

hankel_transform(::Nothing; kwargs...)::Nothing = nothing
function hankel_transform(
    quantity::NumericPsPQuantity{T,S};
    qs::AbstractVector{TT}=range(start=T(0), stop=T(30), length=3001),
    quadrature_method::QuadratureMethod=Simpson()
) where {T<:Real,S<:RealSpace,TT<:Real}
    n = length(quantity.r)
    work_weights=Vector{T}(undef, n)
    work_integrand=Vector{TT}(undef, n)
    work_f=Vector{T}(undef, n)
    return hankel_transform(quantity, qs, quadrature_method, work_weights, work_integrand, work_f)
end

cubic_spline(::Real, ::Nothing)::Nothing = nothing
function cubic_spline(
    maximum_spacing::T,
    q::NumericPsPQuantity{T,S}
)::NumericPsPQuantity{T,S} where {T<:Real, S<:EvaluationSpace}
    r = UniformMesh(extrema_radii(q)..., maximum_spacing)
    return cubic_spline(r, q)
end

function (quantity::NumericPsPQuantity{T,S})(
    x::TT
)::TT where {T<:Real,S<:EvaluationSpace,TT<:Real}
    return quantity.itp(x)
end

function (quantity::NumericPsPQuantity{T,S})(
    y::AbstractVector{TT}, x::AbstractVector{TT}
) where {T<:Real,S<:EvaluationSpace,TT<:Real}
    return quantity.itp(y, x)
end

struct NumericLocalPotential{
    T<:Real,
    S<:EvaluationSpace
} <: NumericPsPQuantity{T,S}
    r::AbstractVector{T}
    f::AbstractVector{T}
    itp::CubicSpline{T}
end

function NumericLocalPotential{T,S}(r::AbstractVector{T}, f::AbstractVector{T}) where {T<:Real,S<:EvaluationSpace}
    itp = CubicSpline(r, f)
    return NumericLocalPotential{T,S}(r, f, itp)
end

function NumericLocalPotential(r::AbstractVector{T}, f::AbstractVector{T}) where {T<:Real}
    return NumericLocalPotential{T,RealSpace}(r, f)
end

function hankel_transform(
    quantity::NumericLocalPotential{T,S},
    correction::LocalPotentialCorrection{T,S};
    qs::AbstractVector{TT}=range(start=T(0), stop=T(30), length=3001),
    quadrature_method::QuadratureMethod=Simpson()
) where {T<:Real,S<:RealSpace,TT<:Real}
    n = length(quantity.r)
    work_weights=Vector{T}(undef, n)
    work_integrand=Vector{TT}(undef, n)
    work_f=Vector{T}(undef, n)
    return hankel_transform(quantity, correction, qs, quadrature_method, work_weights, work_integrand, work_f)
end

function hankel_transform(
    quantity::NumericLocalPotential{T,S},
    correction_r::LocalPotentialCorrection{T,S},
    qs::AbstractVector{TT},
    quadrature_method::QuadratureMethod,
    work_weights::AbstractVector{T},
    work_integrand::AbstractVector{TT},
    work_f::AbstractVector{T}
) where {T<:Real,S<:RealSpace,TT<:Real}
    mesh = quantity.r
    r = collect(mesh)
    f_r = quantity.f
    work_f .= r .* f_r .- correction_r.(r)

    work_weights = @view work_weights[eachindex(r)]
    work_integrand = @view work_integrand[eachindex(r)]

    integration_weights!(work_weights, mesh, quadrature_method)
    correction_q = hankel_transform(correction_r)

    f_q = map(qs) do q
        work_integrand .= work_f .* sin.(q .* r)
        integral = dot(work_weights, work_integrand)
        return 4TT(π) * (integral / q + correction_q(q))
    end

    return NumericLocalPotential{TT,FourierSpace}(qs, f_q)
end

function cubic_spline(
    r::AbstractVector{T},
    Vloc_r::NumericLocalPotential{T,S},
) where {T<:Real,S<:EvaluationSpace}
    f = Vloc_r.(r)
    return NumericLocalPotential{T,S}(r, f, Vloc_r.itp)
end

struct NumericProjector{
    T<:Real,
    S<:EvaluationSpace
} <: NumericPsPQuantity{T,S}
    n::Int
    l::Int
    j::T
    r::AbstractVector{T}
    f::AbstractVector{T}
    itp::CubicSpline{T}
end

function NumericProjector{T,S}(n::Int, l::Int, j::T, r::AbstractVector{T}, f::AbstractVector{T}) where {T<:Real,S<:EvaluationSpace}
    itp = CubicSpline(r, f)
    return NumericProjector{T,S}(n, l, j, r, f, itp)
end

function NumericProjector(n::Int, l::Int, j::T, r::AbstractVector{T}, f::AbstractVector{T}) where {T<:Real}
    return NumericProjector{T,RealSpace}(n, l, j, r, f)
end

function hankel_transform(
    quantity::NumericProjector{T,S},
    qs::AbstractVector{TT},
    quadrature_method::QuadratureMethod,
    work_weights::AbstractVector{T},
    work_integrand::AbstractVector{TT},
    ::AbstractVector{T}
) where {T<:Real,S<:RealSpace,TT<:Real}
    n = quantity.n
    l = quantity.l
    j = TT(quantity.j)
    mesh = quantity.r
    r = collect(mesh)
    r²f_r = quantity.f

    work_weights = @view work_weights[eachindex(r)]
    work_integrand = @view work_integrand[eachindex(r)]

    integration_weights!(work_weights, mesh, quadrature_method)
    jₗ = fast_sphericalbesselj(l)

    f_q = map(qs) do q
        work_integrand .= r²f_r .* jₗ.(q .* r)
        integral = dot(work_weights, work_integrand)
        return 4TT(π) * integral
    end

    return NumericProjector{TT,FourierSpace}(n, l, j, qs, f_q)
end

function cubic_spline(
    r::AbstractVector{T},
    quantity::NumericProjector{T,S},
) where {T<:Real,S<:EvaluationSpace}
    f = quantity.(r)
    return NumericProjector{T,S}(quantity.n, quantity.l, quantity.j, r, f, quantity.itp)
end


struct NumericDensity{
    T<:Real,
    S<:EvaluationSpace
} <: NumericPsPQuantity{T,S}
    r::AbstractVector{T}
    f::AbstractVector{T}
    itp::CubicSpline{T}
end

function NumericDensity{T,S}(r::AbstractVector{T}, f::AbstractVector{T}) where {T<:Real,S<:EvaluationSpace}
    itp = CubicSpline(r, f)
    return NumericDensity{T,S}(r, f, itp)
end

function NumericDensity(r::AbstractVector{T}, f::AbstractVector{T}) where {T<:Real}
    return NumericDensity{T,RealSpace}(r, f)
end

function hankel_transform(
    quantity::NumericDensity{T,S},
    qs::AbstractVector{TT},
    quadrature_method::QuadratureMethod,
    work_weights::AbstractVector{T},
    work_integrand::AbstractVector{TT},
    ::AbstractVector{T}
) where {T<:Real,S<:RealSpace,TT<:Real}
    mesh = quantity.r
    r = collect(mesh)
    r²f_r = quantity.f

    work_weights = @view work_weights[eachindex(r)]
    work_integrand = @view work_integrand[eachindex(r)]

    integration_weights!(work_weights, mesh, quadrature_method)

    f_q = map(qs) do q
        work_integrand .= r²f_r .* fast_sphericalbesselj0.(q .* r)
        integral = dot(work_weights, work_integrand)
        return 4TT(π) * integral
    end

    return NumericDensity{TT,FourierSpace}(qs, f_q)
end

function cubic_spline(
    r::AbstractVector{T},
    quantity::NumericDensity{T,S},
) where {T<:Real,S<:EvaluationSpace}
    f = quantity.(r)
    return NumericDensity{T,S}(r, f, quantity.itp)
end
