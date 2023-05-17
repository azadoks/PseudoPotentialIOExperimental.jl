# TODO: interpolate_onto and hankel_transform are heavily duplicated

abstract type NumericPsPQuantity{T<:Real,S<:EvaluationSpace} <: AbstractPsPQuantity end

minimum_radius(q::NumericPsPQuantity) = minimum(q.r; init=Inf)
maximum_radius(q::NumericPsPQuantity) = maximum(q.r; init=0)
extrema_radii(q::NumericPsPQuantity) = (minimum_radius(q), maximum_radius(q))

hankel_transform(::Nothing, args...; kwargs...)::Nothing = nothing
function hankel_transform(quantity::NumericPsPQuantity{T,RealSpace},
                          qs::AbstractVector{TT}=range(T(0), T(30, 3001)); kwargs...) where {T,TT<:Real}
    n = length(quantity.r)
    work_weights = Vector{T}(undef, n)
    work_integrand = Vector{TT}(undef, n)
    work_f = Vector{T}(undef, n)
    return hankel_transform(quantity, qs, work_weights, work_integrand, work_f; kwargs...)
end

interpolate_onto(::Nothing, ::Any)::Nothing = nothing
function interpolate_onto(q::NumericPsPQuantity{T,S},
                          maximum_spacing::T)::NumericPsPQuantity{T,S} where {T<:Real,S<:EvaluationSpace}
    r = UniformMesh(extrema_radii(q)..., maximum_spacing)
    return interpolate_onto(q, r)
end

function (quantity::NumericPsPQuantity{T,S})(x::TT)::TT where {T<:Real,S<:EvaluationSpace,TT<:Real}
    return quantity.itp(x)
end

struct NumericLocalPotential{T<:Real,
                             S<:EvaluationSpace} <: NumericPsPQuantity{T,S}
    Zval::T
    r::AbstractVector{T}
    f::AbstractVector{T}
    itp::BSplineKit.SplineInterpolation
end

function (quantity::NumericLocalPotential{T,FourierSpace})(x::TT)::TT where {T<:Real,TT<:Real}
    !iszero(x) && return quantity.itp(x)
    return TT(-Inf)
end

function NumericLocalPotential{T,S}(Zval::T, r::AbstractVector{T}, f::AbstractVector{T};
                                    spline_order=BSplineOrder(4)) where {T<:Real,S<:EvaluationSpace}
    # The local potential at q==0 is undefined; the value will be garbage and will mess up
    # the interpolation, so we cut it out.
    itp = BSplineKit.interpolate(r[r .> 0], f[r .> 0], spline_order, Natural())
    return NumericLocalPotential{T,S}(Zval, r, f, itp)
end

function NumericLocalPotential(Zval::T, r::AbstractVector{T}, f::AbstractVector{T}; kwargs...) where {T<:Real}
    return NumericLocalPotential{T,RealSpace}(Zval, r, f; kwargs...)
end

function hankel_transform(quantity::NumericLocalPotential{T,RealSpace},
                          qs::AbstractVector{TT},
                          work_weights::AbstractVector{T},
                          work_integrand::AbstractVector{TT},
                          work_f::AbstractVector{T};
                          quadrature_method=Simpson(),
                          local_potential_correction=ErfCoulombCorrection(quantity.Zval)) where {T<:Real,TT<:Real}
    mesh = quantity.r
    r = collect(mesh)
    f_r = quantity.f
    work_f .= r .* f_r .- local_potential_correction.(r)

    work_weights = @view work_weights[eachindex(r)]
    work_integrand = @view work_integrand[eachindex(r)]

    integration_weights!(work_weights, mesh, quadrature_method)
    correction_q = hankel_transform(local_potential_correction)

    f_q = map(qs) do q
        work_integrand .= work_f .* sin.(q .* r)
        integral = dot(work_weights, work_integrand)
        return 4TT(π) * (integral / q + correction_q(q))
    end

    return NumericLocalPotential{TT,FourierSpace}(quantity.Zval, qs, f_q)
end

function interpolate_onto(Vloc_r::NumericLocalPotential{T,S},
                          r::AbstractVector{T}) where {T<:Real,S<:EvaluationSpace}
    f = Vloc_r.(r)
    return NumericLocalPotential{T,S}(Vloc_r.Zval, r, f, Vloc_r.itp)
end

function energy_correction(TT::Type{<:Real}, Vloc::NumericLocalPotential{T,RealSpace};
                           quadrature_method::QuadratureMethod=Simpson(),
                           kwargs...)::TT where {T<:Real}
    integrand = Vloc.r .* (Vloc.r .* Vloc.f .+ Vloc.Zval)
    weights = integration_weights(Vloc.r, quadrature_method)
    return 4TT(π) * dot(weights, integrand)
end

function energy_correction(Vloc::NumericLocalPotential{T,RealSpace}; kwargs...)::T where {T<:Real}
    return energy_correction(T, Vloc; kwargs...)
end

struct NumericProjector{T<:Real,
                        S<:EvaluationSpace} <: NumericPsPQuantity{T,S}
    n::Int
    l::Int
    j::T
    r::AbstractVector{T}
    f::AbstractVector{T}
    itp::BSplineKit.SplineInterpolation
end

function NumericProjector{T,S}(n::Int, l::Int, j::T, r::AbstractVector{T}, f::AbstractVector{T};
                               spline_order=BSplineOrder(4)) where {T<:Real,S<:EvaluationSpace}
    itp = BSplineKit.interpolate(r, f, spline_order, Natural())
    return NumericProjector{T,S}(n, l, j, r, f, itp)
end

function NumericProjector(n::Int, l::Int, j::T, r::AbstractVector{T}, f::AbstractVector{T};
                          kwargs...) where {T<:Real}
    return NumericProjector{T,RealSpace}(n, l, j, r, f; kwargs...)
end

function hankel_transform(quantity::NumericProjector{T,RealSpace},
                          qs::AbstractVector{TT},
                          work_weights::AbstractVector{T},
                          work_integrand::AbstractVector{TT},
                          ::AbstractVector{T};
                          kwargs...)::NumericProjector{TT,FourierSpace} where {T<:Real,TT<:Real}
    f_q = hankel_transform(quantity.f, quantity.r, qs, work_weights, work_integrand, quantity.l; kwargs...)

    return NumericProjector{TT,FourierSpace}(quantity.n, quantity.l, quantity.j, qs, f_q)
end

function interpolate_onto(quantity::NumericProjector{T,S},
                          r::AbstractVector{T}) where {T<:Real,S<:EvaluationSpace}
    f = quantity.(r)
    return NumericProjector{T,S}(quantity.n, quantity.l, quantity.j, r, f, quantity.itp)
end

struct NumericDensity{T<:Real,
                      S<:EvaluationSpace} <: NumericPsPQuantity{T,S}
    r::AbstractVector{T}
    f::AbstractVector{T}
    itp::BSplineKit.SplineInterpolation
end

function NumericDensity{T,S}(r::AbstractVector{T}, f::AbstractVector{T};
                             spline_order=BSplineOrder(4)) where {T<:Real,S<:EvaluationSpace}
    itp = BSplineKit.interpolate(r, f, spline_order, Natural())
    return NumericDensity{T,S}(r, f, itp)
end

function NumericDensity(r::AbstractVector{T}, f::AbstractVector{T}; kwargs...) where {T<:Real}
    return NumericDensity{T,RealSpace}(r, f; kwargs...)
end

function hankel_transform(quantity::NumericDensity{T,RealSpace}, qs::AbstractVector{TT},
                          work_weights::AbstractVector{T}, work_integrand::AbstractVector{TT},
                          ::AbstractVector{T};
                          kwargs...)::NumericDensity{TT,FourierSpace} where {T<:Real,TT<:Real}
    f_q = hankel_transform(quantity.f, quantity.r, qs, work_weights, work_integrand; kwargs...)
    return NumericDensity{TT,FourierSpace}(qs, f_q)
end

function interpolate_onto(quantity::NumericDensity{T,S},
                          r::AbstractVector{T}) where {T<:Real,S<:EvaluationSpace}
    f = quantity.(r)
    return NumericDensity{T,S}(r, f, quantity.itp)
end

struct NumericAugmentation{T<:Real,
                           S<:EvaluationSpace} <: NumericPsPQuantity{T,S}
    n::Int
    m::Int
    l::Int
    r::AbstractVector{T}
    f::AbstractVector{T}
    itp::BSplineKit.SplineInterpolation
end

function NumericAugmentation{T,S}(n::Int, m::Int, l::Int, r::AbstractVector{T}, f::AbstractVector{T};
                                  spline_order=BSplineOrder(4)) where {T<:Real,S<:EvaluationSpace}
    itp = BSplineKit.interpolate(r, f, spline_order, Natural())
    return NumericAugmentation{T,S}(n, m, l, r, f, itp)
end

function NumericAugmentation(n::Int, m::Int, l::Int, r::AbstractVector{T}, f::AbstractVector{T};
                             kwargs...) where {T<:Real}
    return NumericAugmentation{T,RealSpace}(n, m, l, r, f; kwargs...)
end

function hankel_transform(quantity::NumericAugmentation{T,RealSpace},
                          qs::AbstractVector{TT},
                          work_weights::AbstractVector{T},
                          work_integrand::AbstractVector{TT},
                          work_f::AbstractVector{T};
                          kwargs...)::NumericAugmentation{TT,FourierSpace} where {T<:Real,TT<:Real}
    work_f = @view work_f[eachindex(quantity.f)]
    work_f .= quantity.r .^ 2 .* quantity.f
    f_q = hankel_transform(quantity.f, quantity.r, qs, work_weights, work_integrand, quantity.l; kwargs...)

    return NumericAugmentation{TT,FourierSpace}(quantity.n, quantity.m, quantity.l, qs, f_q)
end

function interpolate_onto(quantity::NumericAugmentation{T,S},
                          r::AbstractVector{T}) where {T<:Real,S<:EvaluationSpace}
    f = quantity.(r)
    return NumericAugmentation{T,S}(quantity.n, quantity.m, quantity.l, r, f, quantity.itp)
end
