struct HghLocalPotential{T<:Real,S<:EvaluationSpace} <: AnalyticalLocalPotential{T,S}
    Zval::T
    cloc::Vector{T}
    rloc::T
end

# [GTH98] (1)
function (Vloc::HghLocalPotential{T,S})(r::TT) where {T<:Real,S<:RealSpace,TT<:Real}
    r::TT += iszero(r) ? eps(TT) : zero(TT)  # quick hack for the division by zero below
    rr::TT = r / Vloc.rloc
    cloc = Vloc.cloc
    return -Vloc.Zval / r * erf(rr / sqrt(TT(2))) + exp(-rr^2 / 2) * (cloc[1] + cloc[2] * rr^2 + cloc[3] * rr^4 + cloc[4] * rr^6)
end

function hankel_transform(Vloc::HghLocalPotential{T,S}) where {T<:Real,S<:RealSpace}
    return HghLocalPotential{T,FourierSpace}(Vloc.Zval, Vloc.cloc, Vloc.rloc)
end

@doc raw"""
The local potential of a HGH pseudopotentials in reciprocal space
can be brought to the form ``Q(t) / (t^2 exp(t^2 / 2))``
where ``x = r_\text{loc} q`` and `Q`
is a polynomial of at most degree 8. This function returns `Q`.
"""
@inline function _local_potential_polynomial_fourier(Vloc::HghLocalPotential{T,S}, x::TT) where {T<:Real,S<:FourierSpace,TT<:Real}
    rloc::TT = Vloc.rloc
    Zval::TT = Vloc.Zval

    # The polynomial prefactor P(t) (as used inside the { ... } brackets of equation
    # (5) of the HGH98 paper)
    P = (Vloc.cloc[1]
         + Vloc.cloc[2] * (3 - x^2)
         + Vloc.cloc[3] * (15 - 10x^2 + x^4)
         + Vloc.cloc[4] * (105 - 105x^2 + 21x^4 - x^6))

    return 4TT(π) * rloc^2 * (-Zval + sqrt(TT(π) / 2) * rloc * x^2 * P)
end

# [GTH98] (6) except they do it with plane waves normalized by 1/sqrt(Ω).
function (Vloc::HghLocalPotential{T,S})(q::TT) where {T<:Real,S<:FourierSpace,TT<:Real}
    x::TT = q * Vloc.rloc
    return _local_potential_polynomial_fourier(Vloc, x) * exp(-x^2 / 2) / x^2
end


struct HghProjector{T<:Real,S<:EvaluationSpace} <: AnalyticalProjector{T,S}
    n::Int
    l::Int
    rnl::T
end

function HghProjector(n::Int, l::Int, rnl::T) where {T<:Real}
    return HghProjector{T,RealSpace}(n, l, rnl)
end

function (β::HghProjector{T,S})(r::TT) where {T<:Real,S<:RealSpace,TT<:Real}
    n = β.n
    l = β.l
    rnl::TT = β.rnl
    ired::TT = (4n - 1) / TT(2)
    return sqrt(TT(2)) * r^(l + 2(n - 1)) * exp(-r^2 / 2rnl^2) / rnl^(l + ired) /
           sqrt(gamma(l + ired))
end

function hankel_transform(β::HghProjector{T,S}) where {T<:Real,S<:RealSpace}
    return HghProjector{T,FourierSpace}(β.n, β.l, β.rnl)
end

@doc raw"""
The nonlocal projectors of a HGH pseudopotentials in reciprocal space
can be brought to the form ``Q(t) e^{-t^2 / 2}`` where ``t = r_l q``
and `Q` is a polynomial. This function returns `Q`.
"""
@inline function _beta_projector_polynomial_fourier(β::HghProjector{T,S}, x::TT) where {T<:Real,S<:FourierSpace,TT<:Real}
    n = β.n
    l = β.l
    rnl::TT = β.rnl
    common::TT = 4TT(π)^(5 / TT(4)) * sqrt(TT(2^(l + 1)) * rnl^3)

    # Note: In the (l == 0 && i == 2) case the HGH paper has an error.
    #       The first 8 in equation (8) should not be under the sqrt-sign
    #       This is the right version (as shown in the GTH paper)
    (l == 0 && n == 1) && return convert(TT, common)
    (l == 0 && n == 2) && return common * 2 / sqrt(TT(15)) * (3 - x^2)
    (l == 0 && n == 3) && return common * 4 / 3sqrt(TT(105)) * (15 - 10x^2 + x^4)
    #
    (l == 1 && n == 1) && return common * 1 / sqrt(TT(3)) * x
    (l == 1 && n == 2) && return common * 2 / sqrt(TT(105)) * x * (5 - x^2)
    (l == 1 && n == 3) && return common * 4 / 3sqrt(TT(1155)) * x * (35 - 14x^2 + x^4)
    #
    (l == 2 && n == 1) && return common * 1 / sqrt(TT(15)) * x^2
    (l == 2 && n == 2) && return common * 2 / 3sqrt(TT(105)) * x^2 * (7 - x^2)
    #
    (l == 3 && n == 1) && return common * 1 / sqrt(TT(105)) * x^3

    throw(ArgumentError("Not implemented for l=$l and i=$n"))
end

function(β::HghProjector{T,S})(q::TT) where {T<:Real,S<:FourierSpace,TT<:Real}
    x::TT = q * β.rnl
    return _beta_projector_polynomial_fourier(β, x) * exp(-x^2 / TT(2))
end
