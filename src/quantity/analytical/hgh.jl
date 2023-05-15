struct HghLocalPotential{T<:Real,S<:EvaluationSpace} <: AnalyticalLocalPotential{T,S}
    Zval::T
    cloc::Vector{T}
    rloc::T
end

function (Vloc::HghLocalPotential{T,S})(r::R) where {T<:Real,S<:RealSpace,R<:Real}
    r::R += iszero(r) ? eps(R) : zero(R)
    rr::R = r / Vloc.rloc
    cloc = Vloc.cloc
    return -Vloc.Zval / r * erf(rr / sqrt(R(2))) + exp(-rr^2 / 2) * (cloc[1] + cloc[2] * rr^2 + cloc[3] * rr^4 + cloc[4] * rr^6)
end

function hankel_transform(Vloc::HghLocalPotential{T,S}) where {T<:Real,S<:RealSpace}
    return HghLocalPotential{T,FourierSpace}(Vloc.Zval, Vloc.cloc, Vloc.rloc)
end

function (Vloc::HghLocalPotential{T,S})(q::F) where {T<:Real,S<:FourierSpace,F<:Real}
    rloc::F = Vloc.rloc
    Zval::F = Vloc.Zval
    x::F = q * rloc
    cloc = Vloc.cloc

    P = (cloc[1] + cloc[2] * (3 - x^2) + cloc[3] * (15 - 10x^2 + x^4) + cloc[4] * (105 - 105x^2 + 21x^4 - x^6))
    P = 4F(π) * rloc^2 * (-Zval + sqrt(T(π) / 2) * rloc * x^2 * P)
    return exp(-x^2 / 2) / x^2 * P
end

struct HghProjector{T<:Real,S<:EvaluationSpace} <: AnalyticalProjector{T,S}
    n::Int
    l::Int
    rnl::T
end

function (β::HghProjector{T,S})(r::TT) where {T<:Real,S<:RealSpace,TT<:Real}
    rnl::TT = β.rnl
end

function hankel_transform(β::HghProjector{T,S}) where {T<:Real,S<:RealSpace}
    return HghProjector{T,FourierSpace}(β.rnl)
end

function(β::HghProjector{T,S})(q::TT) where {T<:Real,S<:FourierSpace,TT<:Real}
end
