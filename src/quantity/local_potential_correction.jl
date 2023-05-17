abstract type LocalPotentialCorrection{T,S} end
Base.Broadcast.broadcastable(corr::LocalPotentialCorrection) = Ref(corr)


struct CoulombCorrection{T<:Real,S<:EvaluationSpace} <: LocalPotentialCorrection{T,S}
    Z::T
end

function CoulombCorrection(Z::T) where {T<:Real}
    return CoulombCorrection{T,RealSpace}(Z)
end

function (corr::CoulombCorrection{T,S})(_::TT)::TT where {T<:Real,TT<:Real,S<:RealSpace}
    return -corr.Z
end

function (corr::CoulombCorrection{T,S})(q::TT)::TT where {T<:Real,TT<:Real,S<:FourierSpace}
    return -corr.Z / q^2
end

function hankel_transform(corr::CoulombCorrection{T,S}) where {T<:Real,S<:RealSpace}
    return CoulombCorrection{T,FourierSpace}(corr.Z)
end

struct ErfCoulombCorrection{T<:Real,S<:EvaluationSpace} <: LocalPotentialCorrection{T,S}
    Z::T
end

function ErfCoulombCorrection(Z::T) where {T<:Real}
    return ErfCoulombCorrection{T,RealSpace}(Z)
end

function (corr::ErfCoulombCorrection{T,S})(r::TT)::TT where {T<:Real,TT<:Real,S<:RealSpace}
    return -corr.Z * erf(r)
end

function (corr::ErfCoulombCorrection{T,S})(q::TT)::TT where {T<:Real,TT<:Real,S<:FourierSpace}
    return -corr.Z * exp(-q^2 / 4) / q^2
end

function hankel_transform(corr::ErfCoulombCorrection{T,S}) where {T<:Real,S<:RealSpace}
    return ErfCoulombCorrection{T,FourierSpace}(corr.Z)
end
