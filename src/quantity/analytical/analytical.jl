abstract type AnalyticPsPQuantity{T,S} <: AbstractPsPQuantity end

abstract type AnalyticalLocalPotential{T,S} <: AnalyticPsPQuantity{T,S} end
abstract type AnalyticalProjector{T,S} <: AnalyticPsPQuantity{T,S} end
abstract type AnalyticalDensity{T,S} <: AnalyticPsPQuantity{T,S} end
