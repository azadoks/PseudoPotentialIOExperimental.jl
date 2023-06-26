abstract type NumericPsP{T,S} <: AbstractPsP{T,S} end

identifier(psp::NumericPsP)::String = psp.identifier
element(psp::NumericPsP)::PeriodicTable.Element = PeriodicTable.elements[round(Int, psp.Zatom)]
max_angular_momentum(psp::NumericPsP)::Int = psp.lmax
valence_charge(psp::NumericPsP) = psp.Zval
atomic_charge(psp::NumericPsP) = psp.Zatom
has_spin_orbit(::NumericPsP)::Bool = false  # This is a current limitation

has_quantity(psp::NumericPsP, q::AtomicQuantityFlag, args...) = !isnothing(get_quantity(psp, q, args...))

get_quantity(psp::NumericPsP, ::LocalPotential, args...) = psp.Vloc
get_quantity(psp::NumericPsP, ::PseudoValenceDensity, args...) = psp.ρval
get_quantity(psp::NumericPsP, ::CoreDensity, args...) = psp.ρcore
get_quantity(psp::NumericPsP, ::PseudoState) = psp.χ
get_quantity(psp::NumericPsP, ::NonLocalProjector) = psp.β
get_quantity(psp::NumericPsP, q::ProjectorFlag, l) = get_quantity(psp, q)[l]
get_quantity(psp::NumericPsP, q::ProjectorFlag, l, n) = get_quantity(psp, q, l)[n]
get_quantity(psp::NumericPsP, ::NonLocalCoupling) = psp.D
get_quantity(psp::NumericPsP, ::NonLocalCoupling, l) = psp.D[l]
get_quantity(psp::NumericPsP, ::NonLocalCoupling, l, n) = psp.D[l][n,n]
get_quantity(psp::NumericPsP, ::NonLocalCoupling, l, n, m) = psp.D[l][n,m]

function n_radials(psp::NumericPsP, quantity::ProjectorFlag, l)
    return length(get_quantity(psp, quantity, l))
end

function energy_correction(TT::Type{<:Real}, psp::NumericPsP{T,RealSpace}; kwargs...)::TT where {T}
    return energy_correction(TT, get_quantity(psp, LocalPotential()), kwargs...)
end

function energy_correction(psp::NumericPsP{T,RealSpace}; kwargs...)::T where {T}
    return energy_correction(get_quantity(psp, LocalPotential()); kwargs...)
end
