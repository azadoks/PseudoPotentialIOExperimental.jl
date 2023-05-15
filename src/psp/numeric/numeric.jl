abstract type NumericPsP{T,S} <: AbstractPsP{T,S} end

identifier(psp::NumericPsP)::String = psp.identifier
element(psp::NumericPsP)::PeriodicTable.Element = PeriodicTable.elements[psp.Zatom]
max_angular_momentum(psp::NumericPsP)::Int = psp.lmax
valence_charge(psp::NumericPsP) = psp.Zval
atomic_charge(psp::NumericPsP) = psp.Zatom
has_spin_orbit(::NumericPsP)::Bool = false  # This is a current limitation

get_quantity(psp::NumericPsP, ::LocalPotential, args...) = psp.Vloc
get_quantity(psp::NumericPsP, ::ValenceChargeDensity, args...) = psp.ρval
get_quantity(psp::NumericPsP, ::CoreChargeDensity, args...) = psp.ρcore
get_quantity(psp::NumericPsP, ::ChiProjector) = psp.χ
get_quantity(psp::NumericPsP, ::BetaProjector) = psp.β
get_quantity(psp::NumericPsP, q::ProjectorFlag, l) = get_quantity(psp, q)[l]
get_quantity(psp::NumericPsP, q::ProjectorFlag, l, n) = get_quantity(psp, q, l)[n]

function n_radials(psp::NumericPsP, quantity::ProjectorFlag, l)
    return length(get_quantity(psp, quantity, l))
end
