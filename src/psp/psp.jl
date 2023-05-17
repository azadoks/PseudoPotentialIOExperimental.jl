abstract type AbstractPsP{T,S} end

function identifier(::AbstractPsP) end
function element(::AbstractPsP) end
function max_angular_momentum(::AbstractPsP) end
function valence_charge(::AbstractPsP) end
function atomic_charge(::AbstractPsP) end
function has_spin_orbit(::AbstractPsP) end
function is_norm_conserving(::AbstractPsP) end
function is_ultrasoft(::AbstractPsP) end
function is_paw(::AbstractPsP) end

function has_quantity(::AbstractPsP, ::PsPQuantityFlag, args...) end
function get_quantity(::AbstractPsP, ::PsPQuantityFlag, args...) end
function n_radials(::AbstractPsP, ::ProjectorFlag, l) end
function hankel_tansform(::AbstractPsP, args...; kwargs...) end
function energy_correction(::AbstractPsP, args...; kwargs...) end

"""
Angular momenta values of the pseudopotential.
"""
angular_momenta(psp::AbstractPsP) = 0:max_angular_momentum(psp)

"""
Type of relativistic treatment (fully relativistic or scalar-relativistic).
"""
relativistic_treatment(psp::AbstractPsP)::Symbol = has_spin_orbit(psp) ? :full : :scalar

"""
Formalism of the pseudopotential (norm-conserving, ultrasoft, projector-augmented wave,
or Coulomb).
"""
function formalism(psp::AbstractPsP)::Type
    is_norm_conserving(psp) && return NormConservingPsP
    is_ultrasoft(psp) && return UltrasoftPsP
    is_paw(psp) && return ProjectorAugmentedWavePsP
end

"""
The number of radial parts corresponding to a quantity at a given angular
momentum, or the total number of radial parts for all angular momenta for
a given quantity.
"""
function n_radials(psp::AbstractPsP, q::ProjectorFlag)
    return sum(l -> n_radials(psp, q, l), angular_momenta(psp); init=0)
end

"""
The number of radial + angular parts corresponding to a quantity at a given angular
momentum, or the total number of radial + angular parts for all angular momenta for
a given quantity.

i.e., count the number of combinations for the quantum numbers ``l`` and ``m`` up to the
maximum angular momentum of the pseudopotential, accounting for multi-projector
pseudopotentials that provide multiple radial parts at an individual angular momentum.
"""
function n_angulars(psp::AbstractPsP, quantity::ProjectorFlag, l)
    return (2l + 1) * n_radials(psp, quantity, l)
end
function n_angulars(psp::AbstractPsP, quantity::ProjectorFlag)
    return sum(l -> n_angulars(psp, quantity, l), angular_momenta(psp); init=0)
end

maximum_radius(::Nothing) = -Inf
maximum_radius(psp::AbstractPsP, q::PsPQuantityFlag, args...) = maximum_radius(get_quantity(psp, q, args...))
function maximum_radius(psp::AbstractPsP, q::ProjectorFlag, l)
    return maximum(maximum_radius(psp, q, l, n) for n in 1:n_radials(psp, q, l); init=-Inf)
end
function maximum_radius(psp::AbstractPsP, q::ProjectorFlag)
    return maximum(maximum_radius(psp, q, l) for l in angular_momenta(psp); init=-Inf)
end
function maximum_radius(psp::AbstractPsP)
    quantities = (LocalPotential(), ValenceChargeDensity(), CoreChargeDensity(),
                  ChiProjector(), BetaProjector())
    return maximum(maximum_radius(psp, q) for q in quantities; init=-Inf)
end

minimum_radius(::Nothing) = Inf
minimum_radius(psp::AbstractPsP, q::PsPQuantityFlag, args...) = minimum_radius(get_quantity(psp, q, args...))
function minimum_radius(psp::AbstractPsP, q::ProjectorFlag, l)
    return minimum(minimum_radius(psp, q, l, n) for n in 1:n_radials(psp, q, l); init=Inf)
end
function minimum_radius(psp::AbstractPsP, q::ProjectorFlag)
    return minimum(minimum_radius(psp, q, l) for l in angular_momenta(psp); init=Inf)
end
function minimum_radius(psp::AbstractPsP)
    quantities = (LocalPotential(), ValenceChargeDensity(), CoreChargeDensity(),
                  ChiProjector(), BetaProjector())
    return minimum(minimum_radius(psp, q) for q in quantities; init=Inf)
end

function extrema_radii(psp::AbstractPsP, q::PsPQuantityFlag, args...)
    return (minimum_radius(psp, q, args...), maximum_radius(psp, q, args...))
end

function CoulombCorrection(psp::AbstractPsP{T,S}) where {T<:Real,S<:EvaluationSpace}
    return CoulombCorrection{T,S}(valence_charge(psp))
end

function ErfCoulombCorrection(psp::AbstractPsP{T,S}) where {T<:Real,S<:EvaluationSpace}
    return ErfCoulombCorrection{T,S}(valence_charge(psp))
end

Base.Broadcast.broadcastable(psp::AbstractPsP) = Ref(psp)

function Base.show(io::IO, psp::AbstractPsP)
    typename = string(typeof(psp))
    el = element(psp)
    z = valence_charge(psp)
    return print(io, "$typename(element=$el, z_valence=$z)")
end

function Base.show(io::IO, ::MIME"text/plain", psp::AbstractPsP)
    println(io, typeof(psp))
    @printf "%032s: %s\n" "identifier" identifier(psp)
    @printf "%032s: %s\n" "formalism" formalism(psp)
    @printf "%032s: %s\n" "element" element(psp)
    @printf "%032s: %f\n" "valence charge" valence_charge(psp)
    @printf "%032s: %s\n" "relativistic treatment" relativistic_treatment(psp)
    @printf "%032s: %s\n" "non-linear core correction" has_quantity(psp,
                                                                    CoreChargeDensity())
    @printf "%032s: %d\n" "maximum angular momentum" max_angular_momentum(psp)
    @printf "%032s: %s\n" "number of beta projectors" n_radials(psp, BetaProjector())
    @printf "%032s: %s" "number of chi projectors" n_radials(psp, ChiProjector())
end
