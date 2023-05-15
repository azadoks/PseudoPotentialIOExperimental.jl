"""
Type representing a numeric projector-augmented wave pseudopotential.
"""
struct ProjectorAugmentedWavePsP{T,S} <: NumericPsP{T,S} end

function ProjectorAugmentedWavePsP(::PsPFile)
    error("Projecter-augmented wave pseudopotentials are not implemented")
end
