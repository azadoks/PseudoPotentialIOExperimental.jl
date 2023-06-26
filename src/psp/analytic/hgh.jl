@doc raw"""
Analytical Hartwigsen-Goedecker-Hutter pseudopotential.

C. Hartwigsen, S. Goedecker, and J. Hutter.
*Pys. Rev. B* **58**, 3641 (1998)](https://doi.org/10.1103/PhysRevB.58.3641)
"""
struct HghPsP{T<:Real,S<:EvaluationSpace} <: AnalyticPsP{T,S}
    "Identifier"
    identifier::String
    "Atomic charge"
    Zatom::Union{Nothing,Int}
    "Valence charge"
    Zval::Int
    "Maximum angular momentum"
    lmax::Int
    "Local part of the potential"
    Vloc::HghLocalPotential{T,S}
    "Nonlocal projectors β[l][n]"
    β::OffsetVector{Vector{HghProjector{T,S}},Vector{Vector{HghProjector{T,S}}}}
    "Nonlocal projector coupling coefficients `D[l][n,m]`"
    D::OffsetVector{Matrix{T},Vector{Matrix{T}}}
end

function HghPsP(file::HghFile)
    el = element(file)
    Zatom = isnothing(el) ? nothing : el.number
    Zval = sum(Float64.(file.zion))

    cloc = file.cloc
    length(cloc) <= 4 || error("length(cloc) > 4 not supported.")
    if length(cloc) < 4
        n_extra = 4 - length(cloc)
        cloc = [cloc; zeros(Float64, n_extra)]
    end

    Vloc = HghLocalPotential{Float64,RealSpace}(Zval, cloc, file.rloc)

    β = map(0:file.lmax) do l
        map(1:size(file.h[l+1], 1)) do n
            return HghProjector(n, l, file.rp[l+1])
        end
    end
    β = OffsetVector(β, 0:file.lmax)  # Ha / √a₀

    D = OffsetVector(file.h, 0:(file.lmax))  # 1 / Ha

    return HghPsP{Float64,RealSpace}(file.identifier, Zatom, Zval, file.lmax, Vloc, β, D)
end

identifier(psp::HghPsP)::String = psp.identifier
function element(psp::HghPsP)
    return isnothing(psp.Zatom) ? nothing : PeriodicTable.elements[Int(psp.Zatom)]
end
max_angular_momentum(psp::HghPsP)::Int = psp.lmax
valence_charge(psp::HghPsP) = psp.Zval
atomic_charge(psp::HghPsP) = psp.Zatom
has_spin_orbit(::HghPsP)::Bool = false
is_norm_conserving(::HghPsP)::Bool = true
is_ultrasoft(::HghPsP)::Bool = false
is_paw(::HghPsP)::Bool = false

has_quantity(::HghPsP, ::AtomicQuantityFlag) = false
has_quantity(::HghPsP, ::LocalPotential) = true
has_quantity(::HghPsP, ::NonLocalProjector) = true
has_quantity(::HghPsP, ::NonLocalCoupling) = true

get_quantity(psp::HghPsP, ::NonLocalCoupling) = psp.D
get_quantity(psp::HghPsP, ::NonLocalCoupling, l) = psp.D[l]
get_quantity(psp::HghPsP, ::NonLocalCoupling, l, n) = psp.D[l][n, n]
get_quantity(psp::HghPsP, ::NonLocalCoupling, l, n, m) = psp.D[l][n, m]
get_quantity(psp::HghPsP, ::NonLocalProjector) = psp.β
get_quantity(psp::HghPsP, ::NonLocalProjector, l) = psp.β[l]
get_quantity(psp::HghPsP, ::NonLocalProjector, l, n) = psp.β[l][n]
get_quantity(psp::HghPsP, ::LocalPotential) = psp.Vloc

n_radials(psp::HghPsP, ::NonLocalProjector, l::Int) = length(psp.β[l])
n_radials(::HghPsP, ::PseudoState, ::Int) = 0

function hankel_transform(psp::HghPsP{T,RealSpace})::HghPsP{T,FourierSpace} where {T<:Real}
    Vloc = hankel_transform(get_quantity(psp, LocalPotential()))

    β = map(get_quantity(psp, NonLocalProjector())) do βl
        map(βl) do βln
            hankel_transform(βln)
        end
    end

    return HghPsP{T,FourierSpace}(psp.identifier, psp.Zatom, psp.Zval, psp.lmax, Vloc, β, psp.D)
end

function energy_correction(TT::Type{<:Real}, psp::HghPsP{T,RealSpace})::TT where {T}
    return energy_correction(TT, get_quantity(psp, LocalPotential()))
end

function energy_correction(psp::HghPsP{T,RealSpace})::T where {T}
    return energy_correction(T, get_quantity(psp, LocalPotential()))
end
