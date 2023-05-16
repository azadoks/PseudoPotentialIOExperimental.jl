"""
Type representing a numeric ultrasoft pseudopotential.
"""
struct UltrasoftPsP{T,S} <: NumericPsP{T,S}
    "Identifier"
    identifier::String
    "Total charge."
    Zatom::Int
    "Valence charge."
    Zval::Int
    "Maximum angular momentum."
    lmax::Int
    "Local part of the potential on the radial mesh (without an r² prefactor)."
    Vloc::NumericLocalPotential{T,S}
    "Nonlocal projectors β[l][n] on the radial mesh (with an r² prefactor)."
    β::OffsetVector{Vector{NumericProjector{T,S}},Vector{Vector{NumericProjector{T,S}}}}
    "Projector coupling coefficients D[l][n,m]."
    D::OffsetVector{Matrix{T},Vector{Matrix{T}}}
    "Pseudo-atomic wavefunctions χ[l][n] on the radial mesh (with an r² prefactor)."
    χ::OffsetVector{Vector{NumericProjector{T,S}},Vector{Vector{NumericProjector{T,S}}}}
    "Augmentation charge density functions Q[l][n,m] on the radial mesh"
    Q::OffsetVector{Matrix{NumericAugmentation{T,S}},Vector{Matrix{NumericAugmentation{T,S}}}}
    "Augmentation charges q[l][n,m]"
    q::OffsetVector{Matrix{T},Vector{Matrix{T}}}
    "Model core charge density on the radial mesh (with an r² prefactor)."
    ρcore::Union{Nothing,NumericDensity{T,S}}
    "Valence charge density on the radial mesh (with an r² prefactor)."
    ρval::Union{Nothing,NumericDensity{T,S}}
end

function UltrasoftPsP(upf::UpfFile)
    if !in(upf.header.pseudo_type, ("US", "USPP"))
        error("Provided `UpfFile` is not an ultrasoft pseudo")
    end
    if (upf.header.relativistic == "full") | upf.header.has_so | (!isnothing(upf.spin_orb))
        error("Fully relativistic pseudos are not supported")
    end
    return _upf_construct_us_internal(upf)
end

function _upf_construct_augmentation_q_with_l(upf::UpfFile)
    Q = OffsetVector([Matrix{NumericAugmentation{Float64,RealSpace}}(undef, upf.header.number_of_proj,
                                                                     upf.header.number_of_proj)
                      for l in 0:(2upf.header.l_max)], 0:(2upf.header.l_max))
    for l in 0:(2upf.header.l_max)
        # Fill Q with zeroed quantities
        for i in 1:(upf.header.number_of_proj), j in 1:(upf.header.number_of_proj)
            r = ArbitraryMesh(upf.mesh.r, upf.mesh.rab)
            f = zero(r)
            Q[l][i, j] = NumericAugmentation(i, j, l, r, f)
        end
        # Replace the zeroed quantities with data from the UPF where present
        Q_upf_l = filter(qijl -> qijl.angular_momentum == l,
                         upf.nonlocal.augmentation.qijls)
        for Q_upf in Q_upf_l
            n1 = Q_upf.first_index
            n2 = Q_upf.second_index
            f = Q_upf.qijl
            r = ArbitraryMesh(upf.mesh.r[eachindex(f)], upf.mesh.rab[eachindex(f)])

            Q[l][n1, n2] = NumericAugmentation(n1, n2, l, r, f)
            Q[l][n2, n1] = NumericAugmentation(n2, n1, l, r, f)
        end
    end
    return Q
end

@views function _upf_construct_augmentation_qfcoef(upf::UpfFile)
    #TODO check correctness
    r = upf.mesh.r
    r2 = upf.mesh.r .^ 2
    nqf = upf.nonlocal.augmentation.nqf
    nqlc = 2upf.header.l_max + 1

    Q = OffsetVector([Matrix{NumericAugmentation{Float64,RealSpace}}(undef, upf.header.number_of_proj,
                                                                     upf.header.number_of_proj)
                      for l in 0:(2upf.header.l_max)], 0:(2upf.header.l_max))
    for l in 0:(2upf.header.l_max), i in 1:(upf.header.number_of_proj),
        j in 1:(upf.header.number_of_proj)
        # Fill Q with zero vectors
        r = ArbitraryMesh(upf.mesh.r, upf.mesh.rab)
        f = zero(r)
        Q[l][i, j] = NumericAugmentation(i, j, l, r, f)
    end
    for (Q_upf, Qfcoef_upf) in
        zip(upf.nonlocal.augmentation.qijs, upf.nonlocal.augmentation.qfcoefs)
        # Replace the zero vectors with datat from the UPF where present
        # It's not worth the effort to make these into OffsetVectors zero-indexed for l.
        qfcoef = reshape(Qfcoef_upf.qfcoef, nqf, nqlc)
        rinner = upf.nonlocal.augmentation.rinner

        i = Q_upf.first_index
        j = Q_upf.second_index

        li = upf.nonlocal.betas[i].angular_momentum
        lj = upf.nonlocal.betas[j].angular_momentum

        for l in abs(li - lj):2:(li + lj)
            qij = copy(Q_upf.qij)
            ircut = findfirst(i -> r[i] > rinner[l + 1], eachindex(r)) - 1
            poly = Polynomial(qfcoef[:, l + 1])
            qij[1:ircut] = r[1:ircut] .^ (l + 2) .* poly.(r2[1:ircut])

            n1 = Q_upf.first_index
            n2 = Q_upf.second_index
            r = ArbitraryMesh(upf.mesh.r[1:ircut], upf.mesh.rab[1:ircut])

            Q[l][n1, n2] = NumericAugmentation(n1, n2, l, r, qij)
            Q[l][n2, n1] = NumericAugmentation(n2, n1, l, r, qij)
        end
    end
    return Q
end

function _upf_construct_us_internal(upf::UpfFile)
    nc = _upf_construct_nc_internal(upf)
    # Number of projectors at each angular momentum
    nβ = OffsetVector(length.(nc.β), 0:(nc.lmax))
    # Find the first/last indices in upf.nonlocal.dij for each angular momentum so the
    # sub-arrays q[l][n,m] can be extracted
    cum_nβ = [0, cumsum(nβ)...]
    q_upf = upf.nonlocal.augmentation.q
    q = OffsetVector(map(i -> collect(q_upf[(cum_nβ[i] + 1):cum_nβ[i + 1],
                                            (cum_nβ[i] + 1):cum_nβ[i + 1]]),
                         1:(length(cum_nβ) - 1)), 0:(nc.lmax))
    # Wrangle the agumentation functions. For UPF v2.0.1, they are given with indeices
    # i ∈ 1:nβ, j ∈ 1:nβ, l ∈ 0:2lmax. In the old format, they are given only at i and j
    # with additional coefficients for a polynomial expansion within a cutoff radius
    # These are used to reconstruct the full Qijl.
    if upf.nonlocal.augmentation.q_with_l
        Q = _upf_construct_augmentation_q_with_l(upf)
    elseif upf.nonlocal.augmentation.nqf > 0
        Q = _upf_construct_augmentation_qfcoef(upf)
    else
        error("q_with_l == false and nqf == 0, unsure what to do...")
    end
    return UltrasoftPsP{Float64,RealSpace}(nc.identifier, nc.Zatom, nc.Zval, nc.lmax,
                                           nc.Vloc, nc.β, nc.D, nc.χ, Q, q,
                                           nc.ρcore, nc.ρval)
end

is_norm_conserving(::UltrasoftPsP)::Bool = false
is_ultrasoft(::UltrasoftPsP)::Bool = true
is_paw(::UltrasoftPsP)::Bool = false

# #TODO test the augmentation functions
# has_quantity(psp::UltrasoftPsP, ::AugmentationCoupling) = true
has_quantity(psp::UltrasoftPsP, ::AugmentationFunction) = true
# get_quantity(psp::UltrasoftPsP, ::AugmentationCoupling, l) = psp.q[l]
# get_quantity(psp::UltrasoftPsP, ::AugmentationCoupling, l, n) = psp.q[l][n, n]
# get_quantity(psp::UltrasoftPsP, ::AugmentationCoupling, l, n, m) = psp.q[l][n, m]
get_quantity(psp::UltrasoftPsP, ::AugmentationFunction) = psp.Q
get_quantity(psp::UltrasoftPsP, ::AugmentationFunction, l) = psp.Q[l]
get_quantity(psp::UltrasoftPsP, ::AugmentationFunction, l, n) = psp.Q[l][n, n]
get_quantity(psp::UltrasoftPsP, ::AugmentationFunction, l, n, m) = psp.Q[l][n, m]

function hankel_transform(psp::UltrasoftPsP{T,S};
                          qs::AbstractVector{TT}=range(; start=T(0), stop=T(30), length=3001),
                          quadrature_method=Simpson(),
                          local_potential_correction=CoulombCorrection(psp)
                          )::UltrasoftPsP{T,FourierSpace} where {T<:Real,S<:RealSpace,TT<:Real}
    n = maximum_mesh_length(psp)
    work_weights = Vector{T}(undef, n)
    work_integrand = Vector{TT}(undef, n)
    work_f = Vector{T}(undef, n)

    Vloc = hankel_transform(get_quantity(psp, LocalPotential()), local_potential_correction, qs,
                            quadrature_method, work_weights, work_integrand, work_f)
    β = map(get_quantity(psp, BetaProjector())) do βl
        map(βl) do βln
            return hankel_transform(βln, qs, quadrature_method, work_weights, work_integrand, work_f)
        end
    end
    χ = map(get_quantity(psp, ChiProjector())) do χl
        map(χl) do χln
            return hankel_transform(χln, qs, quadrature_method, work_weights, work_integrand, work_f)
        end
    end
    Q = map(get_quantity(psp, AugmentationFunction())) do Ql
        map(Ql) do Qijl
            return hankel_transform(Qijl, qs, quadrature_method, work_weights, work_integrand, work_f)
        end
    end
    ρcore = hankel_transform(get_quantity(psp, CoreChargeDensity()), qs, quadrature_method, work_weights,
                             work_integrand, work_f)
    ρval = hankel_transform(get_quantity(psp, ValenceChargeDensity()), qs, quadrature_method, work_weights,
                            work_integrand, work_f)

    return UltrasoftPsP{TT,FourierSpace}(psp.identifier, psp.Zatom, psp.Zval, psp.lmax, Vloc, β, psp.D, χ, Q,
                                         psp.q, ρcore, ρval)
end

function interpolate_onto(maximum_spacing::T,
                          psp::UltrasoftPsP{T,S}) where {T<:Real,S<:EvaluationSpace}
    Vloc = interpolate_onto(maximum_spacing, get_quantity(psp, LocalPotential()))
    β = map(get_quantity(psp, BetaProjector())) do βl
        map(βl) do βln
            return interpolate_onto(maximum_spacing, βln)
        end
    end
    χ = map(get_quantity(psp, ChiProjector())) do χl
        map(χl) do χln
            return interpolate_onto(maximum_spacing, χln)
        end
    end
    Q = map(get_quantity(psp, AugmentationFunction())) do Ql
        map(Ql) do Qijl
            return interpolate_onto(maximum_spacing, Qijl)
        end
    end
    ρcore = interpolate_onto(maximum_spacing, get_quantity(psp, CoreChargeDensity()))
    ρval = interpolate_onto(maximum_spacing, get_quantity(psp, ValenceChargeDensity()))
    return UltrasoftPsP{T,S}(psp.identifier, psp.Zatom, psp.Zval, psp.lmax, Vloc, β, psp.D, χ, Q, psp.q, ρcore,
                             ρval)
end
