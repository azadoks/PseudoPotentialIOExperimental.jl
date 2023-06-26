Pkg.activate("play");
begin
    using PseudoPotentialIOExperimental
    import PseudoPotentialIOExperimental: n_angulars
    import DFTK: DFTK
    using LinearAlgebra
    using Unitful
end;

begin  # Structure factors
    function compute_structure_factors!(A::AbstractArray{Complex{T},3}, basis::DFTK.PlaneWaveBasis{T},
                                        position::AbstractVector) where {T}
        qs_frac = DFTK.G_vectors(basis)
        map!(Base.Fix1(dot, -position), A, qs_frac)  # Compute -(G ⋅ r) for all G
        A .= DFTK.cis2pi.(A)
        return A
    end

    function compute_structure_factors!(A::AbstractVector{Complex{T}}, basis::DFTK.PlaneWaveBasis{T},
                                        kpt::DFTK.Kpoint, position::AbstractVector) where {T}
        qs_frac = DFTK.Gplusk_vectors(basis, kpt)
        A .= -dot.(qs_frac, Ref(position))  # Compute -(G ⋅ r) for all G
        A .= DFTK.cis2pi.(A)
        return A
    end
end;

begin  # Density unique species
    function atomic_density_superposition(basis::DFTK.PlaneWaveBasis{T}, species::AbstractVector,
                                          species_positions::Union{AbstractVector{Vector{Vector{Float64}}}, AbstractVector{Vector{DFTK.Vec3{Float64}}}},
                                          quantity::PseudoPotentialIOExperimental.DensityFlag;
                                          species_coefficients::AbstractVector{Vector{T}}=[ones(T, length(positions)) for positions in species_positions]) where {T}
        ρ = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
        ρ_atom = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
        structure_factor_work = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})

        qnorms_cart = norm.(G_vectors_cart(basis))
        for (specie, positions, coefficients) in zip(species, species_positions, species_coefficients)
            atomic_density = get_quantity(specie, quantity)
            ρ_atom .= atomic_density.itp.(qnorms_cart)
            for (position, coefficient) in zip(positions, coefficients)
                compute_structure_factors!(structure_factor_work, basis, position)
                ρ .+= coefficient .* structure_factor_work .* ρ_atom
            end
        end
        ρ ./= sqrt(basis.model.unit_cell_volume)

        DFTK.enforce_real!(basis, ρ)
        return DFTK.irfft(basis, ρ)
    end
end;

begin  # Local potential unique species
    function local_potential_superposition(basis::DFTK.PlaneWaveBasis{T}, species::AbstractVector,
                                           species_positions::Union{AbstractVector{Vector{Vector{Float64}}}, AbstractVector{Vector{DFTK.Vec3{Float64}}}},) where {T}
        V = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
        V_atom = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
        structure_factor_work = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})

        qnorms_cart = norm.(G_vectors_cart(basis))
        for (specie, positions) in zip(species, species_positions)
            local_potential = get_quantity(specie, PseudoPotentialIOExperimental.LocalPotential())
            V_atom .= local_potential.itp.(qnorms_cart)
            for position in positions
                compute_structure_factors!(structure_factor_work, basis, position)
                V .+= structure_factor_work .* V_atom
            end
        end
        V ./= sqrt(basis.model.unit_cell_volume)

        DFTK.enforce_real!(basis, V)
        return DFTK.irfft(basis, V)
    end
end;

# begin  # Density every atom
#     function atomic_density!(ρ::AbstractArray{Complex{T}}, basis::DFTK.PlaneWaveBasis{T}, atom,
#                              position::AbstractVector{T},
#                              quantity::PseudoPotentialIOExperimental.DensityFlag;
#                              coefficient::T=one(T)) where {T}
#         qnorms_cart = norm.(DFTK.G_vectors_cart(basis))

#         atomic_density = get_quantity(atom, quantity)

#         compute_structure_factors!(ρ, basis, position)  # Fill ρ with structure factors
#         ρ .*= atomic_density.itp.(qnorms_cart)  # Form factors
#         ρ .*= coefficient  # Often magnetic moments
#         ρ ./= sqrt(basis.model.unit_cell_volume)  # Normalization
#         return ρ
#     end

#     function atomic_density(basis::DFTK.PlaneWaveBasis{T}, atom, position::AbstractVector{T},
#                             quantity::PseudoPotentialIOExperimental.DensityFlag;
#                             coefficient::T=one(T)) where {T}
#         ρ = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
#         return atomic_density!(ρ, basis, atom, position, quantity; coefficient)
#     end

#     function atomic_density_superposition(basis::DFTK.PlaneWaveBasis{T}, atoms, positions,
#                                           quantity::PseudoPotentialIOExperimental.DensityFlag;
#                                           coefficients::AbstractVector{T}=ones(T, length(atoms))) where {T}
#         ρ = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
#         ρ_atom = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})

#         for (atom, position, coefficient) in zip(atoms, positions, coefficients)
#             atomic_density!(ρ_atom, basis, atom, position, quantity; coefficient)
#             ρ .+= ρ_atom
#         end

#         DFTK.enforce_real!(basis, ρ)
#         return DFTK.irfft(basis, ρ)
#     end
# end;

# begin  # Local potential every atom
#     function local_potential!(V::AbstractArray{Complex{T}}, basis::DFTK.PlaneWaveBasis{T}, atom,
#                               position::AbstractVector{T}) where {T}
#         qnorms_cart = norm.(DFTK.G_vectors_cart(basis))

#         local_potential = get_quantity(atom, PseudoPotentialIOExperimental.LocalPotential())

#         compute_structure_factors!(V, basis, position)  # Fill V with structure factors
#         V .*= local_potential.itp.(qnorms_cart)  # Form factors
#         V ./= sqrt(basis.model.unit_cell_volume)  # Normalization
#         return V
#     end

#     function local_potential(basis::DFTK.PlaneWaveBasis{T}, atom,
#                              position::AbstractVector{T}) where {T}
#         V = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
#         return local_potential!(V, basis, atom, position)
#     end

#     function local_potential_superposition(basis::DFTK.PlaneWaveBasis{T}, atoms, positions) where {T}
#         V = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
#         V_atom = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})

#         for (atom, position) in zip(atoms, positions)
#             local_potential!(ρ_atom, basis, atom, position)
#             V .+= V_atom
#         end

#         DFTK.enforce_real!(basis, V)
#         return DFTK.irfft(basis, V)
#     end
# end;

begin  # Projection unique species
    function build_projection_vectors!(P::AbstractArray{Complex{T}},
                                       structure_factor_work::AbstractVector{Complex{T}},
                                       radial_work::AbstractArray{Complex{T}},
                                       angular_work::AbstractVector{Complex{T}},
                                       basis::DFTK.PlaneWaveBasis{T},
                                       kpt::DFTK.Kpoint, atom, position,
                                       quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        qs_cart = DFTK.Gplusk_vectors_cart(basis, kpt)

        # Precompute structure factors so they aren't re-computed for each l, m, n
        compute_structure_factors!(structure_factor_work, basis, kpt, position)

        i_proj = 1
        for l in angular_momenta(atom)
            @views for m in (-l):(+l)
                # Precompute angular parts so they aren't re-computed for each n
                angular_work .= (-im)^l .* DFTK.ylm_real.(l, m, qs_cart)
                for n in 1:n_radials(atom, quantity, l)
                    P[:,i_proj] .= structure_factor_work  # Copy over the structure factors
                    P[:,i_proj] .*= angular_work  # Angular part of the form factors
                    P[:,i_proj] .*= radial_work[:, n, l+1]  # Radial part of the form factors
                    P[:,i_proj] ./= sqrt(basis.model.unit_cell_volume)  # Normalization
                    i_proj += 1
                end
            end
        end
        return P
    end

    function compute_form_factor_radials!(radial_work::AbstractMatrix{Complex{T}}, qs_cart::AbstractVector{T},
                                          specie,
                                          quantity::PseudoPotentialIOExperimental.ProjectorFlag,
                                          l::Integer) where {T}
        @views for n in 1:n_radials(specie, quantity, l)
            atomic_projector = get_quantity(specie, quantity, l, n)
            radial_work[:, n] .= atomic_projector.itp.(qs_cart)
        end
        return radial_work
    end

    function compute_form_factor_radials!(radial_work::AbstractArray{Complex{T},3},
                                          basis::DFTK.PlaneWaveBasis{T}, kpt::DFTK.Kpoint,
                                          species,
                                          quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        qs_cart = norm.(DFTK.Gplusk_vectors_cart(basis, kpt))
        for l in angular_momenta(species)
            compute_form_factor_radials!(@view(radial_work[:,:,l+1]), qs_cart, species, quantity, l)
        end
        return radial_work
    end

    function build_projection_vectors!(P::AbstractArray{Complex{T}}, basis::DFTK.PlaneWaveBasis{T},
                                       kpt::DFTK.Kpoint, species::AbstractVector,
                                       species_positions::AbstractVector,
                                       quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        # Allocate working arrays for the angular and radial parts of the form factors
        max_l = maximum(max_angular_momentum, species)
        max_n_radials = maximum(species) do specie
            maximum(angular_momenta(specie)) do l
                n_radials(specie, quantity, l)
            end
        end
        radial_work = DFTK.zeros_like(P, size(P, 1), max_n_radials, max_l + 1)
        angular_work = DFTK.zeros_like(P, size(P, 1))
        structure_factor_work = DFTK.zeros_like(P, size(P, 1))

        i_proj = 1
        for (specie, positions) in zip(species, species_positions)
            n_specie_projs = n_angulars(specie, quantity)
            # Precompute the radial parts of the form factors for all l, n
            compute_form_factor_radials!(radial_work, basis, kpt, specie, quantity)
            for position in positions
                build_projection_vectors!(@view(P[:,i_proj:(i_proj + n_specie_projs - 1)]),
                                          structure_factor_work, radial_work, angular_work,
                                          basis, kpt, specie, position,
                                          quantity)
                i_proj += n_specie_projs
            end
        end
        return P
    end

    function build_projection_vectors(basis::DFTK.PlaneWaveBasis{T}, kpt::DFTK.Kpoint, species::AbstractVector,
                                      species_positions::Union{AbstractVector{Vector{Vector{Float64}}}, AbstractVector{Vector{DFTK.Vec3{Float64}}}},
                                      quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        n_q = length(DFTK.Gplusk_vectors(basis, kpt))
        n_proj = sum(zip(species, species_positions)) do (specie, positions)
            n_angulars(specie, quantity) * length(positions)
        end
        P = DFTK.zeros_like(DFTK.Gplusk_vectors(basis, kpt), Complex{T}, n_q, n_proj)
        return build_projection_vectors!(P, basis, kpt, species, species_positions, quantity)
    end
end;

function construct_nonlocal(basis, species_PPIO, species_positions_PPIO, psps, psp_positions)
    ops = map(basis.kpoints) do kpt
        P = build_projection_vectors(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector()) ./ 2
        D = DFTK.build_projection_coefficients_(Float64, psps, psp_positions)
        DFTK.NonlocalOperator(basis, kpt, P, D)
    end
    return DFTK.TermAtomicNonlocal(ops)
end

function construct_local(basis, species_PPIO, species_positions_PPIO)
    V = local_potential_superposition(basis, species_PPIO, species_positions_PPIO)
    return DFTK.TermAtomicLocal(V)
end

function construct_xc(basis::PlaneWaveBasis{T}, xc, species_PPIO, species_positions_PPIO) where {T}
    ρcore = atomic_density_superposition(basis, species_PPIO, species_positions_PPIO, PseudoPotentialIOExperimental.CoreDensity())
    ρcore = DFTK.ρ_from_total(basis, ρcore)
    functionals = map(xc.functionals) do fun
        # Strip duals from functional parameters if needed
        newparams = convert_dual.(T, parameters(fun))
        change_parameters(fun, newparams; keep_identifier=true)
    end
    return TermXc(convert(Vector{Functional}, functionals),
                  DFTK.convert_dual(T, xc.scaling_factor),
                  T(xc.potential_threshold), ρcore)
end

begin  # System
    a = 10.26  # Silicon lattice constant in Bohr
    lattice = a / 2 * [[0 1 1.0];
                       [1 0 1.0];
                       [1 1 0.0]]

    family_path = PseudoPotentialIOExperimental.resolve_family("pd_nc_sr_pbe_standard_0.4.1_upf")
    Si_DFTK = DFTK.ElementPsp(:Si; psp=DFTK.load_psp(family_path, "Si.upf"))
    Ge_DFTK = DFTK.ElementPsp(:Ge; psp=DFTK.load_psp(family_path, "Ge.upf"))

    atoms_DFTK = [Si_DFTK, Ge_DFTK]
    positions = [ones(3) / 8, -ones(3) / 8]

    model = DFTK.Model(lattice, atoms_DFTK, positions; terms=[DFTK.Kinetic(), Ewald(), PspCorrection()])
    # Ecut is 128 Ha / 256 Ry -- as big as you would ever go (for psp testing, basically)
    # kgrid is ~2x QE input generator's "very fine 0.15 1/Å"
    # basis = DFTK.PlaneWaveBasis(model; Ecut=128, kgrid=[28, 28, 28])
    basis = DFTK.PlaneWaveBasis(model; Ecut=24, kgrid=[12, 12, 12])

    psp_groups = [group for group in model.atom_groups if model.atoms[first(group)] isa DFTK.ElementPsp]
    psps = [model.atoms[first(group)].psp for group in psp_groups]
    psp_positions = [model.positions[group] for group in psp_groups]

    ########################################################################################
    qmax_Gpk = maximum(basis.kpoints) do kpoint
        return maximum(norm.(DFTK.Gplusk_vectors_cart(basis, kpoint)))
    end
    qmax_G = maximum(norm.(DFTK.G_vectors_cart(basis)))
    qmax = max(qmax_Gpk, qmax_G)

    Si_r = PseudoPotentialIOExperimental.load_psp(family_path, "Si.upf")
    Si_r = interpolate_onto(Si_r, 0.01)
    Si_q = hankel_transform(Si_r, 0.0:0.01:(qmax + 0.1))
    Ge_r = PseudoPotentialIOExperimental.load_psp(family_path, "Ge.upf")
    Ge_r = interpolate_onto(Ge_r, 0.01)
    Ge_q = hankel_transform(Ge_r, 0.0:0.01:(qmax + 0.1))
    atoms_PPIO = [Si_q, Ge_q]
    species_PPIO = atoms_PPIO
    species_positions_PPIO = [[positions[i]] for i in eachindex(positions)]
end;
