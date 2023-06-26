Pkg.activate("play");
begin
    using PseudoPotentialIOExperimental
    import PseudoPotentialIOExperimental: n_angulars
    import DFTK: DFTK
    using LinearAlgebra
    using Unitful
end;

begin  # Density every atom
    function compute_structure_factors!(ρ::AbstractArray{Complex{T},3}, basis::DFTK.PlaneWaveBasis{T},
                                        position::AbstractVector) where {T}
        qs_frac = DFTK.G_vectors(basis)
        map!(Base.Fix1(dot, -position), ρ, qs_frac)  # Compute -(G ⋅ r) for all G
        ρ .= DFTK.cis2pi.(ρ)
        return ρ
    end

    function atomic_density!(ρ::AbstractArray{Complex{T}}, basis::DFTK.PlaneWaveBasis{T}, atom,
                             position::AbstractVector{T},
                             quantity::PseudoPotentialIOExperimental.DensityFlag;
                             coefficient::T=one(T)) where {T}
        qnorms_cart = norm.(DFTK.G_vectors_cart(basis))

        atomic_density = get_quantity(atom, quantity)

        compute_structure_factors!(ρ, basis, position)  # Fill ρ with structure factors
        ρ .*= atomic_density.itp.(qnorms_cart)  # Form factors
        ρ .*= coefficient  # Often magnetic moments
        ρ ./= sqrt(basis.model.unit_cell_volume)  # Normalization
        return ρ
    end

    function atomic_density(basis::DFTK.PlaneWaveBasis{T}, atom, position::AbstractVector{T},
                            quantity::PseudoPotentialIOExperimental.DensityFlag;
                            coefficient::T=one(T)) where {T}
        ρ = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
        return atomic_density!(ρ, basis, atom, position, quantity; coefficient)
    end

    function atomic_density_superposition(basis::DFTK.PlaneWaveBasis{T}, atoms, positions,
                                          quantity::PseudoPotentialIOExperimental.DensityFlag;
                                          coefficients::AbstractVector{T}=ones(T, length(atoms))) where {T}
        ρ = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})
        ρ_atom = DFTK.zeros_like(DFTK.G_vectors(basis), Complex{T})

        for (atom, position, coefficient) in zip(atoms, positions, coefficients)
            atomic_density!(ρ_atom, basis, atom, position, quantity; coefficient)
            ρ .+= ρ_atom
        end

        DFTK.enforce_real!(basis, ρ)
        return DFTK.irfft(basis, ρ)
    end
end;

begin  # Projection every atom
    function build_projection_vector!(P::AbstractVector{Complex{T}}, basis::DFTK.PlaneWaveBasis{T},
                                      kpt::DFTK.Kpoint, atom, position::AbstractVector{T},
                                      quantity::PseudoPotentialIOExperimental.ProjectorFlag,
                                      l::Integer, m::Integer, n::Integer) where {T}
        qs_frac = DFTK.Gplusk_vectors(basis, kpt)
        qs_cart = DFTK.Gplusk_vectors_cart(basis, kpt)

        atomic_projector = get_quantity(atom, quantity, l, n)

        map!(Base.Fix1(dot, -position), P, qs_frac)  # Compute -(G ⋅ r) for all G
        map!(DFTK.cis2pi, P, P)  # Structure factors
        P .*= (-im)^l .* DFTK.ylm_real.(l, m, qs_cart) .* atomic_projector.itp.(norm.(qs_cart))  # Form factors
        P ./= sqrt(basis.model.unit_cell_volume)  # Normalization
        return P
    end

    # function build_projection_vector!(P::AbstractVector{Complex{T}}, basis::DFTK.PlaneWaveBasis{T},
    #                                 kpt::DFTK.Kpoint, atom, position::AbstractVector{T},
    #                                 quantity::PseudoPotentialIOExperimental.ProjectorFlag,
    #                                 l::Integer, n::Integer, m::Integer) where {T}
    #     qs_cart = DFTK.Gplusk_vectors_cart(basis, kpt)
    #     atomic_projector = get_quantity(atom, quantity, l, n)
    #     radial = atomic_projector.itp.(norm.(qs_cart))
    #     return build_projection_vector!(P, basis, kpt, atom, position, radial, l, m)
    # end

    function build_projection_vector(basis::DFTK.PlaneWaveBasis{T}, kpt::DFTK.Kpoint, args...) where {T}
        P = DFTK.zeros_like(Gplusk_vectors(basis, kpt), Complex{T})
        return build_projection_vector!(P, basis, kpt, args...)
    end

    # function build_projection_vectors!(P::AbstractArray{Complex{T}}, basis::DFTK.PlaneWaveBasis{T},
    #                                    kpt::DFTK.Kpoint, atom, position,
    #                                    quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
    #     i_proj = 1
    #     for l in angular_momenta(atom), m in (-l):(+l), n in 1:n_radials(atom, quantity, l)
    #         build_projection_vector!(@view(P[:,i_proj]), basis, kpt, atom, position, quantity, l, m, n)
    #         i_proj += 1
    #     end
    # end

    function build_projection_vectors!(P::AbstractArray{Complex{T}},
                                                 radial_work::AbstractArray{Complex{T}},
                                                 angular_work::AbstractVector{Complex{T}},
                                                 basis::DFTK.PlaneWaveBasis{T},
                                                 kpt::DFTK.Kpoint, atom, position,
                                                 quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        qs_frac = DFTK.Gplusk_vectors(basis, kpt)
        qs_cart = DFTK.Gplusk_vectors_cart(basis, kpt)
        normalization = sqrt(basis.model.unit_cell_volume)

        # Precompute structure factors so they aren't re-computed for each l, m, n
        map!(Base.Fix1(dot, -position), @view(P[:,1]), qs_frac)  # Compute -(G ⋅ r) for all G
        map!(DFTK.cis2pi, @view(P[:,1]), @view(P[:,1]))  # Structure factors
        @views for i in axes(P, 2)
            P[:,i] .= P[:,1]
        end

        i_proj = 1
        for l in angular_momenta(atom)
            # Precompute radial parts so they aren't re-computed for each m
            # In testing, this gives ~2x speedup with large Ecut and high k-density
            @views for n in 1:n_radials(atom, quantity, l)
                atomic_projector = get_quantity(atom, quantity, l, n)
                radial_work[:, n] .= atomic_projector.itp.(norm.(qs_cart))
            end
            @views for m in (-l):(+l)
                # Precompute angular parts so they aren't re-computed for each n
                angular_work[:] .= (-im)^l .* DFTK.ylm_real.(l, m, qs_cart)
                for n in 1:n_radials(atom, quantity, l)
                    P[:,i_proj] .*= angular_work  # Angular part of the form factors
                    P[:,i_proj] .*= radial_work[:, n]  # Radial part of the form factors
                    P[:,i_proj] ./= normalization  # Normalization
                    i_proj += 1
                end
            end
        end
        return P
    end

    function build_projection_vectors(basis::DFTK.PlaneWaveBasis{T}, kpt::DFTK.Kpoint, atom, position::DFTK.Vec3{T},
                                      quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        n_q = length(DFTK.Gplusk_vectors(basis, kpt))
        n_proj = n_angulars(atom, quantity)
        P = DFTK.zeros_like(DFTK.Gplusk_vectors(basis, kpt), Complex{T}, n_q, n_proj)
        return build_projection_vectors!(P, basis, kpt, atom, position, quantity)
    end

    function build_projection_vectors!(P::AbstractArray{Complex{T}}, basis::DFTK.PlaneWaveBasis{T},
                                                 kpt::DFTK.Kpoint, atoms::AbstractVector, positions::AbstractVector,
                                                 quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        # Allocate working arrays for the angular and radial parts of the form factors
        max_n_radials = maximum(atoms) do atom
            maximum(angular_momenta(atom)) do l
                n_radials(atom, quantity, l)
            end
        end
        radial_work = DFTK.zeros_like(P, size(P, 1), max_n_radials)
        angular_work = DFTK.zeros_like(P, size(P, 1))

        i_proj = 1
        for (atom, position) in zip(atoms, positions)
            n_atom_projs = n_angulars(atom, quantity)
            build_projection_vectors!(@view(P[:,i_proj:(i_proj + n_atom_projs - 1)]),
                                      radial_work, angular_work, basis, kpt, atom, position,
                                      quantity)
            i_proj += n_atom_projs
        end
        return P
    end

    function build_projection_vectors(basis::DFTK.PlaneWaveBasis{T}, kpt::DFTK.Kpoint, atoms::AbstractVector,
                                      positions::AbstractVector,
                                      quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        n_q = length(DFTK.Gplusk_vectors(basis, kpt))
        n_proj = sum(Base.Fix2(n_angulars, quantity), atoms)
        P = DFTK.zeros_like(DFTK.Gplusk_vectors(basis, kpt), Complex{T}, n_q, n_proj)
        return build_projection_vectors!(P, basis, kpt, atoms, positions, quantity)
    end
end;

begin  # Projection unique species
    function compute_structure_factors!(structure_factor_work::AbstractVector{Complex{T}}, basis::DFTK.PlaneWaveBasis{T},
                                        kpt::DFTK.Kpoint, position::AbstractVector) where {T}
        qs_frac = DFTK.Gplusk_vectors(basis, kpt)
        structure_factor_work .= -dot.(qs_frac, Ref(position))
        structure_factor_work .= DFTK.cis2pi.(structure_factor_work)
        return structure_factor_work
    end

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

begin  # Projection unique species conservative
    function build_projection_vectors_conservative(basis::DFTK.PlaneWaveBasis{T}, kpt::DFTK.Kpoint, species::AbstractVector,
                                      species_positions::Union{AbstractVector{Vector{Vector{Float64}}}, AbstractVector{Vector{DFTK.Vec3{Float64}}}},
                                      quantity::PseudoPotentialIOExperimental.ProjectorFlag) where {T}
        qs_frac = DFTK.Gplusk_vectors(basis, kpt)
        qs_cart = DFTK.Gplusk_vectors_cart(basis, kpt)
        qnorms_cart = norm.(qs_cart)

        n_q = length(qs_frac)
        n_proj = sum(zip(species, species_positions)) do (specie, positions)
            n_angulars(specie, quantity) * length(positions)
        end
        P = DFTK.zeros_like(qs_frac, Complex{T}, n_q, n_proj)

        i_proj = 1
        for (specie, positions) in zip(species, species_positions)
            for position in positions
                for l in angular_momenta(specie)
                    for m in (-l):(+l)
                        for n in 1:n_radials(specie, quantity, l)
                            atomic_projector = get_quantity(specie, quantity, l, n)
                            P[:,i_proj] .= DFTK.cis2pi.(-dot.(qs_frac, Ref(position)))  # Copy over the structure factors
                            P[:,i_proj] .*= (-im)^l .* DFTK.ylm_real.(l, m, qs_cart)  # Angular part of the form factors
                            P[:,i_proj] .*= atomic_projector.(qnorms_cart)  # Radial part of the form factors
                            P[:,i_proj] ./= sqrt(basis.model.unit_cell_volume)  # Normalization
                            i_proj += 1
                        end
                    end
                end
            end
        end
        return P
    end
end;

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

    model = DFTK.Model(lattice, atoms_DFTK, positions; terms=[DFTK.Kinetic()])
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

begin  # Big system
    lattice = [
        [ 8.173000000,  0.2832058330,  1.9828117054];;
        [ 0.000000000, 12.8658833920,  0.9056813287];;
        [ 0.000000000,  0.0000000000, 12.7086008360]
    ] .* u"Å"

    family_path = PseudoPotentialIOExperimental.resolve_family("pd_nc_sr_pbe_standard_0.4.1_upf")
    elements_DFTK = Dict(
        :Al => DFTK.ElementPsp(:Al; psp=DFTK.load_psp(family_path, "Al.upf")),
        :Ca => DFTK.ElementPsp(:Ca; psp=DFTK.load_psp(family_path, "Ca.upf")),
        :O => DFTK.ElementPsp(:O; psp=DFTK.load_psp(family_path, "O.upf")),
        :Si => DFTK.ElementPsp(:Si; psp=DFTK.load_psp(family_path, "Si.upf"))
    )

    atoms_DFTK = [
        repeat([elements_DFTK[:Ca]], 8)...,
        repeat([elements_DFTK[:Al]], 16)...,
        repeat([elements_DFTK[:Si]], 16)...,
        repeat([elements_DFTK[:O]], 64)...
    ]
    positions = [
        [0.821600000000, 0.986400000000, 0.913300000000],
        [0.178400000000, 0.013600000000, 0.086700000000],
        [0.274300000000, 0.031200000000, 0.456500000000],
        [0.725700000000, 0.968800000000, 0.543500000000],
        [0.767500000000, 0.535900000000, 0.458800000000],
        [0.232500000000, 0.464100000000, 0.541200000000],
        [0.311300000000, 0.505200000000, 0.925300000000],
        [0.688700000000, 0.494800000000, 0.074700000000],
        [0.604600000000, 0.161000000000, 0.388800000000],
        [0.395400000000, 0.839000000000, 0.611200000000],
        [0.614400000000, 0.665800000000, 0.887200000000],
        [0.385600000000, 0.334200000000, 0.112800000000],
        [0.126400000000, 0.815200000000, 0.882400000000],
        [0.873600000000, 0.184800000000, 0.117600000000],
        [0.113900000000, 0.314500000000, 0.378800000000],
        [0.886100000000, 0.685500000000, 0.621200000000],
        [0.467400000000, 0.113000000000, 0.848100000000],
        [0.532600000000, 0.887000000000, 0.151900000000],
        [0.476700000000, 0.611000000000, 0.332600000000],
        [0.523300000000, 0.389000000000, 0.667400000000],
        [0.991600000000, 0.871900000000, 0.327500000000],
        [0.008400000000, 0.128100000000, 0.672500000000],
        [0.996400000000, 0.377500000000, 0.818400000000],
        [0.003600000000, 0.622500000000, 0.181600000000],
        [0.095200000000, 0.159200000000, 0.895600000000],
        [0.904800000000, 0.840800000000, 0.104400000000],
        [0.098000000000, 0.656000000000, 0.395800000000],
        [0.902000000000, 0.344000000000, 0.604200000000],
        [0.607400000000, 0.815400000000, 0.386500000000],
        [0.392600000000, 0.184600000000, 0.613500000000],
        [0.605800000000, 0.320400000000, 0.890100000000],
        [0.394200000000, 0.679600000000, 0.109900000000],
        [0.983200000000, 0.103400000000, 0.335400000000],
        [0.016800000000, 0.896600000000, 0.664600000000],
        [0.978200000000, 0.606700000000, 0.850500000000],
        [0.021800000000, 0.393300000000, 0.149500000000],
        [0.513400000000, 0.882900000000, 0.812400000000],
        [0.486600000000, 0.117100000000, 0.187600000000],
        [0.497200000000, 0.378900000000, 0.326600000000],
        [0.502800000000, 0.621100000000, 0.673400000000],
        [0.969100000000, 0.124200000000, 0.004000000000],
        [0.030900000000, 0.875800000000, 0.996000000000],
        [0.502300000000, 0.125700000000, 0.516500000000],
        [0.497700000000, 0.874300000000, 0.483500000000],
        [0.999300000000, 0.624100000000, 0.513200000000],
        [0.000700000000, 0.375900000000, 0.486800000000],
        [0.479700000000, 0.624700000000, 0.003400000000],
        [0.520300000000, 0.375300000000, 0.996600000000],
        [0.569000000000, 0.991300000000, 0.856600000000],
        [0.431000000000, 0.008700000000, 0.143400000000],
        [0.065900000000, 0.989700000000, 0.362100000000],
        [0.934100000000, 0.010300000000, 0.637900000000],
        [0.562300000000, 0.487500000000, 0.364600000000],
        [0.437700000000, 0.512500000000, 0.635400000000],
        [0.065200000000, 0.493200000000, 0.861400000000],
        [0.934800000000, 0.506800000000, 0.138600000000],
        [0.265200000000, 0.101800000000, 0.919400000000],
        [0.734800000000, 0.898200000000, 0.080600000000],
        [0.793300000000, 0.096800000000, 0.394300000000],
        [0.206700000000, 0.903200000000, 0.605700000000],
        [0.272200000000, 0.595700000000, 0.395300000000],
        [0.727800000000, 0.404300000000, 0.604700000000],
        [0.794400000000, 0.603400000000, 0.920200000000],
        [0.205600000000, 0.396600000000, 0.079800000000],
        [0.326800000000, 0.855400000000, 0.855700000000],
        [0.673200000000, 0.144600000000, 0.144300000000],
        [0.792100000000, 0.851800000000, 0.396600000000],
        [0.207900000000, 0.148200000000, 0.603400000000],
        [0.312800000000, 0.355900000000, 0.388500000000],
        [0.687200000000, 0.644100000000, 0.611500000000],
        [0.791400000000, 0.358700000000, 0.866700000000],
        [0.208600000000, 0.641300000000, 0.133300000000],
        [0.121000000000, 0.279600000000, 0.864900000000],
        [0.879000000000, 0.720400000000, 0.135100000000],
        [0.626900000000, 0.290900000000, 0.352600000000],
        [0.373100000000, 0.709100000000, 0.647400000000],
        [0.125000000000, 0.776900000000, 0.365600000000],
        [0.875000000000, 0.223100000000, 0.634400000000],
        [0.641800000000, 0.796500000000, 0.849000000000],
        [0.358200000000, 0.203500000000, 0.151000000000],
        [0.103600000000, 0.680600000000, 0.895600000000],
        [0.896400000000, 0.319400000000, 0.104400000000],
        [0.592400000000, 0.689900000000, 0.398700000000],
        [0.407600000000, 0.310100000000, 0.601300000000],
        [0.093600000000, 0.178800000000, 0.389900000000],
        [0.906400000000, 0.821200000000, 0.610100000000],
        [0.590400000000, 0.196300000000, 0.902500000000],
        [0.409600000000, 0.803700000000, 0.097500000000],
        [0.009100000000, 0.105900000000, 0.808300000000],
        [0.990900000000, 0.894100000000, 0.191700000000],
        [0.469200000000, 0.102500000000, 0.315300000000],
        [0.530800000000, 0.897500000000, 0.684700000000],
        [0.980100000000, 0.607900000000, 0.321000000000],
        [0.019900000000, 0.392100000000, 0.679000000000],
        [0.511100000000, 0.604300000000, 0.798100000000],
        [0.488900000000, 0.395700000000, 0.201900000000],
        [0.006900000000, 0.874000000000, 0.789300000000],
        [0.993100000000, 0.126000000000, 0.210700000000],
        [0.548800000000, 0.856400000000, 0.280300000000],
        [0.451200000000, 0.143600000000, 0.719700000000],
        [0.044800000000, 0.362800000000, 0.266800000000],
        [0.955200000000, 0.637200000000, 0.733200000000],
        [0.496400000000, 0.369700000000, 0.803000000000],
        [0.503600000000, 0.630300000000, 0.197000000000],
    ]

    model = DFTK.Model(lattice, atoms_DFTK, positions; terms=[DFTK.Kinetic()])
    # Ecut is 128 Ha / 256 Ry -- as big as you would ever go (for psp testing, basically)
    # kgrid is ~2x QE input generator's "very fine 0.15 1/Å"
    # basis = DFTK.PlaneWaveBasis(model; Ecut=128, kgrid=[28, 28, 28])
    basis = DFTK.PlaneWaveBasis(model; Ecut=42, kgrid=[6, 4, 4])

    psp_groups = [group for group in model.atom_groups if model.atoms[first(group)] isa DFTK.ElementPsp]
    psps = [model.atoms[first(group)].psp for group in psp_groups]
    psp_positions = [model.positions[group] for group in psp_groups]

    ########################################################################################
    qmax_Gpk = maximum(basis.kpoints) do kpoint
        return maximum(norm.(DFTK.Gplusk_vectors_cart(basis, kpoint)))
    end
    qmax_G = maximum(norm.(DFTK.G_vectors_cart(basis)))
    qmax = max(qmax_Gpk, qmax_G)

    elements_PPIO_r = Dict(
        :Al => PseudoPotentialIOExperimental.load_psp(family_path, "Al.upf"),
        :Ca => PseudoPotentialIOExperimental.load_psp(family_path, "Ca.upf"),
        :O => PseudoPotentialIOExperimental.load_psp(family_path, "O.upf"),
        :Si => PseudoPotentialIOExperimental.load_psp(family_path, "Si.upf")
    )
    elements_PPIO_r = Dict(
        key => interpolate_onto(value, 0.01) for (key, value) in elements_PPIO_r
    )
    elements_PPIO_q = Dict(
        key => hankel_transform(value, 0.0:0.01:(qmax + 0.1)) for (key, value) in elements_PPIO_r
    )

    atoms_PPIO = [
        repeat([elements_PPIO_q[:Ca]], 8)...,
        repeat([elements_PPIO_q[:Al]], 16)...,
        repeat([elements_PPIO_q[:Si]], 16)...,
        repeat([elements_PPIO_q[:O]], 64)...
    ]
    species_PPIO = [elements_PPIO_q[:Al], elements_PPIO_q[:Ca], elements_PPIO_q[:O], elements_PPIO_q[:Si]]
    species_positions_PPIO = psp_positions
end;

let  # Density (valence)
    @time ρ_ref = DFTK.guess_density(basis, DFTK.ValenceDensityPseudo())
    N_ref = sum(ρ_ref) * model.unit_cell_volume / prod(basis.fft_size)

    @time ρ = atomic_density_superposition(basis, atoms_PPIO, positions, PseudoValenceDensity())
    N = sum(ρ) * basis.model.unit_cell_volume / prod(basis.fft_size)
    if !isnothing(basis.model.n_electrons) && (N > 0)
        ρ .*= basis.model.n_electrons / N  # Renormalize to the correct number of electrons
    end
    N = sum(ρ) * model.unit_cell_volume / prod(basis.fft_size)

    println("$(N_ref) $(N)")
    @assert isapprox(sum(ρ), sum(ρ_ref); rtol=1e-4)
end

begin  # Fast projection test all atoms
    ikpt = 1
    kpt = basis.kpoints[ikpt]

    psp_groups = [group for group in model.atom_groups if model.atoms[first(group)] isa DFTK.ElementPsp]
    psps = [model.atoms[first(group)].psp for group in psp_groups]
    psp_positions = [model.positions[group] for group in psp_groups]

    P_ppio = build_projection_vectors(basis, kpt, atoms_PPIO, positions, NonLocalProjector())
    P_ref = DFTK.build_projection_vectors_(basis, kpt, psps, psp_positions)

    extrema(real(P_ppio ./ 2 .- P_ref)), extrema(imag(P_ppio ./ 2 .- P_ref))
end

let  # Projection all atoms
    psp_groups = [group for group in model.atom_groups if model.atoms[first(group)] isa DFTK.ElementPsp]
    psps = [model.atoms[first(group)].psp for group in psp_groups]
    psp_positions = [model.positions[group] for group in psp_groups]

    for (ikpt, kpt) in enumerate(basis.kpoints)
        P_ppio = build_projection_vectors(basis, kpt, atoms_PPIO, positions, NonLocalProjector())
        P_ref = DFTK.build_projection_vectors_(basis, kpt, psps, psp_positions)

        for atom in atoms_PPIO
            i_proj = 1
            for l in angular_momenta(atom), m in (-l):(+l), n in 1:n_radials(atom, NonLocalProjector(), l)
                @assert all(isapprox.(P_ppio[:, i_proj] ./ 2, P_ref[:, i_proj]; rtol=1e-6, atol=1e-6)) "$(ikpt) $(element(atom)) l=$l m=$m n=$n i_proj=$i_proj"
                i_proj += 1
            end
        end
    end
end

let  # Projection unique species
    @time let
        qmax_Gpk = maximum(basis.kpoints) do kpoint
            return maximum(norm.(DFTK.Gplusk_vectors_cart(basis, kpoint)))
        end
        qmax_G = maximum(norm.(DFTK.G_vectors_cart(basis)))
        qmax = max(qmax_Gpk, qmax_G)

        Si_r = PseudoPotentialIOExperimental.load_psp(family_path, "Si.upf")
        Si_r = interpolate_onto(Si_r, UniformMesh(0.0:0.01:5.99))
        Si_q = hankel_transform(Si_r, 0.0:0.01:(qmax + 0.1))
        Ge_r = PseudoPotentialIOExperimental.load_psp(family_path, "Ge.upf")
        Ge_r = interpolate_onto(Ge_r, UniformMesh(0.0:0.01:5.99))
        Ge_q = hankel_transform(Ge_r, 0.0:0.01:(qmax + 0.1))
        species_PPIO = [Si_q, Ge_q]

        for kpt in basis.kpoints
            build_projection_vectors(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector())
        end
    end
    @time for kpt in basis.kpoints
        DFTK.build_projection_vectors_(basis, kpt, psps, psp_positions)
    end
end;

begin  # Fast projection test unique species
    ikpt = 1
    kpt = basis.kpoints[ikpt]

    P_ppio = build_projection_vectors(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector()) ./ 2
    # P_ppio = build_projection_vectors_conservative(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector()) ./ 2
    P_ref = DFTK.build_projection_vectors_(basis, kpt, psps, psp_positions)

    extrema(real(P_ppio .- P_ref)), extrema(imag(P_ppio .- P_ref))
end

let  # Fast projection timing unique species
    ikpt = 1
    kpt = basis.kpoints[ikpt]

    @time build_projection_vectors(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector()) ./ 2
    # @time build_projection_vectors_conservative(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector()) ./ 2
    @time DFTK.build_projection_vectors_(basis, kpt, psps, psp_positions)
end;

let  # Slow projection timing unique species
    @time for kpt in basis.kpoints
        build_projection_vectors(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector()) ./ 2
    end
    # @time build_projection_vectors_conservative(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector()) ./ 2
    @time for kpt in basis.kpoints
        DFTK.build_projection_vectors_(basis, kpt, psps, psp_positions)
    end
end;

begin  # Plot projectors PPIO vs DFTK
    fig = Figure()

    qgrid = 0:0.1:10
    natom = 1
    for (atom_PPIO, atom_DFTK) in zip(atoms_PPIO, atoms_DFTK)
        for l in angular_momenta(atom_PPIO)
            ax = Axis(fig[l+1,natom], title="$(element(atom_PPIO)) $l")
            for n in 1:n_radials(atom_PPIO, NonLocalProjector(), l)
                β_PPIO = atom_PPIO.β[l][n].(qgrid)
                β_DFTK = [DFTK.eval_psp_projector_fourier(atom_DFTK.psp, n, l, q) for q in qgrid]

                lines!(ax, qgrid, (β_PPIO ./ 2 .- β_DFTK) ./ β_DFTK, label="$n")
            end
            axislegend(ax)
        end
        natom += 1
    end

    fig
end

let  # Projection
    for (ikpt, kpt) in enumerate(basis.kpoints)
        P_ppio = build_projection_vectors(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector())
        P_ref = DFTK.build_projection_vectors_(basis, kpt, psps, psp_positions)

        for atom in atoms_PPIO
            i_proj = 1
            for l in angular_momenta(atom), m in (-l):(+l), n in 1:n_radials(atom, NonLocalProjector(), l)
                @assert all(isapprox.(P_ppio[:, i_proj] ./ 2, P_ref[:, i_proj]; rtol=1e-6, atol=1e-6)) "$(ikpt) $(element(atom)) l=$l m=$m n=$n i_proj=$i_proj"
                i_proj += 1
            end
        end
    end
end

let
    ops = map(basis.kpoints) do kpt
        P = build_projection_vectors(basis, kpt, species_PPIO, species_positions_PPIO, NonLocalProjector()) ./ 2
        D = DFTK.build_projection_coefficients_(Float64, psps, psp_positions)
        DFTK.NonlocalOperator(basis, kpt, P, D)
    end
    push!(basis.terms, DFTK.TermAtomicNonlocal(ops))
end;
