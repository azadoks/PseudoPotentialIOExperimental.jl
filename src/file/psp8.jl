"""
ABINIT PSeudoPotential format 8 header block.
"""
struct Psp8Header
    "Description"
    title::String
    "Atomic charge"
    zatom::Float64
    "Pseudo-atomic (valence) charge"
    zion::Float64
    "Generation date `ddmmyy`"
    pspd::Int
    "PSP format version"
    pspcod::Int
    "ABINIT code for exchange-correlation if positive, LibXC code if negative"
    pspxc::Int
    "Maximum angular momentum channel of the Kleinman-Bylander projectors"
    lmax::Int
    "Angular momentum channel of the local potential (`>lmax` if none)"
    lloc::Int
    "Number of points on the radial mesh"
    mmax::Int
    "Unused"
    r2well::Float64
    "Radius beyond which the model core charge (if present) is zero or negligible"
    rchrg::Float64
    "Pseudopotential has a model core charge if > 0"
    fchrg::Float64
    "Unused"
    qchrg::Float64
    "Number of projectors for each angular momentum channel"
    nproj::Vector{Int}
    "Signals presence of format extensions: 0 (none), 1 (valence charge), 2 (spin-orbit), 3 (spin-orbit and valence charge)"
    extension_switch::Int
    "Number of spin-orbit projectors for each angular momentum (if present)"
    nprojso::Union{Nothing,Vector{Int}}
end

function psp8_parse_header(io::IO)
    # Line 1 contains the title / comment describing the pseudopotential
    title = readline(io)
    # Line 2 contains the atomic charge `zatom`, pseudoatomic charge `zion`, and creation
    # date `pspd` (ddmmyy)
    s = split(readline(io))
    zatom = _parse_fortran(Float64, s[1])
    zion = _parse_fortran(Float64, s[2])
    pspd = parse(Int, s[3])
    # Line 3 contains the format version `pspcod`, exchange-correlation functional `pspxc`,
    # maximum angular momentum channel `lmax`, equivalent angular momentum (pseudoindex of
    # angular momentum block) where the local potential is listed `lloc`, number of points
    # on the radial mesh `mmax`, and one deprecated quantity `r2well`
    s = split(readline(io))
    pspcod = parse(Int, s[1])
    pspxc = parse(Int, s[2])
    lmax = parse(Int, s[3])
    lloc = parse(Int, s[4])
    mmax = parse(Int, s[5])
    r2well = parse(Int, s[6])
    @assert pspcod == 8
    # Line 4 contains three integrated charge quantities: `rchrg``, `fchrg``, and `qchrg``.
    # `fchrg`` is used to signal the presence of non-linear core-correction when greater
    # than zero. `rchrg` and `qchrg` are not used.
    s = split(readline(io))
    rchrg = _parse_fortran(Float64, s[1])
    fchrg = _parse_fortran(Float64, s[2])
    qchrg = _parse_fortran(Float64, s[3])
    # Line 5 contains the number of KB non-local projectors per angular momentum channel
    s = split(readline(io))
    nproj = [parse(Int, s[l + 1]) for l in 0:lmax]
    if lloc <= lmax
        @assert nproj[lloc + 1] == 0
    end
    # Line 6 contains the "extension switch" flags which are used to signal the presence of
    # spin-orbit coupling (and theoretically any other extension to the format)
    s = split(readline(io))
    extension_switch = parse(Int, s[1])
    # Line 7 contains the number of spin-orbit coupling KB non-local projectors per angular
    # momentum channel when the extension switch has a value of 2 or 3
    if extension_switch in (2, 3)
        s = split(readline(io))
        nprojso = [0, [parse(Int, s[i]) for i in 1:lmax]...]
    else
        nprojso = nothing
    end

    return Psp8Header(title, zatom, zion, pspd, pspcod, pspxc, lmax, lloc, mmax, r2well,
                      rchrg, fchrg, qchrg, nproj, extension_switch, nprojso)
end

function psp8_parse_beta_projector_block(io, nproj, mmax)
    # The first line of each projector block contains first the angular momentum of the
    # projectors then the "KB energies" (projector coupling constants) among the projectors
    # in the block
    header_line = split(readline(io))
    l = parse(Int, header_line[1])
    nproj_l = nproj[l + 1]
    ekb = [_parse_fortran(Float64, header_line[i + 1]) for i in 1:nproj_l]

    # The block itself has two metadata columns (index and radial coordinate in Bohr) and
    # nproj[l + 1] data columns (the projectors).
    rgrid = Vector{Float64}(undef, mmax)
    projectors = [Vector{Float64}(undef, mmax) for _ in 1:(nproj_l)]
    for i in 1:mmax
        s = split(readline(io))
        rgrid[i] = _parse_fortran(Float64, s[2])
        for j in 1:nproj_l
            # Vector-of-vectors is used because each projector will be separately truncated
            # where it goes to zero at a later time
            projectors[j][i] = _parse_fortran(Float64, s[2 + j])
        end
    end

    return (; l, rgrid, projectors, ekb)
end

function psp8_parse_v_local_block(io, mmax)
    header_line = split(readline(io))
    l = parse(Int, header_line[1])

    rgrid = Vector{Float64}(undef, mmax)
    v_local = Vector{Float64}(undef, mmax)
    for i in 1:mmax
        s = split(readline(io))
        rgrid[i] = _parse_fortran(Float64, s[2])
        v_local[i] = _parse_fortran(Float64, s[3])
    end

    return (; l, rgrid, v_local)
end

function psp8_parse_main_blocks(io, mmax, nproj, lmax, lloc)
    beta_projector_blocks = []
    v_local_block = ()

    if lmax < lloc  # The local potential does not replace an angular momentum channel
        n_blocks = lmax + 1
    else  # The local potential takes the place of an angular momentum channel in the file
        n_blocks = lmax
    end

    for _ in 0:n_blocks
        # Record the position at the start of the block so we can read in the first line
        # and go back to it after reconaissance
        pos = position(io)
        # Read the block header
        header_line = split(readline(io))
        # Go back to the start of the block
        seek(io, pos)
        # Parse the block's angular momentum
        block_l = parse(Int, header_line[1])
        if block_l == lloc
            v_local_block = psp8_parse_v_local_block(io, mmax)
            if lloc <= lmax
                # If the local potential is mixed in with the projectors, add an empty block
                # to maintain proper angular momentum indexing
                push!(beta_projector_blocks, (; l=block_l, rgrid=nothing, projectors=[], ekb=[]))
            end
        else
            block = psp8_parse_beta_projector_block(io, nproj, mmax)
            push!(beta_projector_blocks, block)
        end
    end

    projectors = [block.projectors for block in beta_projector_blocks]
    ekb = [block.ekb for block in beta_projector_blocks]

    return (; v_local_block.rgrid, v_local_block.v_local, projectors, ekb)
end

function psp8_parse_spin_orbit_blocks(io, mmax, nprojso, lmax)
    # Spin-orbit coupling projector blocks have the same shape as "normal" projector blocks,
    # but there is no spin-orbit local potential.
    beta_projector_blocks = []
    for _ in 1:lmax
        block = psp8_parse_beta_projector_block(io, nprojso, mmax)
        push!(beta_projector_blocks, block)
    end

    projectors = [block.projectors for block in beta_projector_blocks]
    ekb = [block.ekb for block in beta_projector_blocks]

    return (; projectors, ekb)
end

function psp8_parse_rhov_block(io, mmax)
    rhov = Vector{Float64}(undef, mmax)
    for i in 1:mmax
        s = split(readline(io))
        # index  r  ρval  ????  ????
        rhov[i] = _parse_fortran(Float64, s[3])  # Has a 4π prefactor
    end
    return rhov
end

function psp8_parse_nlcc_block(io, mmax)
    # The non-linear core correction block contains 7 columns: index, radial grid, core
    # charge density (including a 4π factor), and the first four derivatives of the core
    # charge density w.r.t. radial coordinate
    radial_grid = Vector{Float64}(undef, mmax)
    rhoc = Vector{Float64}(undef, mmax)
    d_rhoc_dr = Vector{Float64}(undef, mmax)
    d2_rhoc_dr2 = Vector{Float64}(undef, mmax)
    d3_rhoc_dr3 = Vector{Float64}(undef, mmax)
    d4_rhoc_dr4 = Vector{Float64}(undef, mmax)
    for i in 1:mmax
        s = split(readline(io))
        radial_grid[i] = _parse_fortran(Float64, s[2])  # r
        rhoc[i] = _parse_fortran(Float64, s[3])         # 4π ρ
        d_rhoc_dr[i] = _parse_fortran(Float64, s[4])    # dρ / dr
        d2_rhoc_dr2[i] = _parse_fortran(Float64, s[5])  # d²ρ / dr²
        d3_rhoc_dr3[i] = _parse_fortran(Float64, s[6])  # d³ρ / dr³
        d4_rhoc_dr4[i] = _parse_fortran(Float64, s[7])  # d⁴ρ / dr⁴
    end
    return (; rhoc, d_rhoc_dr, d2_rhoc_dr2, d3_rhoc_dr3, d4_rhoc_dr4)
end

"""
ABINIT PSeudoPotential format 8 file contents. Information on the file format specification
and the meaning of the quantities within the file can be found on the
["psp8" page](https://docs.abinit.org/developers/psp8_info/) of the ABINIT documentation.
"""
struct Psp8File <: PsPFile
    "Identifier"
    identifier::String
    "SHA1 Checksum"
    checksum::Vector{UInt8}
    "Various pseudopotential metadata"
    header::Psp8Header
    "Uniform radial grid starting at `r = 0.0`"
    rgrid::Vector{Float64}
    "Local part of the pseudopotential"
    v_local::Vector{Float64}
    "Radial part of the Kleinman-Bylander projectors for each angular momentum"
    projectors::Vector{Vector{Vector{Float64}}}
    "Kleinman-Bylander energies for each angular momentum"
    ekb::Vector{Vector{Float64}}
    "Radial part of the spin-orbit Kleinman-Bylander projectors for each angular momentum"
    projectors_so::Union{Nothing,Vector{Vector{Vector{Float64}}}}
    "Spin-orbit Kleinman-Bylander energies for each angular momentum"
    ekb_so::Union{Nothing,Vector{Vector{Float64}}}
    "Model core charge density with 4π prefactor"
    rhoc::Union{Nothing,Vector{Float64}}
    "First derivative of the model core charge density"
    d_rhoc_dr::Union{Nothing,Vector{Float64}}
    "Second derivative of the model core charge density"
    d2_rhoc_dr2::Union{Nothing,Vector{Float64}}
    "Third derivative of the model core charge density"
    d3_rhoc_dr3::Union{Nothing,Vector{Float64}}
    "Fourth derivative of the model core charge density"
    d4_rhoc_dr4::Union{Nothing,Vector{Float64}}
    "Valence charge density with 4π prefactor"
    rhov::Union{Nothing, Vector{Float64}}
end

function Psp8File(io::IO; identifier="")
    checksum = SHA.sha1(io)
    seek(io, 0)
    # NOTE: parsing _must_ be done in order because it is done by reading the file
    # incrementally and depends on the order of the lines in the file
    header = psp8_parse_header(io)
    main_blocks = psp8_parse_main_blocks(io, header.mmax, header.nproj, header.lmax,
                                         header.lloc)
    if header.extension_switch in (2, 3)
        spin_orbit = psp8_parse_spin_orbit_blocks(io, header.mmax, header.nprojso,
                                                  header.lmax)
        # Add empty vectors for l=0 to maintain angular momentum <-> index correspondence
        spin_orbit = (projectors=[[], spin_orbit.projectors...],
                      ekb=[[], spin_orbit.ekb...])
    else
        spin_orbit = (projectors=nothing, ekb=nothing)
    end
    if header.fchrg > 0
        nlcc = psp8_parse_nlcc_block(io, header.mmax)
    else
        nlcc = (rhoc=nothing, d_rhoc_dr=nothing, d2_rhoc_dr2=nothing, d3_rhoc_dr3=nothing,
                d4_rhoc_dr4=nothing)
    end
    if header.extension_switch in (1, 3)
        rhov = psp8_parse_rhov_block(io, header.mmax)
    else
        rhov = nothing
    end
    return Psp8File(identifier, checksum, header, main_blocks.rgrid, main_blocks.v_local,
                    main_blocks.projectors, main_blocks.ekb, spin_orbit.projectors,
                    spin_orbit.ekb, nlcc.rhoc, nlcc.d_rhoc_dr, nlcc.d2_rhoc_dr2,
                    nlcc.d3_rhoc_dr3, nlcc.d4_rhoc_dr4, rhov)
end

function Psp8File(path::AbstractString)
    open(path, "r") do io
        return Psp8File(io; identifier=splitpath(path)[end])
    end
end

function _parse_fortran(::Type{T}, x::AbstractString) where {T<:Real}
    return parse(T, replace(lowercase(x), "d" => "e"))
end

identifier(psp::Psp8File)::String = psp.identifier
format(::Psp8File)::String = "PSP8"
function element(file::Psp8File)
    return PeriodicTable.elements[Int(file.header.zatom)]
end
has_spin_orbit(file::Psp8File)::Bool = file.header.extension_switch in (2, 3)
has_core_density(file::Psp8File)::Bool = file.header.fchrg > 0
is_norm_conserving(file::Psp8File)::Bool = true
is_ultrasoft(file::Psp8File)::Bool = false
is_paw(file::Psp8File)::Bool = false
valence_charge(file::Psp8File) = file.header.zion
max_angular_momentum(file::Psp8File)::Int = file.header.lmax
# n_radials(::NumericProjector, file::Psp8File, l::Int)::Int = file.header.nproj[l + 1]
# n_radials(::NumericState, ::Psp8File, l::Int)::Int = 0
# has_quantity(::AbstractPsPQuantity, ::Psp8File) = false
# has_quantity(::NumericLocalPotential, ::Psp8File) = true
# has_quantity(::NumericProjector, ::Psp8File) = true
# has_quantity(::BetaCoupling, ::Psp8File) = true
# has_quantity(::CoreChargeDensity, file::Psp8File) = !isnothing(file.rhoc)
# has_quantity(::ValenceChargeDensity, file::Psp8File) = !isnothing(file.rhov)
