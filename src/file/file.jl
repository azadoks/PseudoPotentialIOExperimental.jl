@doc raw"""
Abstract type representing a pseudopotential file.

The structure of the data should closely mirror the format of the file, and the values
of quantities should be exactly those found in the file.
"""
abstract type PsPFile end

#!!! Required functions !!!#
"""
Identifying data (preferably unique).
"""
function identifier(file::PsPFile) end

"""
SHA1 checksum of file contents.
"""
function checksum(file::PsPFile) end

"""
Pseudopotential file format.
"""
function format(file::PsPFile) end

"""
Exchange-correlation functional in the file's format.
"""
function functional(file::PsPFile) end

"""
The element which the pseudopotential was constructed to reproduce.
"""
function element(file::PsPFile) end

"""
Maximum angular momentum channel in the local part of the pseudopotential.
"""
function max_angular_momentum(file::PsPFile) end

"""
Pseudo-atomic valence charge.
"""
function ionic_charge(file::PsPFile) end

"""
Whether the pseudopotential is of the norm-conserving kind.
"""
function is_norm_conserving(file::PsPFile) end

"""
Whether the pseudopotential is of the ultrasoft kind.
"""
function is_ultrasoft(file::PsPFile) end

"""
Whether the pseudopotential is of the plane-augmented wave kind.
"""
function is_paw(file::PsPFile) end

"""
Whether the pseudopotential contains relativistic spin-orbit coupling data.
"""
function has_spin_orbit(file::PsPFile) end

"""
Whether the pseudopotential contains non-linear core corrections.
"""
function has_nlcc(file::PsPFile) end

#!!! Convenience functions !!!#
"""
Formalism of the pseudopotential.
"""
function formalism(file::PsPFile)::String
    is_paw(file) && return "projector-augmented wave"
    is_ultrasoft(file) && return "ultrasoft"
    is_norm_conserving(file) && return "norm-conserving"
end

"""
Type of relativistic treatment (fully relativistic or scalar-relativistic).
"""
relativistic_treatment(file::PsPFile)::Symbol = has_spin_orbit(file) ? :full : :scalar

Base.Broadcast.broadcastable(file::PsPFile) = Ref(file)

Base.:(==)(file1::PsPFile, file2::PsPFile) = file1.checksum == file2.checksum
Base.hash(file::PsPFile) = hash(file.checksum)

function Base.show(io::IO, file::PsPFile)
    typename = string(typeof(file))
    el = element(file)
    z = ionic_charge(file)
    return print(io, "$typename(element=$el, z_valence=$z)")
end

function Base.show(io::IO, ::MIME"text/plain", file::PsPFile)
    println(io, typeof(file))
    @printf "%032s: %s\n" "identifier" identifier(file)
    @printf "%032s: %s\n" "SHA1 checksum" checksum(file)
    @printf "%032s: %s\n" "format" format(file)
    @printf "%032s: %s\n" "element" element(file)
    @printf "%032s: %f\n" "valence charge" ionic_charge(file)
    @printf "%032s: %d\n" "maximum angular momentum" max_angular_momentum(file)
    @printf "%032s: %s\n" "relativistic treatment" relativistic_treatment(file)
    @printf "%032s: %s\n" "exchange-correlation functional" functional(file)
    @printf "%032s: %s\n" "non-linear core correction" has_nlcc(file)
    @printf "%032s: %s\n" "formalism" formalism(file)
    @printf "%032s: %s\n" "number of beta projectors" n_projector_radials(file)
    @printf "%032s: %s" "number of chi projectors" n_orbital_radials(file)
end
