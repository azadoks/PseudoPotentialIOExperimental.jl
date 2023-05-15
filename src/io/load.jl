_FILE_EXT_LOADERS = Dict(".upf" => UpfFile,
                         ".psp8" => Psp8File,
                         ".hgh" => HghFile)

"""
Parse a pseudopotential file into a `PsPFile` struct.
"""
function load_psp_file(path::AbstractString)
    isfile(path) || throw(SystemError("$path"))
    _, ext = lowercase.(splitext(path))
    ext in keys(_FILE_EXT_LOADERS) ||
        throw(ArgumentError("$ext: Unknown pseudopotential file extension"))
    return _FILE_EXT_LOADERS[ext](path)
end
function load_psp_file(family_name_or_dir::AbstractString, filename::AbstractString)
    dir = resolve_family(family_name_or_dir)
    return load_psp_file(joinpath(dir, filename))
end

"""
Load a pseudopotential file into its corresponding `AbstractPsP` subtype.
"""
load_psp(file::PsPFile) = formalism(file)(file)
load_psp(file::HghFile) = HghPsP(file)
load_psp(path::AbstractString) = load_psp(load_psp_file(path))
function load_psp(family_name_or_dir::AbstractString, filename::AbstractString)
    dir = resolve_family(family_name_or_dir)
    return load_psp(joinpath(dir, filename))
end

"""
Load all pseudopotentials in a family into `PsPFile` structs.
"""
function load_family_psp_files(filepaths::AbstractVector{T}) where {T<:AbstractString}
    return load_psp_file.(filepaths)
end
function load_family_psp_files(family_name_or_dir::AbstractString)
    dir = resolve_family(family_name_or_dir)
    psps = list_family_psps(dir; with_info=false)
    filepaths = map(filename -> joinpath(dir, filename), psps)
    return load_family_psp_files(filepaths)
end

"""
Load and standardize all pseudopotentials into `AbstractPsP` structs.
"""
function load_family_psps(filepaths::AbstractVector{T}) where {T<:AbstractString}
    return load_psp.(filepaths)
end
function load_family_psps(family_name_or_dir::AbstractString)
    dir = resolve_family(family_name_or_dir)
    psps = list_family_psps(dir; with_info=false)
    filepaths = map(filename -> joinpath(dir, filename), psps)
    return load_family_psps(filepaths)
end
