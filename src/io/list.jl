"""
List all known pseudopotential families with additional information about the number and
format(s) of the pseudopotentials within if requested (slow!).
"""
function list_families()
    artifacts_toml_path = Artifacts.find_artifacts_toml(pathof(PseudoPotentialIO))
    artifacts = Artifacts.load_artifacts_toml(artifacts_toml_path)
    return collect(keys(artifacts))
end

"""
List all files in the directory of a given family which have a supported pseudopotential
file extension (.upf, .psp8, or .hgh).
"""
function list_family_psps(family_name_or_dir::AbstractString; with_info=false)
    dir = resolve_family(family_name_or_dir)
    (_, _, filenames) = first(walkdir(dir))
    psp_filenames = filter(f -> !startswith(f, "."), filenames)
    psp_filenames = filter(f -> lowercase(splitext(f)[2]) in keys(_FILE_EXT_LOADERS),
                           psp_filenames)
    return psp_filenames
end

# Resolve the name of a PsP family artifact into its directory
# or return the input if it's already a directory
function resolve_family(family_name_or_dir::AbstractString)
    isdir(family_name_or_dir) && return family_name_or_dir
    try
        return @artifact_str(family_name_or_dir)
    catch
        error("PsP family $family_name_or_dir does not exist")
    end
end
