"""
Show a brief summary of a pseudopotential family, listing its name, directory on disk,
the pseudopotential formalisms and file-formats it contains, as well as the number of
pseudopotentials it provides.
"""
function show_family_summary(family_name_dir_or_data; name="", dir="")
    name, dir = determine_name_dir(family_name_dir_or_data, name, dir)
    return show_family_summary(Core.stdout, family_name_dir_or_data; name, dir)
end
function show_family_summary(io::IO, family_dir_or_name::AbstractString; name="", dir="")
    name, dir = determine_name_dir(family_dir_or_name, name, dir)
    return show_family_summary(io, load_family_psp_files(family_dir_or_name); name, dir)
end
function show_family_summary(io::IO, family_psp_files::AbstractVector{T};
                             name="", dir="") where {T<:PsPFile}
    data = get_pseudos_data(family_psp_files)
    return show_family_summary(io::IO, data; name, dir)
end
function show_family_summary(io::IO, family_data::OrderedDict; name="", dir="")
    @printf io "%12s: %s\n" "Family" name
    @printf io "%12s: %s\n" "Directory" dir
    @printf io "%12s: %s\n" "Formalism(s)" join(string.(unique(family_data["Formalism"])),
                                                ", ")
    @printf io "%12s: %s\n" "Format(s)" join(string.(unique(family_data["Format"])), ", ")
    @printf io "%12s: %d\n" "Pseudo(s)" length(family_data["Filename"])
    return nothing
end

"""
List the pseudopotentials in a pseudopotential family in a pretty table.
The elements for which pseudos are shown can be restricted by passing a list of strings,
e.g. ["Ag"].
"""
function show_family_table(family_name_dir_or_data; kwargs...)
    return show_family_table(Core.stdout, family_name_dir_or_data; kwargs...)
end
function show_family_table(io::IO, family_name_or_dir::AbstractString; kwargs...)
    return show_family_table(io, load_family_psp_files(family_name_or_dir); kwargs...)
end
function show_family_table(io::IO, family_psp_files; elements=[], crop=:none,
                           vlines=:all, hlines=:elements, kwargs...)
    data = get_pseudos_data(family_psp_files)
    # Filter the family to find pseudopotentials of the given elements
    if !isempty(elements)
        indices = filter(i -> data["Element"][i] in elements,
                         eachindex(data["Element"]))
        data = OrderedDict(k => v[indices] for (k, v) in data)
    end

    # Insert horizontal lines between groups of elements
    if hlines == :elements
        hlines = [0, 1]
        for i in 2:length(data["Element"])
            if data["Element"][i - 1] != data["Element"][i]
                push!(hlines, i)
            end
        end
        push!(hlines, length(data["Element"]) + 1)
    end

    header = vec(collect((keys(data))))  # keys => vector
    # Convert the values of the ordered dictionary to a matrix with the same number
    # of columns as values in the dictionary
    data = permutedims(mapreduce(permutedims, vcat, collect(values(data))), (2, 1))

    pretty_table(io, data; header, crop, vlines, hlines, kwargs...)
    return nothing
end

"""
Show the elements contained in a pseudopotential family as a periodic table.
"""
function show_family_periodic_table(family_name_dir_or_data)
    return show_family_periodic_table(Core.stdout, family_name_dir_or_data)
end
function show_family_periodic_table(io::IO, family_name_or_dir::AbstractString)
    return show_family_periodic_table(io, load_family_psp_files(family_name_or_dir))
end
function show_family_periodic_table(io::IO,
                                    family_psp_files::AbstractVector{T}) where {T<:PsPFile}
    data = get_pseudos_data(family_psp_files)
    elements = unique([PeriodicTable.elements[Symbol(element)]
                       for element in data["Element"]])
    show(io, "text/plain", PeriodicTable.Elements(elements))
    return nothing
end

function determine_name_dir(family_name_dir_or_data::AbstractString, name, dir)
    if isempty(name)
        name = isdir(family_name_dir_or_data) ? "" : family_name_dir_or_data
    end
    if isempty(dir)
        dir = isdir(family_name_dir_or_data) ? family_name_dir_or_data : ""
    end
    return name, dir
end

function get_pseudos_data(psps)
    return OrderedDict{String,Any}("Element" => map(el -> el.symbol,
                                                    PseudoPotentialIO.element.(psps)),
                                   "Filename" => identifier.(psps),
                                   "Valence Charge" => valence_charge.(psps),
                                   "NLCC" => has_core_density.(psps),
                                   "Spin Orbit" => has_spin_orbit.(psps),
                                   "Format" => format.(psps),
                                   "Formalism" => formalism.(psps),
                                   "File" => psps)
end
