module PseudoPotentialIOExperimental
using Artifacts
using EzXML
using LazyArtifacts
using LinearAlgebra
using OffsetArrays
using Printf
using Statistics
using SHA
using PrettyTables
using OrderedCollections

using PeriodicTable: PeriodicTable
import Base.Broadcast.broadcastable

## DocStringExtensions Templates
using DocStringExtensions
@template (FUNCTIONS, METHODS, MACROS) = """
                                         $(TYPEDSIGNATURES)
                                         $(DOCSTRING)
                                         $(METHODLIST)
                                         """

@template TYPES = """
                  $(TYPEDEF)
                  $(DOCSTRING)
                  $(TYPEDFIELDS)
                  """

## File datastructures and interface
# (mirror file layout)
export PsPFile
export format
export element
export is_norm_conserving
export is_ultrasoft
export is_paw
export formalism
export has_spin_orbit
export relativistic_treatment
export valence_charge
export max_angular_momentum
include("file/file.jl")
export UpfFile
include("file/upf.jl")
include("file/upf1.jl")
include("file/upf2.jl")
export Psp8File
include("file/psp8.jl")
export HghFile
include("file/hgh.jl")

## Deprecated loaders
export load_upf
export load_psp8
include("deprecated/upf.jl")
include("deprecated/psp8.jl")

## Loading/listing/showing functions
export list_families
export list_family_psps
include("io/list.jl")
export load_family_psp_files
export load_psp_file
include("io/load.jl")
export show_family_summary
export show_family_table
export show_family_periodic_table
include("io/show.jl")
end
