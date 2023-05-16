module PseudoPotentialIOExperimental
using Artifacts
using BSplineKit
using EzXML
using LazyArtifacts
using LinearAlgebra
using OffsetArrays
using Polynomials
using Printf
using Statistics
using SHA
using PrettyTables
using OrderedCollections

using PeriodicTable: PeriodicTable
import Base.Broadcast.broadcastable
import Bessels: gamma, sphericalbesselj
import SpecialFunctions: erf

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

include("quantity/mesh.jl")
include("common/quadrature.jl")
include("common/fast_sphericalbesselj.jl")
include("common/hankel_transform.jl")

include("quantity/evaluation_space.jl")
export CoulombCorrection
export ErfCoulombCorrection
include("quantity/quantity.jl")
export maximum_radius
export minimum_radius
export extrema_radius
export interpolate_onto
export hankel_transform
include("quantity/numeric.jl")
include("quantity/analytical/analytical.jl")
include("quantity/analytical/hgh.jl")
export PsPQuantityFlag
export LocalPotential
export ValenceChargeDensity
export CoreChargeDensity
export ChiProjector
export BetaProjector
include("quantity/flag.jl")

## File datastructures and interface
export PsPFile
export format
export element
export is_norm_conserving
export is_ultrasoft
export is_paw
export formalism
export has_spin_orbit
export relativistic_treatment
export has_core_density
export valence_charge
export max_angular_momentum
export n_radials
include("file/file.jl")
export UpfFile
include("file/upf.jl")
include("file/upf1.jl")
include("file/upf2.jl")
export Psp8File
include("file/psp8.jl")
export HghFile
include("file/hgh.jl")

export AbstractPsP
export identifier
export element
export max_angular_momentum
export angular_momenta
export valence_charge
export atomic_charge
export has_spin_orbit
export get_quantity
include("psp/psp.jl")
export NumericPsP
include("psp/numeric/numeric.jl")
export NormConservingPsP
include("psp/numeric/norm_conserving.jl")
export UltrasoftPsP
include("psp/numeric/ultrasoft.jl")
export ProjectorAugmentedWavePsP
include("psp/numeric/paw.jl")
export AnalyticPsP
include("psp/analytic/analytic.jl")
export HghPsP
include("psp/analytic/hgh.jl")

## Loading/listing/showing functions
export list_families
export list_family_psps
include("io/list.jl")
export load_family_psp_files
export load_family_psps
export load_psp_file
export load_psp
include("io/load.jl")
export show_family_summary
export show_family_table
export show_family_periodic_table
include("io/show.jl")

end
