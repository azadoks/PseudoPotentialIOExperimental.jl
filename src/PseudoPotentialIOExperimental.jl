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

## Flagging structs
# (used in checking for and retrieving quantities from PsPFile and AbstractPsP structs)
export PsPQuantityFlag
export LocalPotential
export ValenceChargeDensity
export CoreChargeDensity
export ProjectorFlag
export ChiProjector
export BetaProjector
export BetaCoupling
export AugmentationFunction
export AugmentationCoupling
include("quantity/flag.jl")
export RealSpace
export FourierSpace
include("quantity/evaluation_space.jl")

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

## Pseudopotential quantities
# (contain physical quantities and describe how to evaluate, interpolate,
#  and Fourier transform them)
include("quantity/quantity.jl")
export RadialMesh
export ArbitraryMesh
export UniformMesh
export LogMeshWithZero
export LogMeshWithoutZero
export deriv
include("quantity/mesh.jl")
export QuadratureMethod
export Trapezoid
export Simpson
export QESimpson
export AbinitCorrectedTrapezoid
export integration_weights
export integration_weights!
include("common/quadrature.jl")
export CoulombCorrection
export ErfCoulombCorrection
include("quantity/local_potential_correction.jl")
export maximum_radius
export minimum_radius
export extrema_radii
export interpolate_onto
export energy_correction
include("quantity/numeric.jl")
include("quantity/analytical/analytical.jl")
include("quantity/analytical/hgh.jl")

## Pseudopotential structs
# (hold summary information and organize physical quantities)
export AbstractPsP
export identifier
export element
export max_angular_momentum
export angular_momenta
export valence_charge
export atomic_charge
export has_spin_orbit
export get_quantity
export has_quantity
export n_radials
export n_angulars
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

## Common utilities
export fast_sphericalbesselj
export fast_sphericalbesselj0
include("common/fast_sphericalbesselj.jl")
export hankel_transform
include("common/hankel_transform.jl")

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
export load_family_psps
export load_psp_file
export load_psp
include("io/load.jl")
export show_family_summary
export show_family_table
export show_family_periodic_table
include("io/show.jl")
end
