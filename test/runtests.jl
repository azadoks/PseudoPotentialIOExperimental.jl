using Test
using PseudoPotentialIOExperimental
using Aqua
using LazyArtifacts
using Random
using PeriodicTable
using Bessels
using QuadGK
using JSON
using Quaternions
using LinearAlgebra
using BSplineKit
import SpecialFunctions: erfi

Random.seed!(0)

TAGS = ARGS
isempty(TAGS) && (TAGS = ["all"])

println("\tRunning tests (TAGS = $(join(TAGS, ", "))).")

include("fixtures.jl")

@testset "PseudoPotentialIOExperimental.jl" begin
    if any(in.(("all", "aqua"), Ref(TAGS)))
        include("aqua.jl")
    end

    ## Common
    # if any(in.(("all", "common", "hankel_transform"), Ref(TAGS)))
    #     include("common/hankel_transform.jl")
    # end

    if any(in.(("all", "common", "mesh"), Ref(TAGS)))
        include("common/mesh.jl")
    end

    if any(in.(("all", "common", "quadrature"), Ref(TAGS)))
        include("common/quadrature.jl")
    end

    if any(in.(("all", "common", "spherical_bessel"), Ref(TAGS)))
        include("common/spherical_bessel.jl")
    end

    ## I/O
    if any(in.(("all", "io", "load"), Ref(TAGS)))
        include("io/load.jl")
    end

    if any(in.(("all", "io", "list"), Ref(TAGS)))
        include("io/list.jl")
    end

    if any(in.(("all", "io", "show"), Ref(TAGS)))
        include("io/show.jl")
    end

    ## File formats
    if any(in.(("all", "file"), Ref(TAGS)))
        include("file/file.jl")
    end

    if any(in.(("all", "file", "upf"), Ref(TAGS)))
        include("file/upf.jl")
    end

    if any(in.(("all", "file", "upf", "upf1"), Ref(TAGS)))
        include("file/upf1.jl")
    end

    if any(in.(("all", "file", "upf", "upf2"), Ref(TAGS)))
        include("file/upf2.jl")
    end

    if any(in.(("all", "file", "psp8"), Ref(TAGS)))
        include("file/psp8.jl")
    end

    if any(in.(("all", "file", "hgh"), Ref(TAGS)))
        include("file/hgh.jl")
    end

    ## Pseudopotentials
    if any(in.(("all", "psp"), Ref(TAGS)))
        include("psp/psp.jl")
    end

    if any(in.(("all", "psp", "numeric"), Ref(TAGS)))
        include("psp/numeric.jl")
    end

    if any(in.(("all", "psp", "numeric", "norm_conserving"), Ref(TAGS)))
        include("psp/norm_conserving.jl")
    end

    if any(in.(("all", "psp", "numeric", "ultrasoft"), Ref(TAGS)))
        include("psp/ultrasoft.jl")
    end

    if any(in.(("all", "psp", "hgh"), Ref(TAGS)))
        include("psp/hgh.jl")
    end

    if any(in.(("all", "psp", "numeric", "norm_conserving", "hgh", "upf_hgh"), Ref(TAGS)))
        include("psp/upf_hgh.jl")
    end

    if any(in.(("all", "psp", "numeric", "norm_conserving", "upf_psp8"), Ref(TAGS)))
        include("psp/upf_psp8.jl")
    end

    ## Deprecated
    if any(in.(("all", "deprecated"), Ref(TAGS)))
        include("deprecated/upf.jl")
        include("deprecated/upf_json.jl")
        include("deprecated/upf_psp8.jl")
    end
end
