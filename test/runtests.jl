using Test
using PseudoPotentialIOExperimental
using Aqua
using LazyArtifacts
using PeriodicTable
using JSON

TAGS = ARGS
isempty(TAGS) && (TAGS = ["all"])

println("\tRunning tests (TAGS = $(join(TAGS, ", "))).")

include("fixtures.jl")

@testset "PseudoPotentialIOExperimental.jl" begin
    if any(in.(("all", "aqua"), Ref(TAGS)))
        include("aqua.jl")
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

    ## Deprecated
    if any(in.(("all", "deprecated"), Ref(TAGS)))
        include("deprecated/mesh.jl")
        include("deprecated/upf.jl")
        include("deprecated/upf_json.jl")
        include("deprecated/upf_psp8.jl")
    end
end
