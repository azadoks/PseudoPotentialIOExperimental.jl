import PseudoPotentialIOExperimental: resolve_family

@testset "List families" begin
    @test eltype(list_families()) <: AbstractString
end

@testset "List pseudopotentials" begin
    for family in TEST_FAMILIES
        @test eltype(list_family_psps(family)) <: AbstractString
    end
end

@testset "Reslove family directory" begin
    # Test resolution of artifact family directories from the artifact name
    for family in TEST_FAMILIES
        family_dir = resolve_family(family)
        @test isdir(family_dir)
    end

    # Test the `isdir` condition of resolve_family
    (root, dirs, _) = first(walkdir("./"))
    for dir in dirs
        family_dir = resolve_family(joinpath(root, dir))
        @test joinpath(root, dir) == family_dir
    end

    # Test the error condition of resolve_family where the argument corresponds to
    # neither a known artifact nor an existing directory
    @test_throws ErrorException resolve_family("./abcdefg/")
end
