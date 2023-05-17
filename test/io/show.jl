_, io = mktemp()

@testset "Show family summary" begin
    for family in TEST_FAMILIES
        @test isnothing(show_family_summary(io, family))
        @test isnothing(show_family_summary(io, load_family_psp_files(family)))
    end
end

@testset "Show family table" begin
    for family in TEST_FAMILIES
        @test isnothing(show_family_table(io, family))
        @test isnothing(show_family_table(io, load_family_psp_files(family)))
        @test isnothing(show_family_table(io, load_family_psp_files(family), elements=["Si"]))
    end
end

@testset "Show family periodic table" begin
    for family in TEST_FAMILIES
        @test isnothing(show_family_periodic_table(io, family))
        @test isnothing(show_family_periodic_table(io, load_family_psp_files(family)))
    end
end

close(io)
