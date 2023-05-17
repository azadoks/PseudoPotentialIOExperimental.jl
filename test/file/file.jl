QUANTITIES = [ValenceChargeDensity(), CoreChargeDensity(), BetaProjector(),
              ChiProjector(), BetaCoupling(), LocalPotential(),
              AugmentationFunction()]

@testset "AbstractPsPFile" begin
    for filepath in TEST_FILEPATHS
        file = load_psp_file(filepath)
        @test isa(identifier(file), AbstractString)
        @test isa(format(file), AbstractString)
        @test isa(element(file), PeriodicTable.Element)
        @test -1 <= max_angular_momentum(file) <= 5
        @test 0 <= n_radials(BetaProjector(), file)
        @test 0 <= n_radials(ChiProjector(), file)
        @test 0 <= valence_charge(file)
        @test isa(is_norm_conserving(file), Bool)
        @test isa(is_ultrasoft(file), Bool)
        @test isa(is_paw(file), Bool)
        @test formalism(file) in (NormConservingPsP, UltrasoftPsP, ProjectorAugmentedWavePsP)
        @test isa(has_spin_orbit(file), Bool)
        @test relativistic_treatment(file) in (:scalar, :full)
        for quantity in QUANTITIES
            @test isa(has_quantity(quantity, file), Bool)
        end

        #TODO: test n_angulars, broadcasting, hashing, show
    end
end
