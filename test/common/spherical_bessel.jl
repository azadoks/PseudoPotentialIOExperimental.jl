import PseudoPotentialIOExperimental: fast_sphericalbesselj

@testset "Fast spherical Bessel" begin
    for l in 0:10
        for x in (rand(100) .* 100)
            @test sphericalbesselj(l, x) ≈ fast_sphericalbesselj(l)(x) atol = 1e-8
        end
    end
end
