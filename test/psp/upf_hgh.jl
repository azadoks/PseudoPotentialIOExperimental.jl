@testset "UPF-HGH agreement" begin
    qs = [0.01, 0.5, 2.5, 5.0, 10.0, 22.0]

    @testset "$(splitpath(hgh_path)[end])" for (upf2_path, hgh_path) in UPF2_HGH_FILEPATHS
        upf2_r = load_psp(upf2_path)
        upf2_q = hankel_transform(upf2_r, qs)
        hgh_r = load_psp(hgh_path)
        hgh_q = hankel_transform(hgh_r)

        @test atomic_charge(upf2_r) == atomic_charge(hgh_r)
        @test valence_charge(upf2_r) == valence_charge(hgh_r)
        @test max_angular_momentum(upf2_r) == max_angular_momentum(hgh_r)

        @testset "LocalPotential()" begin
            upf2_Vloc_r = get_quantity(upf2_r, LocalPotential())
            hgh_Vloc_r = get_quantity(hgh_r, LocalPotential())
            rgrid = collect(upf2_Vloc_r.r)
            @test all(isapprox.(upf2_Vloc_r.(rgrid), hgh_Vloc_r.(rgrid), rtol=1e-7, atol=1e-7))

            upf2_Vloc_q = get_quantity(upf2_q, LocalPotential())
            hgh_Vloc_q = get_quantity(hgh_q, LocalPotential())
            @test all(isapprox.(upf2_Vloc_q.(qs), hgh_Vloc_q.(qs), rtol=1e-6, atol=1e-6))
        end

        @testset "BetaProjector()" begin
            @testset "l=$(l)" for l in angular_momenta(hgh_r)
                @test n_radials(upf2_r, BetaProjector(), l) == n_radials(hgh_r, BetaProjector(), l)
                @testset "n=$(n)" for n in 1:n_radials(hgh_r, BetaProjector(), l)
                    upf2_β_r = get_quantity(upf2_r, BetaProjector(), l, n)
                    hgh_β_r = get_quantity(hgh_r, BetaProjector(), l, n)
                    rgrid = collect(upf2_β_r.r)
                    # HGH provides β in units of Ha / √a₀, while UPF has been converted to 1 / √a₀
                    @test all(isapprox.(upf2_β_r.(rgrid), rgrid .^ 2 .* hgh_β_r.(rgrid) .* 2, rtol=1e-7,
                                        atol=1e-7))

                    upf2_β_q = get_quantity(upf2_q, BetaProjector(), l, n)
                    hgh_β_q = get_quantity(hgh_q, BetaProjector(), l, n)
                    # HGH provides β in units of Ha / √a₀, while UPF has been converted to 1 / √a₀
                    @test all(isapprox.(upf2_β_q.(qs), hgh_β_q.(qs) .* 2, rtol=1e-6, atol=1e-6))
                end
            end
        end

        @testset "BetaCoupling()" begin
            for l in angular_momenta(hgh_r)
                # HGH provides D in units of 1 / Ha, while UPF has been converted to Ha
                @test all(isapprox.(get_quantity(upf2_r, BetaCoupling(), l),
                                    get_quantity(hgh_r, BetaCoupling(), l) ./ 4; rtol=1e-8, atol=1e-8))
            end
        end

        @testset "energy_correction" begin
            @test energy_correction(upf2_r) .≈ energy_correction(hgh_r) rtol = 1e-5 atol = 1e-5
        end
    end
end
