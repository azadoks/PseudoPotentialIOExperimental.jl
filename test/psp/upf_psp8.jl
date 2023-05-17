@testset "UPF-PSP8 agreement" begin
    qs = [0.01, 0.5, 2.5, 5.0, 10.0, 22.0]
    r_min = 0.0  # True for all the PseudoDojo UPF-PSP8 pairs
    Δr = 0.01  # True for all the PseudoDojo UPF-PSP8 pairs

    #* You _must_ evaluate the UPF and PSP8 files up to the lower radial cutoff to get
    #* results that agree. Often, it seems that this trunction is _before_
    #* the quantity actually decays to zero!

    #* For some elements, e.g. Li, the quantities seem to be non-exact matches between the
    #* formats; the tolerances are artifically inflated for these cases.

    @testset "$(splitpath(psp8_path)[end])" for (upf2_path, psp8_path) in UPF2_PSP8_FILEPATHS
        psp8_r = load_psp(psp8_path)
        upf2_r = load_psp(upf2_path)

        @test atomic_charge(upf2_r) == atomic_charge(psp8_r)
        @test valence_charge(upf2_r) == valence_charge(psp8_r)
        @test max_angular_momentum(upf2_r) == max_angular_momentum(psp8_r)

        @testset "$(flag)" for flag in (LocalPotential(), CoreChargeDensity(), ValenceChargeDensity())
            @test has_quantity(psp8_r, flag) == has_quantity(upf2_r, flag)
            if has_quantity(psp8_r, flag) && has_quantity(upf2_r, flag)
                upf2_quant = get_quantity(upf2_r, flag)
                psp8_quant = get_quantity(psp8_r, flag)

                r_max = min(maximum_radius(upf2_quant), maximum_radius(psp8_quant))
                r_grid = UniformMesh(r_min:Δr:r_max)

                upf2_quant_q = hankel_transform(interpolate_onto(upf2_quant, r_grid), qs)
                psp8_quant_q = hankel_transform(interpolate_onto(psp8_quant, r_grid), qs)

                @test all(isapprox.(upf2_quant.(r_grid), psp8_quant.(r_grid); atol=1e-6, rtol=1e-6))
                @test all(isapprox.(upf2_quant_q.f, psp8_quant_q.f; atol=1e-6, rtol=1e-6))
            end
        end

        @testset "BetaProjector()" begin
            @testset "l=$(l)" for l in angular_momenta(psp8_r)
                @test n_radials(upf2_r, BetaProjector(), l) == n_radials(psp8_r, BetaProjector(), l)
                @testset "n=$(n)" for n in 1:n_radials(psp8_r, BetaProjector(), l)
                    upf2_βln = get_quantity(upf2_r, BetaProjector(), l, n)
                    psp8_βln = get_quantity(psp8_r, BetaProjector(), l, n)

                    r_max = min(maximum_radius(upf2_βln), maximum_radius(psp8_βln))
                    r_grid = UniformMesh(r_min:Δr:r_max)

                    upf2_βln_q = hankel_transform(interpolate_onto(upf2_βln, r_grid), qs)
                    psp8_βln_q = hankel_transform(interpolate_onto(psp8_βln, r_grid), qs)

                    @test all(isapprox.(upf2_βln.(r_grid), psp8_βln.(r_grid); atol=1e-6, rtol=1e-6))
                    @test all(isapprox.(upf2_βln_q.f, psp8_βln_q.f; atol=1e-6, rtol=1e-6))
                end
            end
        end

        @testset "BetaCoupling()" begin
            for l in angular_momenta(psp8_r)
                @test all(isapprox.(get_quantity(psp8_r, BetaCoupling(), l), get_quantity(upf2_r, BetaCoupling(), l); atol=1e-6, rtol=1e-6))
            end
        end

        @testset "energy_correction" begin
            upf2_Vloc = get_quantity(upf2_r, LocalPotential())
            psp8_Vloc = get_quantity(psp8_r, LocalPotential())

            r_max = min(maximum_radius(upf2_Vloc), maximum_radius(psp8_Vloc))
            r_grid = UniformMesh(r_min:Δr:r_max)

            upf2_ene_corr = energy_correction(interpolate_onto(upf2_Vloc, r_grid))
            psp8_ene_corr = energy_correction(interpolate_onto(psp8_Vloc, r_grid))

            @test all(isapprox.(upf2_ene_corr, psp8_ene_corr; atol=1e-6, rtol=1e-6))
        end
    end
end
