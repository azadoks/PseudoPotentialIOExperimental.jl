@testset "HGH" begin
    r_test = 200:-0.1:0
    i_test = randperm(length(HGH_FILEPATHS))[1:30]  # Takes ≈3s/pseudo; 237 possible HGH pseudos
    qs = [0.01, 0.5, 2.5, 5.0, 10.0, 22.0]

    @testset "$(splitpath(filepath)[end])" for filepath in HGH_FILEPATHS[i_test]
        psp_r = load_psp(filepath)
        psp_q = hankel_transform(psp_r)

        @testset "NonLocalProjector()" begin
            @testset "l=$(l)" for l in angular_momenta(psp_r)
                @testset "n=$(n)" for n in 1:n_radials(psp_r, NonLocalProjector(), l)
                    β_r = get_quantity(psp_r, NonLocalProjector(), l, n)
                    β_q = hankel_transform(β_r)
                    psp_β_q = get_quantity(psp_q, NonLocalProjector(), l, n)

                    integrand(r, q) = 4π * r^2 * β_r(r) * sphericalbesselj(l, q * r)

                    rmax = r_test[findfirst(r -> !isapprox(0, β_r(r)), r_test)]
                    ref = map(qs) do q
                        quadgk(Base.Fix2(integrand, q), 0.0, rmax)[1]
                    end

                    @testset "q=$q" for (i, q) in enumerate(qs)
                        @test β_q(q) ≈ ref[i] rtol=1e-12 atol=1e-12
                        @test psp_β_q(q) ≈ ref[i] rtol=1e-12 atol=1e-12
                    end
                end
            end
        end

        @testset "LocalPotential()" begin
            Vloc_r = get_quantity(psp_r, LocalPotential())
            Vloc_q = hankel_transform(Vloc_r)
            psp_Vloc_q = get_quantity(psp_q, LocalPotential())

            correction_r = CoulombCorrection(valence_charge(psp_r))
            correction_q = hankel_transform(correction_r)
            integrand(r, q) = 4π * (r * Vloc_r(r) - correction_r(r)) * sin(q * r)

            rmax = r_test[findfirst(r -> !isapprox(0, Vloc_r(r)), r_test)]
            ref = map(qs) do q
                quadgk(Base.Fix2(integrand, q), 0.0, rmax)[1] / q + 4π * correction_q(q)
            end

            @testset "q=$q" for (i, q) in enumerate(qs)
                @test Vloc_q(q) ≈ ref[i] rtol=1e-9 atol=1e-9
                @test psp_Vloc_q(q) ≈ ref[i] rtol=1e-9 atol=1e-9
            end
        end

        @testset "energy_correction" begin
            q_small = 1e-5
            Zval = valence_charge(psp_q)
            ref = get_quantity(psp_q, LocalPotential())(q_small) + 4π * Zval / q_small^2
            @test ref ≈ energy_correction(psp_r) atol=1e-3
        end

    end
end
