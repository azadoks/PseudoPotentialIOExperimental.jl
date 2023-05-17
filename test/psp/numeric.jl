#TODO The reference values are computed by using adaptive Gaussian quadrature to integrate
#TODO the order-4 (cubic) natural B-spline interpolators for each quantity.
#TODO
#TODO In theory, this should be equivalent to performing Simpson's method integration
#TODO (the default quadrature method) on the raw data.
#TODO
#TODO However, the results do _not_ agree super well. Interpolating the data onto a finer
#TODO mesh before performing the Hankel-Fourier transform using Simpson's method quadrature
#TODO improves agreement, but the improvement saturates quickly.
#TODO
#TODO Possible sources of error:
#TODO  - Mesh implementation
#TODO  - Fast spherical Bessel `j` implementation
#TODO  - Simpson's method quadrature implementation

#* Note: q-values are equivalent to the largest q at an Ecut via q^2 = 2Ecut
#*       Ecut up to ≈250 Ha is considered.

@testset "Numeric" begin
    qs = [0.01, 0.5, 2.5, 5.0, 10.0, 22.0]

    interp_sets = [
        (Δr=nothing, density_tol=1e-5, proj_tol=1e-1, vloc_tol=1e-3),
        (Δr=0.01, density_tol=1e-5, proj_tol=1e-3, vloc_tol=1e-3),
        (Δr=0.001, density_tol=1e-6, proj_tol=1e-4, vloc_tol=1e-4)
    ]

    @testset "$(name)" for (name, filepath) in NUMERIC_CASE_FILEPATHS
        psp = load_psp(filepath)
        @testset "Interp. Δr=$(Δr)" for (Δr, density_tol, proj_tol, vloc_tol) in interp_sets
            if !isnothing(Δr)
                psp_r = interpolate_onto(psp, Δr)
            else
                psp_r = psp
            end
            psp_q = hankel_transform(psp_r, qs)

            @testset "$(flag)" for flag in (CoreChargeDensity(), ValenceChargeDensity())
                if has_quantity(psp_r, flag)
                    quantity_r = get_quantity(psp_r, flag)
                    rmin = minimum_radius(psp_r, flag)
                    rmax = maximum_radius(psp_r, flag)

                    psp_quantity_q = get_quantity(psp_q, flag)
                    quantity_q = hankel_transform(quantity_r, qs)

                    integrand(r, q) = 4π * quantity_r(r) * sphericalbesselj(0, q * r)
                    ref = map(qs) do q
                        quadgk(Base.Fix2(integrand, q), rmin, rmax)[1]
                    end

                    @testset "q=$q" for (i, q) in enumerate(qs)
                        @test psp_quantity_q(q) ≈ ref[i] rtol=density_tol atol=density_tol
                        @test psp_quantity_q.f[i] ≈ ref[i] rtol=density_tol atol=density_tol
                        @test quantity_q(q) ≈ ref[i] rtol=density_tol atol=density_tol
                        @test quantity_q.f[i] ≈ ref[i] rtol=density_tol atol=density_tol
                    end
                end
            end

            @testset "$(flag)" for flag in (BetaProjector(), ChiProjector())
                if has_quantity(psp_r, flag)
                    @testset "l=$(l)" for l in angular_momenta(psp_r)
                        @testset "n=$(n)" for n in 1:n_radials(psp_r, flag, l)
                            quantity_r = get_quantity(psp_r, flag, l, n)

                            rmin = minimum_radius(psp_r, flag, l, n)
                            rmax = maximum_radius(psp_r, flag, l, n)

                            psp_quantity_q = get_quantity(psp_q, flag, l, n)
                            quantity_q = hankel_transform(quantity_r, qs)

                            integrand(r, q) = 4π * quantity_r(r) * sphericalbesselj(l, q * r)
                            ref = map(qs) do q
                                quadgk(Base.Fix2(integrand, q), rmin, rmax)[1]
                            end

                            @testset "q=$q" for (i, q) in enumerate(qs)
                                @test psp_quantity_q(q) ≈ ref[i] rtol=proj_tol atol=proj_tol
                                @test psp_quantity_q.f[i] ≈ ref[i] rtol=proj_tol atol=proj_tol
                                @test quantity_q(q) ≈ ref[i] rtol=proj_tol atol=proj_tol
                                @test quantity_q.f[i] ≈ ref[i] rtol=proj_tol atol=proj_tol
                            end
                        end
                    end
                end
            end

            @testset "LocalPotential()" begin
                flag = LocalPotential()
                quantity_r = get_quantity(psp_r, flag)
                rmin = minimum_radius(psp_r, flag)
                rmax = maximum_radius(psp_r, flag)

                psp_quantity_q = get_quantity(psp_q, flag)
                quantity_q = hankel_transform(quantity_r, qs)

                correction_r = ErfCoulombCorrection(valence_charge(psp_r))
                correction_q = hankel_transform(correction_r)
                function integrand(r, q)
                    return 4π * (r * quantity_r(r) - correction_r(r)) * sin(q * r)
                end
                ref = map(qs) do q
                    quadgk(Base.Fix2(integrand, q), rmin, rmax)[1] / q + 4π * correction_q(q)
                end

                @testset "q=$q" for (i, q) in enumerate(qs)
                    @test psp_quantity_q(q) ≈ ref[i] rtol=vloc_tol atol=vloc_tol
                    @test psp_quantity_q.f[i] ≈ ref[i] rtol=vloc_tol atol=vloc_tol
                    @test quantity_q(q) ≈ ref[i] rtol=vloc_tol atol=vloc_tol
                    @test quantity_q.f[i] ≈ ref[i] rtol=vloc_tol atol=vloc_tol
                end
            end

            # TODO: This is a bit fucked
            # @testset "energy_correction" begin
            #     q_small = 1e-5
            #     qs = [q_small, 2q_small, 3q_small, 4q_small, 5q_small]
            #     Zval = valence_charge(psp_r)
            #     Vloc_q = hankel_transform(get_quantity(psp_r, LocalPotential()), qs)
            #     ref = Vloc_q.f[1] + 4π * Zval / q_small^2
            #     @test ref ≈ energy_correction(psp_r) atol=1e-3
            # end

        end
    end
end
