@testset "Quadrature" begin
    function sin_info(mesh::RadialMesh)
        y = sin.(mesh)
        I_true = cos(first(mesh)) - cos(last(mesh))
        return (; y, I_true)
    end

    function exp_sin_info(mesh::RadialMesh)
        xn = last(mesh)
        y = exp.(-(mesh .^ 2)) .* sin.(mesh)
        I_true = real(-(√π * (-2erfi(1 / 2) + erfi(1 / 2 - im * xn) + erfi(1 / 2 + im * xn))) /
                      (4 * exp(1 / 4)))
        return (; y, I_true)
    end

    quadrature_methods = (Trapezoid(), Simpson(), QESimpson(), AbinitCorrectedTrapezoid())
    @testset "$(quadrature_method)" for quadrature_method in quadrature_methods
        @testset "$(test_case)" for test_case in (sin_info, exp_sin_info)
            @testset "∫(0,$(xn))[f(x) dx]" for xn in (Float64(2π), Float64(5π / 2))
                @testset "even intervals $(isodd(n))" for n in (10001, 10002)
                    meshes = ((UniformMesh(0.0, xn, n), 1e-7),
                              (LogMeshWithoutZero(1e-6, xn, n), 1e-5),
                              (ArbitraryMesh(collect(range(0.0, xn, n))), 1e-7))
                    @testset "$(typeof(mesh))" for (mesh, atol) in meshes
                        if (quadrature_method == QESimpson()) && !isodd(n)
                            atol = 1e-1  # the QE method is _very_ bad in this case
                        end
                        # the ABINIT method doesn't support non-uniform meshes
                        if (quadrature_method == AbinitCorrectedTrapezoid()) && !(typeof(mesh) <: UniformMesh)
                            continue
                        end
                        info = test_case(mesh)
                        weights = integration_weights(mesh, quadrature_method)
                        @test isapprox(dot(weights, info.y), info.I_true; atol)
                    end
                end
            end
        end
    end
end
