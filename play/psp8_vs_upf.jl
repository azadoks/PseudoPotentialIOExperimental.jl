Pkg.activate("play")
begin
    using PseudoPotentialIOExperimental
    import PseudoPotentialIOExperimental: ArbitraryMesh, UniformMesh, LogMesh, deriv
    import PseudoPotentialIOExperimental: NumericLocalPotential, NumericProjector, NumericDensity
    import PseudoPotentialIOExperimental: hankel_transform
    import PseudoPotentialIOExperimental: Simpson, Trapezoid, QESimpson, AbinitCorrectedTrapezoid, integration_weights
    import PseudoPotentialIOExperimental: CoulombCorrection, ErfCoulombCorrection
    import PseudoPotentialIOExperimental: RealSpace, FourierSpace
    import PseudoPotentialIOExperimental: fast_sphericalbesselj0, fast_sphericalbesselj
    using CairoMakie
    using CubicSplines
    using LinearAlgebra
    using QuadGK
    using NumericalIntegration
end

begin
    upf_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf")
    upf = load_psp(upf_file)
    psp8_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_psp8", "Si.psp8")
    psp8 = load_psp(psp8_file)
end

begin
    upf_quantity = upf.Vloc
    psp8_quantity = psp8.Vloc

    new_mesh = UniformMesh(0.0, 5.99, 6001)
    new_f = CubicSpline(upf_quantity.r, upf_quantity.f).(new_mesh)
    new_quantity = NumericLocalPotential(new_mesh, new_f)

    upf_quantity_q = hankel_transform(upf_quantity, CoulombCorrection{Float64, RealSpace}(upf.Zval))
    psp8_quantity_q = hankel_transform(psp8_quantity, CoulombCorrection{Float64, RealSpace}(upf.Zval))
    new_quantity_q = hankel_transform(new_quantity, CoulombCorrection{Float64, RealSpace}(upf.Zval))

    let
        fig = Figure()
        ax = Axis(fig[1,1])
        # lines!(ax, upf_quantity_q.r, upf_quantity_q.f .- psp8_quantity_q.f)
        lines!(ax, upf_quantity_q.r, new_quantity_q.f .- psp8_quantity_q.f)
        ax = Axis(fig[2,1])
        lines!(ax, new_quantity_q.r, new_quantity_q.f)
        lines!(ax, psp8_quantity_q.r, psp8_quantity_q.f)
        ylims!(ax, -0.1, 0.1)
        fig
    end

end

# TODO: end goal features
begin
    upf_real = load_psp("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf")

    maximum_spacing = 0.001  # Upper bound on the resolution of the quantities in real-space
    @time upf_real_splined = cubic_spline(maximum_spacing, upf_real)
end;

begin
    upf_real = load_psp("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf")

    maximum_spacing = 0.001  # Upper bound on the resolution of the quantities in real-space
    upf_real_splined = cubic_spline(maximum_spacing, upf_real)

    qgrid = 0:0.01:30
    @time upf_splined_fourier = hankel_transform(
        upf_real_splined;
        qs=qgrid,
        quadrature_method=Simpson(),
        local_potential_correction=CoulombCorrection(upf_real)
    )
end;

begin
    upf_real = load_psp("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf")
    maximum_spacing = 0.001  # Upper bound on the resolution of the quantities in real-space
    qgrid = range(0, 10, 1001)
    fine_qgrid = range(0, 10, 10_000_000)
    Vloc_fine_qgrid = similar(fine_qgrid)
    @time let
        upf_real_splined = cubic_spline(maximum_spacing, upf_real)

        upf_splined_fourier = hankel_transform(
            upf_real_splined;
            qs=qgrid,
            quadrature_method=Simpson(),
            local_potential_correction=CoulombCorrection(upf_real)
        )

        get_quantity(upf_splined_fourier, LocalPotential())(Vloc_fine_qgrid, fine_qgrid)
    end;
end;
