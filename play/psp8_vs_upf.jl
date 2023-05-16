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
    using LinearAlgebra
    using QuadGK
    using NumericalIntegration
    using BSplineKit
end

let
    upf_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf")
    upf = load_psp(upf_file)
    psp8_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_psp8", "Si.psp8")
    psp8 = load_psp(psp8_file)

    upf_quantity = upf.Vloc
    psp8_quantity = psp8.Vloc

    new_mesh = UniformMesh(0.0, 5.99, 6001)

    @time b_spline_f = BSplineKit.interpolate(upf_quantity.r, upf_quantity.f, BSplineOrder(4), Natural()).(new_mesh)
    b_spline_quantity = NumericLocalPotential(new_mesh, b_spline_f)

    # upf_quantity_q = hankel_transform(upf_quantity, CoulombCorrection{Float64, RealSpace}(upf.Zval))
    # psp8_quantity_q = hankel_transform(psp8_quantity, CoulombCorrection{Float64, RealSpace}(upf.Zval))
    # new_quantity_q = hankel_transform(new_quantity, CoulombCorrection{Float64, RealSpace}(upf.Zval))

    let
        fig = Figure()
        ax = Axis(fig[1,1])

        scatter!(upf_quantity.r, upf_quantity.f)
        lines!(new_mesh, b_spline_f)

        ylims!(ax, -4.77, -4.76)
        xlims!(ax, -0.001, 0.05)
        fig
    end

end

begin
    upf_real = load_psp("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf")
    maximum_spacing = 0.001  # Upper bound on the resolution of the quantities in real-space
    qgrid = range(0, 10, 1001)
    fine_qgrid = range(0, 10, 100_000_000)
    Vloc_fine_qgrid = similar(fine_qgrid)
    @time let
        upf_real_splined = interpolate_onto(maximum_spacing, upf_real)

        upf_splined_fourier = hankel_transform(
            upf_real_splined;
            qs=qgrid,
            quadrature_method=Simpson(),
            local_potential_correction=CoulombCorrection(upf_real)
        )

        quantity = get_quantity(upf_splined_fourier, LocalPotential())
        Vloc_fine_qgrid .= quantity.(fine_qgrid)
    end;
end;
