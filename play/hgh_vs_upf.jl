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

begin
    hgh_file = load_psp_file("hgh_lda_hgh", "si-q4.hgh")
    hgh = load_psp(hgh_file)
    upf_file = load_psp_file("hgh_lda_upf", "Si.pz-hgh.UPF")
    upf = load_psp(upf_file)
end

begin
    qgrid = 0.0:0.01:10

    upf_q = hankel_transform(upf; local_potential_correction=ErfCoulombCorrection(upf))
    hgh_q = hankel_transform(hgh)

    fig = Figure()
    ax = Axis(fig[1,1])

    # lines!(qgrid, upf_q.Vloc.f .- hgh_q.Vloc.(qgrid))
    # lines!(qgrid, upf_q.Vloc.(qgrid))
    # lines!(qgrid, hgh_q.Vloc.(qgrid))
    # fig

    upf_q.Vloc.(qgrid) .- hgh_q.Vloc.(qgrid)
end
