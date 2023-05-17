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
    using Bessels
end

begin
    hgh_file = load_psp_file("hgh_lda_hgh", "si-q4.hgh")
    hgh = load_psp(hgh_file)
    hgh_q = hankel_transform(hgh)

    upf_file = load_psp_file("hgh_lda_upf", "Si.pz-hgh.UPF")
    upf = load_psp(upf_file)
    upf_q = hankel_transform(upf)
end

begin
    hgh = load_psp("hgh_pbe_hgh", "ir-q9.hgh")
    hgh_q = hankel_transform(hgh)

    l = 0
    n = 2
    qgrid = 0.0:0.1:10.0
    r_test = 200:-0.1:0

    β_r = hgh.β[l][n]
    β_q = hgh_q.β[l][n].(qgrid)
    integrand(r, q) = 4π * r^2 * β_r(r) * sphericalbesselj(l, q * r)
    rmax = r_test[findfirst(r -> !isapprox(0, β_r(r)), r_test)]
    β_q_quadgk = [quadgk(Base.Fix2(integrand, q), 0.0, rmax)[1] for q in qgrid]

    diff_q = β_q .- β_q_quadgk

    let
        fig = Figure()
        ax = Axis(fig[1,1])
        lines!(ax, qgrid, β_q)
        lines!(ax, qgrid, β_q_quadgk)
        ax = Axis(fig[2,1])
        lines!(ax, qgrid, diff_q)
        fig
    end
end


begin
    l = 0
    n = 2
    qgrid = 0.0:0.1:10.0
    rgrid = upf.β[l][n].r

    diff_r = upf.β[l][n].(rgrid) .- rgrid.^2 .* hgh.β[l][n].(rgrid)
    diff_q = upf_q.β[l][n].(qgrid) .- hgh_q.β[l][n].(qgrid)

    let
        fig = Figure()

        ax = Axis(fig[1,1])
        lines!(ax, rgrid, upf.β[l][n].(rgrid))
        lines!(ax, rgrid, rgrid.^2 .* hgh.β[l][n].(rgrid))
        ax = Axis(fig[2,1])
        lines!(ax, rgrid, diff_r)


        ax = Axis(fig[1,2])
        lines!(ax, qgrid, upf_q.β[l][n].(qgrid))
        lines!(ax, qgrid, hgh_q.β[l][n].(qgrid))
        ax = Axis(fig[2,2])
        lines!(ax, qgrid, diff_q)
        fig
    end
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
