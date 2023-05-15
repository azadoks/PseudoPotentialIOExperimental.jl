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
    # psp_file = load_psp_file("hgh_lda_upf", "Si.pz-hgh.UPF")
    psp_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf")
    r = psp_file.mesh.r
    rab = psp_file.mesh.rab
    z_valence = Float64(psp_file.header.z_valence)
end

begin
    local_ = psp_file.local_ ./ 2
    r_sub = r[eachindex(local_)]
    rab_sub = rab[eachindex(local_)]
    r_mesh = ArbitraryMesh(r_sub, rab_sub)

    Vloc_r = NumericLocalPotential{Float64, RealSpace}(r_mesh, local_, CubicSpline(r_sub, local_))
    corr_r = CoulombCorrection{Float64,RealSpace}(z_valence)
    Vloc_q = hankel_transform(Vloc_r, Trapezoid(), corr_r)

    let
        fig = Figure()
        ax = Axis(fig[1,1])
        lines!(ax, Vloc_q.r[50:end], Vloc_q.f[50:end])
        # lines!(ax, Vloc_q.r[50:end], Vloc_q.itp(Vloc_q.r[50:end]))
        fig
    end
end

begin
    i = 3

    beta = psp_file.nonlocal.betas[i].beta
    n = 1
    l = psp_file.nonlocal.betas[i].angular_momentum
    j = 0.0

    ircut = psp_file.nonlocal.betas[i].cutoff_radius_index
    rmin = r[1]
    rmax = r[ircut]
    beta = beta[1:ircut]
    r_sub = r[eachindex(beta)]
    rab_sub = rab[eachindex(beta)]
    r_mesh = ArbitraryMesh(r_sub, rab_sub)

    beta_r = NumericProjector{Float64, Vector{Float64}, RealSpace}(n, l, j, r_mesh, beta, CubicSpline(r_sub, beta))

    jₗ = fast_sphericalbesselj(l)
    beta_q_quad = map(0:0.01:30) do q
        ignd(r) = beta_r.itp(r) * jₗ(q * r)
        4π * quadgk(ignd, rmin, rmax)[1]
    end

    itp = CubicSpline(r_sub, beta)
    r_mesh = UniformMesh(rmin, rmax, 6000)
    r_sub = collect(r_mesh)
    beta = itp.(r_mesh)

    beta_r_itp = NumericProjector{Float64, Vector{Float64}, RealSpace}(n, l, j, r_mesh, beta, CubicSpline(r_sub, beta))
    beta_q = hankel_transform(beta_r_itp, Simpson())

    let
        fig = Figure()
        ax = Axis(fig[1,1])
        lines!(ax, beta_r.r, beta_r.f)
        ax = Axis(fig[2,1])
        lines!(ax, beta_q.r, beta_q.f .- beta_q_quad)
        fig
    end

end

begin
    rhoatom = psp_file.rhoatom
    r_sub = r[eachindex(rhoatom)]
    rab_sub = rab[eachindex(rhoatom)]
    r_mesh = ArbitraryMesh(r_sub, rab_sub)
    rmin = first(r_sub)
    rmax = last(r_sub)

    rhoatom_r = NumericDensity{Float64, Vector{Float64}, RealSpace}(r_mesh, rhoatom, CubicSpline(r_sub, rhoatom))

    jₗ = fast_sphericalbesselj(l)
    rhoatom_q_quad = map(0:0.01:30) do q
        function ignd(r::T)::Float64 where {T<:Real}
            rhoatom_r(r) * fast_sphericalbesselj0(q * r)
        end
        4π * quadgk(ignd, rmin, rmax)[1]
    end

    r_mesh_u = UniformMesh(rmin, rmax, 3000)
    r_sub_u = collect(r_mesh_u)
    rhoatom_u = rhoatom_r.(r_mesh_u)

    rhoatom_r_itp = NumericDensity{Float64, Vector{Float64}, RealSpace}(r_mesh_u, rhoatom_u, CubicSpline(r_sub_u, rhoatom_u))
    rhoatom_q = hankel_transform(rhoatom_r_itp, Simpson())

    let
        fig = Figure()
        ax = Axis(fig[1,1])
        lines!(ax, rhoatom_r.r, rhoatom_r.f)
        ax = Axis(fig[2,1])
        lines!(ax, rhoatom_q.r, rhoatom_q.f .- rhoatom_q_quad)
        fig
    end

end
