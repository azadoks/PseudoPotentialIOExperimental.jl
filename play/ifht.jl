using Pkg
Pkg.activate("/Users/azadoks/Source/PseudoPotentialIOExperimental.jl/play")
begin
    using PseudoPotentialIOExperimental
    using LinearAlgebra
    using Makie, CairoMakie
end

begin
    psp_r = psp = load_psp("pd_nc_sr_pbe_standard_0.4.1_psp8", "Si.psp8")
    psp_q = hankel_transform(psp_r, 0.0:0.01:30.0; quadrature_method=Simpson())

    r = psp_r.ρval.r
    r²f = psp_r.ρval.f

    q = psp_q.ρval.r
    F = psp_q.ρval.f

    r′ = r
end

begin
    jₗ = fast_sphericalbesselj(0)
    weights = integration_weights(RadialMesh(q), Simpson())
    q²F = q.^2 .* F

    f′ = map(r′) do ri
        4π * dot(weights, jₗ.(q .* ri) .* q²F) / (2π)^3
    end
    r²f′ = r′.^2 .* f′
end

let
    fig = Figure()
    ax = Axis(fig[1,1])
    lines!(ax, r, r²f)
    lines!(ax, r′, r²f′)
    lines!(Axis(fig[2,1]), r′, psp_r.ρval.itp.(r′) .- r²f′)
    fig
end
