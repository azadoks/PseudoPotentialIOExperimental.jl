Pkg.activate("play")
begin
    using PseudoPotentialIOExperimental
    import PseudoPotentialIOExperimental: NumericPsPQuantity
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

function find_cutoff_index(x::AbstractVector{T}; atol=zero(T)) where {T}
    i = findfirst(fi -> abs(fi) > atol, @view x[end:-1:begin])
    return isnothing(i) ? lastindex(x) : lastindex(x) - i + 1
end

function cutoff(quantity::NumericPsPQuantity{T,S}; atol=zero(T)) where {T,S}
    i_cut = find_cutoff_index(quantity.f; atol)
    r_new = quantity.r[begin:i_cut]
    return interpolate_onto(quantity, r_new)
end

begin
    ele = "Si"

    psp8_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_psp8", "$(ele).psp8")
    psp8 = load_psp(psp8_file)
    upf_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_upf", "$(ele).upf")
    upf = load_psp(upf_file)

    l = 0

    let
        fig = Figure()
        ax1 = Axis(fig[1,1], title=ele)
        ax2 = Axis(fig[2,1])
        for n in 1:n_radials(upf, NonLocalProjector(), l)
            rmax = min(cutoff(psp8.β[l][n]).r[end], cutoff(upf.β[l][n]).r[end])
            rgrid = 0.0:0.01:rmax
            lines!(ax1, rgrid, psp8.β[l][n].(rgrid), label="psp8 β[$l][$n]")
            lines!(ax1, rgrid, upf.β[l][n].(rgrid), label="upf2 β[$l][$n]", linestyle=:dash)
            lines!(ax2, rgrid, psp8.β[l][n].(rgrid) .- upf.β[l][n].(rgrid), label="psp8-upf2 β[$l][$n]")
        end
        axislegend(ax1)
        axislegend(ax2)
        fig
    end
end

begin
    ele = "Li"

    psp8_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_psp8", "$(ele).psp8")
    psp8 = load_psp(psp8_file)
    upf_file = load_psp_file("pd_nc_sr_pbe_standard_0.4.1_upf", "$(ele).upf")
    upf = load_psp(upf_file)

    l = 0
    qgrid = 0:0.1:10

    let
        fig = Figure()
        ax1 = Axis(fig[1,1], title=ele)
        ax2 = Axis(fig[2,1])
        for n in 1:n_radials(upf, NonLocalProjector(), l)
            rmax = min(cutoff(psp8.β[l][n]).r[end], cutoff(upf.β[l][n]).r[end])
            rgrid = 0.0:0.01:rmax

            psp8_β = hankel_transform(interpolate_onto(psp8.β[l][n], UniformMesh(rgrid)), qgrid)
            upf_β = hankel_transform(interpolate_onto(upf.β[l][n], UniformMesh(rgrid)), qgrid)

            lines!(ax1, qgrid, psp8_β.(qgrid), label="psp8 β[$l][$n]")
            lines!(ax1, qgrid, upf_β.(qgrid), label="upf2 β[$l][$n]", linestyle=:dash)
            lines!(ax2, qgrid, psp8_β.(qgrid) .- upf_β.(qgrid), label="psp8-upf2 β[$l][$n]")
        end
        axislegend(ax1)
        axislegend(ax2)
        fig
    end
end

let
    fig = Figure()
    ax = Axis(fig[1,1], title="ρcore")
    lines!(ax, upf.ρcore.r, upf.ρcore.f)
    lines!(ax, psp8.ρcore.r, psp8.ρcore.f)
    ax = Axis(fig[2,1], title="ρval")
    lines!(ax, upf.ρval.r, upf.ρval.f)
    lines!(ax, psp8.ρval.r, psp8.ρval.f)
    fig
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
        upf_real_splined = interpolate_onto(upf_real, maximum_spacing)

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
