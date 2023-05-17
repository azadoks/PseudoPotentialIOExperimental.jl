function load_psp8(filepath::AbstractString)
    open(filepath, "r") do io
        return load_psp8(io)
    end
end

function load_psp8(io::IO)
    @warn "load_psp8 is deprecated, use load_psp_file and the new PspFile structures"
    new_psp8 = Psp8File(io)
    old_psp8 = Dict()

    old_psp8["header"] = Dict("z_atom" => new_psp8.header.zatom,
                              "z_valence" => new_psp8.header.zion,
                              "generation_day" => parse(Int,
                                                        string(new_psp8.header.pspd)[1:2]),
                              "generation_month" => parse(Int,
                                                          string(new_psp8.header.pspd)[3:4]),
                              "generation_year" => parse(Int,
                                                         string(new_psp8.header.pspd)[5:6]),
                              "format_version" => new_psp8.header.pspcod,
                              "xc" => new_psp8.header.pspxc,
                              "l_max" => new_psp8.header.lmax,
                              "l_local" => new_psp8.header.lloc,
                              "mesh_size" => new_psp8.header.mmax,
                              "r2well" => new_psp8.header.r2well,
                              "rchrg" => new_psp8.header.rchrg,
                              "fchrg" => new_psp8.header.fchrg,
                              "qchrg" => new_psp8.header.qchrg,
                              "core_correction" => new_psp8.header.fchrg > 0,
                              "number_of_proj" => new_psp8.header.nproj,
                              "extension_switch" => new_psp8.header.extension_switch,
                              "has_so" => !isnothing(new_psp8.header.nprojso))
    if old_psp8["header"]["has_so"]
        old_psp8["header"]["number_of_proj_so"] = new_psp8.header.nprojso
    end
    @assert old_psp8["header"]["format_version"] == 8

    old_psp8["radial_grid"] = new_psp8.rgrid
    old_psp8["local_potential"] = new_psp8.v_local

    old_psp8["beta_projectors"] = new_psp8.projectors
    old_psp8["ekb"] = map(new_psp8.ekb) do ekb_l
        ekb_mat = zeros(length(ekb_l), length(ekb_l))
        for i in eachindex(ekb_l)
            ekb_mat[i, i] = ekb_l[i]
        end
        return ekb_mat
    end

    if !isnothing(new_psp8.header.nprojso)
        old_psp8["spin_orbit"] = Dict("beta_projectors" => new_psp8.projectors_so,
                                      "ekb" => new_psp8.ekb_so)
    else
        old_psp8["spin_orbit"] = Dict()
    end

    if new_psp8.header.rchrg > 0
        old_psp8["nlcc"] = Dict("core_charge_density" => new_psp8.rhoc,
                                "first_derivative" => new_psp8.d_rhoc_dr,
                                "second_derivative" => new_psp8.d2_rhoc_dr2,
                                "third_derivative" => new_psp8.d3_rhoc_dr3,
                                "fourth_derivative" => new_psp8.d4_rhoc_dr4)
    else
        old_psp8["nlcc"] = Dict()
    end

    return old_psp8
end