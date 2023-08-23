@testset "UPF" begin
    VECTOR_KEYS = ["radial_grid", "core_charge_density", "local_potential", "total_charge_density"]
    HEADER_KEYS = ["number_of_proj", "core_correction", "element", "pseudo_type", "z_valence",
                   "mesh_size", "number_of_wfc"]

    pseudos = load_upf.(UPF_FILEPATHS)

    @testset "mesh" begin
        for upf in pseudos
            @test length(upf["radial_grid"]) == upf["header"]["mesh_size"]
            @test length(upf["radial_grid_derivative"]) == upf["header"]["mesh_size"]
            r_params = upf["radial_grid_parameters"]
            if r_params["mesh_type"] == "log_1"
                r = PseudoPotentialIOExperimental.LogMeshWithoutZero{Float64}(
                    r_params["a"], r_params["b"], upf["radial_grid"][1],
                    upf["radial_grid"][end], length(upf["radial_grid"]))
                @test all(upf["radial_grid"] .≈ r)
            elseif r_params["mesh_type"] == "log_2"
                r = PseudoPotentialIOExperimental.LogMeshWithZero{Float64}(
                    r_params["a"], r_params["b"], upf["radial_grid"][2],
                    upf["radial_grid"][3], length(upf["radial_grid"]))
                @test all(upf["radial_grid"] .≈ r)
            elseif r_params["mesh_type"] == "linear"
                @test all(isapprox.(round.(diff(upf["radial_grid"]), digits=4),
                                    upf["radial_grid_derivative"][1:(end - 1)]))
            end
        end
    end

    @testset "nlcc" begin
        for upf in pseudos
            if upf["header"]["core_correction"]
                @test length(upf["core_charge_density"]) == upf["header"]["mesh_size"]
            else
                @test !haskey(upf, "core_charge_density")
            end
        end
    end

    @testset "local" begin
        for upf in pseudos
            @test length(upf["local_potential"]) == upf["header"]["mesh_size"]
        end
    end

    @testset "beta" begin
        for upf in pseudos
            @test length(upf["beta_projectors"]) == upf["header"]["number_of_proj"]
            for proj in upf["beta_projectors"]
                @test length(proj["radial_function"]) == proj["cutoff_radius_index"]
                @test proj["index"] <= upf["header"]["number_of_proj"]
                @test proj["angular_momentum"] <= upf["header"]["l_max"]
                if upf["header"]["has_so"]
                    @test haskey(proj, "total_angular_momentum")
                end
            end
        end
    end

    @testset "Dij" begin
        for upf in pseudos
            @test ndims(upf["D_ion"]) == 2
            @test size(upf["D_ion"], 1) == size(upf["D_ion"], 2) == upf["header"]["number_of_proj"]
        end
    end

    @testset "pswfc" begin
        for upf in pseudos
            @test length(upf["atomic_wave_functions"]) == upf["header"]["number_of_wfc"]
            @test (sum(wfc -> wfc["occupation"], upf["atomic_wave_functions"]; init=0) <=
                   upf["header"]["z_valence"])
            for wfc in upf["atomic_wave_functions"]
                @test length(wfc["radial_function"]) == upf["header"]["mesh_size"]
                # @test wfc["angular_momentum"] <= clamp(upf["header"]["l_max"], 0:5)  Doesn't seem to be true?
                if upf["header"]["has_so"]
                    @test haskey(wfc, "principal_quantum_number")
                    @test haskey(wfc, "total_angular_momentum")
                end
            end
        end
    end

    @testset "rhoatom" begin
        for upf in pseudos
            @test length(upf["total_charge_density"]) == upf["header"]["mesh_size"]
        end
    end
end
