@testset "PSP8 File" begin
    @testset "Internal data consistency" begin
        for filepath in values(PSP8_FILEPATHS)
            file = load_psp_file(filepath)

            @test format(file) == "PSP8"

            @test file.header.zatom >= file.header.zion
            @test file.header.pspcod == 8

            @test length(file.projectors) == file.header.lmax + 1
            if file.header.lloc < file.header.lmax
                @test isempty(file.projectors[lloc + 1])
                @test isempty(file.ekb[lloc + 1])
            end
            for projectors_l in file.projectors, projector in projectors_l
                @test length(projector) == file.header.mmax
            end
            for (l_plus_one, nproj_l) in enumerate(file.header.nproj)
                @test length(file.projectors[l_plus_one]) == nproj_l
                @test length(file.ekb[l_plus_one]) == nproj_l
            end

            @test length(file.v_local) == file.header.mmax

            core_densities = [file.rhoc, file.d_rhoc_dr, file.d2_rhoc_dr2, file.d3_rhoc_dr3,
                              file.d4_rhoc_dr4]
            if file.header.fchrg > 0
                @test !any(isnothing, core_densities)
                for density in core_densities
                    @test length(density) == file.header.mmax
                end
            else
                @test all(isnothing, core_densities)
            end

            if file.header.extension_switch in (2, 3)
                @test length(file.projectors_so) == file.header.lmax + 1
                @test isempty(file.projectors_so[0 + 1])
                for projectors_l in file.projectors_so, projector in projectors_l
                    @test length(projector) == file.header.mmax
                end
                @test length(file.ekb_so) == file.header.lmax + 1
                @test isempty(file.projectors_so[0 + 1])
                for (l_plus_one, nprojso_l) in enumerate(file.header.nprojso)
                    @test length(file.projectors_so[l_plus_one]) == nprojso_l
                    @test length(file.ekb_so[l_plus_one]) == nprojso_l
                end
            else
                @test isnothing(file.projectors_so)
                @test isnothing(file.ekb_so)
            end
        end
    end

    @testset "H.psp8" begin
        filename = "H.psp8"
        file = load_psp_file(PSP8_CASE_FILEPATHS[filename])

        @test file.header.zatom == 1.0000
        @test file.header.zion == 1.0000
        @test file.header.pspd == 170501
        @test file.header.pspcod == 8
        @test file.header.pspxc == 11
        @test file.header.lmax == 1
        @test file.header.lloc == 4
        @test file.header.mmax == 300
        @test file.header.r2well == 0.0
        @test file.header.rchrg == 2.99000000
        @test file.header.fchrg == 0.00000000
        @test file.header.qchrg == 0.00000000
        @test file.header.nproj == [2, 1]
        @test file.header.extension_switch == 1

        @test file.ekb[0 + 1][1] == -1.7011234716222E+00
        @test file.ekb[0 + 1][2] == -5.3448391737199E-01
        @test file.projectors[0 + 1][1][5] == 1.6367055221032E-01
        @test file.projectors[0 + 1][2][5] == 3.9125025705919E-01

        @test file.v_local[5] == -3.1387511173374E+00
    end

    @testset "Ti.psp8" begin
        filename = "Ti.psp8"
        file = load_psp_file(PSP8_CASE_FILEPATHS[filename])

        @test file.header.zatom == 22.0000
        @test file.header.zion == 12.0000
        @test file.header.pspd == 180429
        @test file.header.pspcod == 8
        @test file.header.pspxc == -116133
        @test file.header.lmax == 2
        @test file.header.lloc == 4
        @test file.header.mmax == 600
        @test file.header.r2well == 0.0
        @test file.header.rchrg == 5.99000000
        @test file.header.fchrg == 2.50000000
        @test file.header.qchrg == 0.00000000
        @test file.header.nproj == [2, 3, 2]
        @test file.header.extension_switch == 3
        @test file.header.nprojso == [0, 4, 3]

        @test file.ekb[0 + 1][1] == 1.2539360126722E+00
        @test file.ekb[0 + 1][2] == 1.1321100966516E+01
        @test file.projectors[0 + 1][1][5] == -1.7889864901133E-01
        @test file.projectors[0 + 1][2][5] == 3.5296270964743E-01

        @test file.v_local[5] == -1.9130317133247E+01

        @test file.ekb_so[1 + 1][1] == 6.9968670783987E-02
        @test file.ekb_so[1 + 1][2] == -3.0732828917905E-02
        @test file.ekb_so[1 + 1][3] == 4.6090926948408E-03
        @test file.ekb_so[1 + 1][4] == -2.0693922436038E-04
        @test file.projectors_so[1 + 1][1][5] == 1.7380194666161E-02
        @test file.projectors_so[1 + 1][2][5] == -6.8327678530662E-02
        @test file.projectors_so[1 + 1][3][5] == 2.3034649777788E-02
        @test file.projectors_so[1 + 1][4][5] == -6.2709563805074E-03
    end

    @testset "Zn.psp8" begin
        filename = "Zn.psp8"
        file = load_psp_file(PSP8_CASE_FILEPATHS[filename])

        @test file.header.zatom == 30.0000
        @test file.header.zion == 20.0000
        @test file.header.pspd == 180112
        @test file.header.pspcod == 8
        @test file.header.pspxc == -116133
        @test file.header.lmax == 2
        @test file.header.lloc == 4
        @test file.header.mmax == 600
        @test file.header.r2well == 0.0
        @test file.header.rchrg == 5.99000000
        @test file.header.fchrg == 2.50000000
        @test file.header.qchrg == 0.00000000
        @test file.header.nproj == [2, 2, 2]
        @test file.header.extension_switch == 1

        @test file.ekb[0 + 1][1] == 2.8308354242488E+00
        @test file.ekb[0 + 1][2] == -1.5191074296972E+00
        @test file.projectors[0 + 1][1][5] == 3.5070387789440E-01
        @test file.projectors[0 + 1][2][5] == 1.2235313804112E-01

        @test file.v_local[5] == -3.1208323049679E+01

        @test file.rhoc[7] == 1.2735292737007E+02
        @test file.d_rhoc_dr[7] == -1.3288155269898E+02
        @test file.d2_rhoc_dr2[7] == -2.0848786910949E+03
        @test file.d3_rhoc_dr3[7] == 6.3716538707879E+03
        @test file.d4_rhoc_dr4[7] == 9.6377063748116E+04
    end
end
