import PseudoPotentialIOExperimental: resolve_family

@testset "Load PsP file" begin
    psp_filepath_tuples = [("hgh_lda_hgh", "si-q4.hgh"),  # HGH
                           ("pd_nc_sr_pbe_standard_0.4.1_upf", "Si.upf"),  # UPF NC SR
                           ("pd_nc_sr_pbe_standard_0.4.1_psp8", "Si.psp8"),  # PSP8 SR
                           ("pd_nc_fr_pbe_standard_0.4_upf", "Si.upf"),  # UPF NC FR
                           ("pd_nc_fr_pbe_standard_0.4_psp8", "Si.psp8"),  # PSP8 FR
                           ("gbrv_lda_1.5_upf", "si_lda_v1.uspp.F.UPF"),  # UPF USPP SR
                           ("sssp_pbe_precision_1.1.2_upf", "La.GGA-PBE-paw-v1.0.UPF")]  # UPF PAW SR
    for (familyname, filename) in psp_filepath_tuples
        familypath = resolve_family(familyname)
        filepath = joinpath(familypath, filename)
        @test load_psp_file(familyname, filename) == load_psp_file(filepath)
    end
    @test_throws ArgumentError load_psp_file("runtests.jl")  # Bad file extension
    @test_throws SystemError load_psp_file("abc.upf")  # File does not exist
end

@testset "Load family PsP files" begin
    for family in TEST_FAMILIES
        family = load_family_psp_files(family)
        @test eltype(family) <: PsPFile
    end

    (root, dirs, _) = first(walkdir("./"))
    for dir in dirs
        family = load_family_psp_files(joinpath(root, dir))
        @test isempty(family)
    end
end
