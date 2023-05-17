#TODO add tests for augmentation real, augmentation fourier-real agreement, etc.
@testset "Ultrasoft" begin
    @testset "From UPF" begin
        @testset "B_pbe_v1.01.uspp.F.UPF" begin
            file = load_psp_file(UPF1_CASE_FILEPATHS["B_pbe_v1.01.uspp.F.UPF"])
            @test isa(load_psp(file), UltrasoftPsP{Float64,RealSpace})
        end
        @testset "Si.pbe-n-rrkjus_psl.1.0.0.UPF" begin
            file = load_psp_file(UPF2_CASE_FILEPATHS["Si.pbe-n-rrkjus_psl.1.0.0.UPF"])
            @test isa(load_psp(file), UltrasoftPsP{Float64,RealSpace})
        end
        @testset "Dy.GGA-PBE-paw-v1.0.UPF" begin
            file = load_psp_file(UPF2_CASE_FILEPATHS["Dy.GGA-PBE-paw-v1.0.UPF"])
            @test_throws ErrorException("Provided `UpfFile` is not an ultrasoft pseudo") UltrasoftPsP(file)
        end
        @testset "Mg_nc-fr-04_pbesol_stringent.upf" begin
            file = load_psp_file(UPF2_CASE_FILEPATHS["Mg.upf"])
            @test_throws ErrorException("Provided `UpfFile` is not an ultrasoft pseudo") UltrasoftPsP(file)
        end
    end
end
