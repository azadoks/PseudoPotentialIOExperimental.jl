@testset "Norm conserving" begin
    @testset "From PSP8" begin
        @testset "H.psp8" begin
            psp = load_psp(PSP8_CASE_FILEPATHS["H.psp8"])
            @test isa(psp, NormConservingPsP{Float64,RealSpace})
        end
        @testset "Ti.psp8" begin
            @test_throws ErrorException("Fully relativistic pseudos are not supported") load_psp(PSP8_CASE_FILEPATHS["Ti.psp8"])
        end
        @testset "Zn.psp8" begin
            psp = load_psp(PSP8_CASE_FILEPATHS["Zn.psp8"])
            @test isa(psp, NormConservingPsP{Float64,RealSpace})
        end
    end

    @testset "From UPF" begin
        @testset "B_pbe_v1.01.uspp.F.UPF" begin
            file = load_psp_file(UPF1_CASE_FILEPATHS["B_pbe_v1.01.uspp.F.UPF"])
            @test_throws ErrorException("Provided `UpfFile` is not a norm-conserving pseudo") NormConservingPsP(file)
        end
        @testset "Dy.GGA-PBE-paw-v1.0.UPF" begin
            file = load_psp_file(UPF2_CASE_FILEPATHS["Dy.GGA-PBE-paw-v1.0.UPF"])
            @test_throws ErrorException("Provided `UpfFile` is not a norm-conserving pseudo") NormConservingPsP(file)
        end
        @testset "Mg.upf" begin
            @test_throws ErrorException("Fully relativistic pseudos are not supported") load_psp(UPF2_CASE_FILEPATHS["Mg.upf"])
        end
    end
end
