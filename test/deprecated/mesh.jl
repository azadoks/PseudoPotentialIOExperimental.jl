import PseudoPotentialIOExperimental: RadialMesh, ArbitraryMesh, UniformMesh
import PseudoPotentialIOExperimental: LogMeshWithZero, LogMeshWithoutZero
import PseudoPotentialIOExperimental: deriv

@testset "Mesh type guess" begin
    @testset "[UniformMesh] Mg.upf" begin
        psp = load_psp_file(UPF2_CASE_FILEPATHS["Mg.upf"])
        mesh = RadialMesh(psp.mesh.r, psp.mesh.rab)
        @test typeof(mesh) <: UniformMesh
        @test mesh.a ≈ 0.01
        @test mesh.b ≈ 0.00
        @test all(mesh .≈ psp.mesh.r)
        @test all(deriv(mesh) .≈ psp.mesh.rab)
        @test all(diff(mesh) .≈ 0.01)
        @test length(mesh) == length(psp.mesh.r)
        @test length(deriv(mesh)) == length(psp.mesh.rab)
        @test length(diff(mesh)) == length(psp.mesh.r) - 1
    end

    @testset "[UniformMesh] Fe.psp8" begin
        psp = load_psp_file(PSP8_CASE_FILEPATHS["Fe.psp8"])
        mesh = RadialMesh(psp.rgrid)
        @test typeof(mesh) <: UniformMesh
        @test mesh.a ≈ 0.01
        @test mesh.b ≈ 0.00
        @test all(mesh .≈ psp.rgrid)
        @test all(deriv(mesh) .≈ 0.01)
        @test all(diff(mesh) .≈ 0.01)
        @test length(mesh) == length(psp.rgrid)
        @test length(deriv(mesh)) == length(psp.rgrid)
        @test length(diff(mesh)) == length(psp.rgrid) - 1
    end

    LogMeshWithoutZero_filenames = ["Al.pbe-n-kjpaw_psl.1.0.0.UPF", "He.pbe-hgh.UPF"]
    @testset "[LogMeshWithoutZero] $(filename)" for filename in LogMeshWithoutZero_filenames
        psp = load_psp_file(UPF2_CASE_FILEPATHS[filename])
        mesh = RadialMesh(psp.mesh.r, psp.mesh.rab)
        @test typeof(mesh) <: LogMeshWithoutZero
        @test mesh.a ≈ psp.mesh.dx
        @test mesh.b ≈ exp(psp.mesh.xmin) / psp.mesh.zmesh
        @test all(mesh .≈ psp.mesh.r)
        @test all(deriv(mesh) .≈ psp.mesh.rab)
        @test all(diff(mesh) .≈ diff(psp.mesh.r))
        @test length(mesh) == length(psp.mesh.r)
        @test length(deriv(mesh)) == length(psp.mesh.rab)
        @test length(diff(mesh)) == length(psp.mesh.r) - 1
    end

    LogMeshWithZero_filenames = ["ag_lda_v1.4.uspp.F.UPF", "B_pbe_v1.01.uspp.F.UPF"]
    @testset "[LogMeshWithZero] $(filename)" for filename in LogMeshWithZero_filenames
        psp = load_psp_file(UPF1_CASE_FILEPATHS[filename])
        mesh = RadialMesh(psp.mesh.r, psp.mesh.rab)
        @test typeof(mesh) <: LogMeshWithZero
        @test all(mesh .≈ psp.mesh.r)
        @test all(deriv(mesh) .≈ psp.mesh.rab)
        @test all(diff(mesh) .≈ diff(psp.mesh.r))
        @test length(mesh) == length(psp.mesh.r)
        @test length(deriv(mesh)) == length(psp.mesh.rab)
        @test length(diff(mesh)) == length(psp.mesh.r) - 1
    end

    @testset "[ArbitraryMesh]" begin
        r = 0:rand():10
        mesh = ArbitraryMesh(r)
        @test typeof(mesh) <: ArbitraryMesh
        @test all(mesh .== r)
        @test all(diff(mesh) .== diff(r))
        @test length(mesh) == length(r)
        @test length(diff(mesh)) == length(diff(r))
    end
end
