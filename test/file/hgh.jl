#TODO test hgh file parser
@testset "HGH File" begin
    @testset "c-q4.hgh" begin
        filename = "c-q4.hgh"
        file = load_psp_file(HGH_CASE_FILEPATHS[filename])

        @test occursin("c", lowercase(file.title))
        @test occursin("pade", lowercase(file.title))
        @test file.zion == [2, 2]
        @test file.rloc == 0.34883045
        @test file.cloc == [-8.51377110, 1.22843203]
        @test file.lmax == 1
        @test file.rp == [0.30455321, 0.2326773]
        @test file.h[1] == 9.52284179 * ones(1, 1)
        @test file.h[2] == zeros(0, 0)
    end

    @testset "ni-q18.hgh" begin
        filename = "ni-q18.hgh"
        file = load_psp_file(HGH_CASE_FILEPATHS[filename])

        @test occursin("ni", lowercase(file.title))
        @test occursin("pade", lowercase(file.title))
        @test file.zion == [4, 6, 8]
        @test file.rloc == 0.35000000
        @test file.cloc == [3.61031072, 0.44963832]
        @test file.lmax == 2
        @test file.rp == [0.24510489, 0.23474136, 0.21494950]
        @test file.h[1] == [[12.16113071, 3.51625420] [3.51625420, -9.07892931]]
        @test file.h[2] == [[-0.82062357, 2.54774737] [2.54774737, -6.02907069]]
        @test file.h[3] == -13.39506212 * ones(1, 1)
    end
end
