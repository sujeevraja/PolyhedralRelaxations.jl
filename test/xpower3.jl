@testset "x^3 test" begin
    model = PR.main()
    @test model.x_index == 1
    @test model.y_index == 2
end