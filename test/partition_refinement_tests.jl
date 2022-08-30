@testset "test partition refinement strategies" begin 

    @testset "bisecting all sub-intervals" begin 
        PR.silence!()
        partition = range(0.0,10.0; step=1.0) |> collect 
        refinement_info = refine_partition!(partition, 0.0; 
            refinement_type = :bisect_all)
        @test refinement_info.num_additional_points == 10 
    end 
end 