@testset "test partition refinement strategies" begin
    @testset "bisecting all sub-intervals" begin
        PR.silence!()
        partition = range(0.0, 10.0; step = 1.0) |> collect
        refinement_info =
            refine_partition!(partition, 0.0; refinement_type = :bisect_all)
        @test refinement_info.num_additional_points == 10
        @test refinement_info.refined_largest_partition == false
    end

    @testset "at_point refinement strategy" begin
        PR.silence!()
        partition = range(0.0, 10.0; step = 1.0) |> collect
        refinement_info =
            refine_partition!(partition, 1E-4; refinement_type = :at_point)
        @test refinement_info.num_additional_points == 0
        @test refinement_info.refined_largest_partition == false

        refinement_info =
            refine_partition!(partition, 0.5; refinement_type = :at_point)
        @test refinement_info.num_additional_points == 1
        @test refinement_info.refined_largest_partition == false
    end

    @testset "bisect refinement strategy" begin
        PR.silence!()
        partition = range(0.0, 10.0; step = 1.0) |> collect
        refinement_info =
            refine_partition!(partition, 1E-4; refinement_type = :bisect)
        @test refinement_info.num_additional_points == 1
        @test refinement_info.refined_largest_partition == false
        @test 0.5 in partition

        refinement_info =
            refine_partition!(partition, 0.5; refinement_type = :bisect)
        @test refinement_info.num_additional_points == 0
        @test refinement_info.refined_largest_partition == false

        partition = [0.0, 1E-4, 3]
        refinement_info =
            refine_partition!(partition, 1E-5; refinement_type = :bisect)
        @test refinement_info.num_additional_points == 1
        @test refinement_info.refined_largest_partition == true

        partition = [0.0, 1E-4]
        refinement_info =
            refine_partition!(partition, 1E-5; refinement_type = :bisect)
        @test refinement_info.num_additional_points == 0
        @test refinement_info.refined_largest_partition == false
    end

    @testset "non_uniform refinement strategy" begin
        PR.silence!()
        partition = range(0.0, 10.0; step = 1.0) |> collect
        refinement_info =
            refine_partition!(partition, 1E-4; refinement_type = :non_uniform)
        @test refinement_info.num_additional_points == 1
        @test refinement_info.refined_largest_partition == false
        @test 9.0E-5 in partition

        partition = range(0.0, 10.0; step = 1.0) |> collect
        refinement_info =
            refine_partition!(partition, 0.5; refinement_type = :non_uniform)
        @test refinement_info.num_additional_points == 2
        @test refinement_info.refined_largest_partition == false
        @test 0.45 in partition
        @test 0.55 in partition

        partition = [0.0, 1E-4, 3]
        refinement_info =
            refine_partition!(partition, 1E-5; refinement_type = :non_uniform)
        @test refinement_info.num_additional_points == 1
        @test refinement_info.refined_largest_partition == true

        partition = [0.0, 1E-4, 3]
        refinement_info = refine_partition!(
            partition,
            1E-5;
            refinement_type = :non_uniform,
            refine_largest = false,
        )
        @test refinement_info.num_additional_points == 0
        @test refinement_info.refined_largest_partition == false
    end
end
