@testset "generate NN potential" begin
    @test make_chiEFTint(;fn_params="parameters/optional_parameters_srg_emn500n4lo_2n3n.jl")
    @test make_chiEFTint(;fn_params="parameters/optional_parameters_nnlosat.jl",deuteron_check=true) 
    @test make_chiEFTint(;fn_params="parameters/optional_parameters.jl")
    @test make_chiEFTint(;fn_params="parameters/optional_parameters_snt.jl")
end