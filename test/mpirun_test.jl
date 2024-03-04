using NuclearToolkit
using Test

@testset "NuclearToolkit.jl" begin
    @testset "mpirun" begin
        @test make_chiEFTint(;nucs=["He4"],itnum=10,optimizer="MCMC",MPIcomm=true,fn_params="parameters/optional_parameters.jl")
    end
end
