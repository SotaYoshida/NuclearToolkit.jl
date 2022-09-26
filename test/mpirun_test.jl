import Pkg; Pkg.add("NuclearToolkit")
using NuclearToolkit
using .NuclearToolkit
using Test

@testset "NuclearToolkit.jl" begin
    @testset "mpirun" begin
        @test make_chiEFTint(;nucs=["He4"],itnum=10,optimizer="MCMC",MPIcomm=true)
    end
end
