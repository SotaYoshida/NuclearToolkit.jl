@testset "2n3n calibration runnable?" begin
    @test make_chiEFTint(;nucs=["He4"],itnum=3,optimizer="MCMC")
    @test make_chiEFTint(;nucs=["He4"],itnum=3,optimizer="LHS")
    @test make_chiEFTint(;nucs=["He4"],itnum=3,optimizer="BayesOpt")
    #@test make_chiEFTint(;nucs=["He4"],itnum=6,optimizer="MCMC",MPIcomm=true)
end
