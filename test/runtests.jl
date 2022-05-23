using NuclearToolkit
using Test

@testset "NuclearToolkit.jl" begin
    
    ### generate NN potential
    @test make_chiEFTint()

    ### HFMBPT & VS-IMSRG calculation 
    hw = 20; emax=2
    nuc = "He4"; core = "He4"; vspace="p-shell"
    sntf = "tbme_em500n3lo_barehw20emax2.snt.bin"
    @test hf_main([nuc],sntf,hw,emax;verbose=false,doIMSRG=true,corenuc=core,ref="nuc",valencespace=vspace)

    ## shell model calculation
    vs_sntf = "vsimsrg_p-shell_coreHe4refHe4_hw20e2_Delta0.0.snt";  n_eigen=5;targetJ=[]
    @test main_sm(vs_sntf,"Be8",n_eigen,targetJ)
    
end
