using NuclearToolkit
using Test

@testset "NuclearToolkit.jl" begin
    
    # generate NN potential
    @test make_chiEFTint()
    @test make_chiEFTint(;fn_params="optional_parameters_snt.jl")
    
    # HFMBPT & VS-IMSRG calculations
    hw = 20; emax=2
    nuc = "He4"; core = "He4"; vspace="p-shell"
    nucs = ["He4"]

    ## HF-MBPT from snt/snt.bin
    sntf = "tbme_em500n3lo_barehw20emax2.snt.bin"
    HFobj1 = hf_main(nucs,sntf,hw,emax;return_HFobj=true)
    Es1 = [HFobj1.E0, HFobj1.EMP2, HFobj1.EMP3]
    sntf = "tbme_em500n3lo_barehw20emax2.snt"
    HFobj2 = hf_main(nucs,sntf,hw,emax;return_HFobj=true)
    Es2 = [HFobj2.E0, HFobj2.EMP2, HFobj2.EMP3]
    @test ((HFobj1.E0-HFobj2.E0)^2 + (HFobj1.EMP2-HFobj2.EMP2)^2 + (HFobj1.EMP3-HFobj2.EMP3)^2) < 1.e-6

    ## IMSRG
    sntf = "tbme_em500n3lo_barehw20emax2.snt.bin"
    @test hf_main(nucs,sntf,hw,emax;doIMSRG=true,corenuc=core,ref="nuc",valencespace=vspace)

    ## shell model calculation
    vs_sntf = "vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt";  n_eigen=2;targetJ=[]
    @test main_sm(vs_sntf,"Be8",n_eigen,targetJ)
    
end
