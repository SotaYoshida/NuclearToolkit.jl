using NuclearToolkit
using Test

@testset  "NuclearToolkit.jl" begin
    @testset "generate NN potential" begin
        @test make_chiEFTint(;fn_params="optional_parameters.jl")
        @test make_chiEFTint(;fn_params="optional_parameters_snt.jl")
    end

    @testset "HFMBPT & VS-IMSRG calculations" begin
        hw = 20; emax=2
        nuc = "He4"; core = "He4"; vspace="p-shell"
        nucs = ["He4"]
        sntf = "tbme_em500n3lo_barehw20emax2.snt.bin"
        ## HF-MBPT from snt/snt.bin
        @testset "HFMBPT results under bare EM500,hw20,e2,nmesh50" begin
            Eref = [1.493, -5.805, 0.395]
            HFobj1 = hf_main(nucs,sntf,hw,emax;return_obj=true,verbose=true)
            Es1 = [HFobj1.E0, HFobj1.EMP2, HFobj1.EMP3]
            println("Eref $Eref")
            println("Es1 $Es1")
            @testset "HF" begin
                @test  (Es1[1] - Eref[1])^2 < 1.e-4
            end
            @testset "EMP2" begin
                @test  (Es1[2] - Eref[2])^2 < 1.e-4
            end
            @testset "EMP3" begin 
                @test  (Es1[3] - Eref[3])^2 < 1.e-4
            end
            @testset "snt & snt.bin must give identical results" begin
                tsntf = replace(sntf,".snt.bin" => ".snt")
                HFobj2 = hf_main(nucs,tsntf,hw,emax;return_obj=true)
                Es2 = [HFobj2.E0, HFobj2.EMP2, HFobj2.EMP3]    
                @test ((HFobj1.E0-HFobj2.E0)^2 + (HFobj1.EMP2-HFobj2.EMP2)^2 + (HFobj1.EMP3-HFobj2.EMP3)^2) < 1.e-6
            end
        end
        Eref = -4.05225276 
        @testset "IMSRG results under bare EM500,hw20,e2,nmesh50" begin
            IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,return_obj=true)
            Es = IMSRGobj.H.zerobody[1]
            @test abs(Eref-Es[1]) < 1.e-6
        end
        @testset "VSIMSRG results under bare EM500,hw20,e2,nmesh50" begin
            IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,corenuc=core,ref="nuc",valencespace=vspace,return_obj=true)
            Es = IMSRGobj.H.zerobody[1]
            @test abs(Eref - Es[1]) < 1.e-6
        end
        @testset "shell model calculation" begin
            Eref = [ -10.720, -8.410]
            vs_sntf = "vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt";  n_eigen=2;targetJ=[]
            Ens = main_sm(vs_sntf,"Be8",n_eigen,targetJ)
            @test ((Eref[1]-Ens[1])^2 + (Eref[2] - Ens[2])^2) < 1.e-6
        end
    end
    @testset "making msnt file for Quantum Calculations" begin
        @test main_trans_msnt("vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt","Li6") == nothing
    end
    # @testset "2n3n calibration runnable?" begin
    #     @test make_chiEFTint(;nucs=["He4"],itnum=5,optimizer="MCMC")
    #     @test make_chiEFTint(;nucs=["He4"],itnum=5,optimizer="LHS")
    #     @test make_chiEFTint(;nucs=["He4"],itnum=5,optimizer="BayesOpt")
    #     @test make_chiEFTint(;nucs=["He4"],itnum=6,optimizer="MCMC",MPIcomm=true)
    # end
    
    @testset "Reading 3BME file from NuHamil code (by T. Miyagi)" begin
        sntf   = "tbme_em500n3lo_barehw20emax2.snt.bin"
        fn_3nf = "small_3BME_N3LOlnl_4_8_4.me3j.gz"
        hw = 16
        emax = 2; e3max = 4; e1max_file=4; e2max_file=8; e3max_file=4
        HFobj = hf_main(["He6"],sntf,hw,emax;return_obj=true,
                        e1max_file=e1max_file,e2max_file=e2max_file,e3max_file=e3max_file,e3max=e3max,fn_3nf=fn_3nf)
        Es = [HFobj.E0, HFobj.EMP2, HFobj.EMP3]
        Eref = [ 8.160870, -9.58767, -0.82059]
        @test (Es[1] - Eref[1])^2 + (Es[2] - Eref[2])^2 + (Es[3] - Eref[3])^2 < 1.e-7
    end
    rm("tbme_em500n3lo_barehw20emax2.snt")
    rm("tbme_em500n3lo_barehw20emax2.snt.bin")
    rm("vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt")
    rm("vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.msnt")
    rm("flowOmega",recursive=true)
end
