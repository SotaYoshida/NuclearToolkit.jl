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
            Eref = [1.477089, -11.88582, -3.38988]
            HFobj1 = hf_main(nucs,sntf,hw,emax;return_obj=true)
            Es1 = [HFobj1.E0, HFobj1.EMP2, HFobj1.EMP3]
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
        @testset "IMSRG results under bare EM500,hw20,e2,nmesh50" begin
            IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,return_obj=true)
            Es = IMSRGobj.H.zerobody[1]
            @test abs(-17.12108879-Es[1]) < 1.e-6
        end
        @testset "VSIMSRG results under bare EM500,hw20,e2,nmesh50" begin
            IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,corenuc=core,ref="nuc",valencespace=vspace,return_obj=true)
            Es = IMSRGobj.H.zerobody[1]
            @test abs(-17.12108879 - Es[1]) < 1.e-6
        end
        @testset "shell model calculation" begin
            vs_sntf = "vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt";  n_eigen=2;targetJ=[]
            Ens = main_sm(vs_sntf,"Be8",n_eigen,targetJ)
            @test ((-23.259-Ens[1])^2 + (-20.113 - Ens[2])^2) < 1.e-6
        end
    end

    @testset "2n3n calibration runnable?" begin
        @test make_chiEFTint(;nucs=["He4"],itnum=5,optimizer="MCMC")
        @test make_chiEFTint(;nucs=["He4"],itnum=5,optimizer="LHS")
        @test make_chiEFTint(;nucs=["He4"],itnum=5,optimizer="BayesOpt")
        @test make_chiEFTint(;nucs=["He4"],itnum=6,optimizer="MCMC",MPIcomm=true)
    end
end
