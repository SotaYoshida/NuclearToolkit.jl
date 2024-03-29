#using LocalCoverage;  generate_coverage("NuclearToolkit"; run_test=true)
using NuclearToolkit
using Test
using Printf

@testset "devDMD" begin

end

# @testset  "NuclearToolkit.jl" begin
#     @testset "generate NN potential" begin
#         @test make_chiEFTint(;fn_params="optional_parameters_srg_emn500n4lo_2n3n.jl")
#         @test make_chiEFTint(;fn_params="optional_parameters_nnlosat.jl",deuteron_check=true) 
#         @test make_chiEFTint(;fn_params="optional_parameters.jl")
#         @test make_chiEFTint(;fn_params="optional_parameters_snt.jl")
#     end

#     @testset "HFMBPT & VS-IMSRG calculations" begin
#         hw = 20; emax=2
#         nuc = "He4"; core = "He4"; vspace="p-shell"
#         nucs = ["He4"]
#         sntf = "tbme_em500n3lo_barehw20emax2.snt.bin"
#         ## HF-MBPT from snt/snt.bin
#         @testset "HFMBPT results under bare EM500,hw20,e2,nmesh50" begin
#             Eref = [1.493, -5.805, 0.395]
#             HFobj1 = hf_main(nucs,sntf,hw,emax;return_obj=true,verbose=true)
#             Es1 = [HFobj1.E0, HFobj1.EMP2, HFobj1.EMP3]
#             println("Eref $Eref")
#             println("Es1 $Es1")
#             @testset "HF" begin
#                 @test  (Es1[1] - Eref[1])^2 < 1.e-4
#             end
#             @testset "EMP2" begin
#                 @test  (Es1[2] - Eref[2])^2 < 1.e-4
#             end
#             @testset "EMP3" begin 
#                 @test  (Es1[3] - Eref[3])^2 < 1.e-4
#             end
#             @testset "snt & snt.bin must give identical results" begin
#                 tsntf = replace(sntf,".snt.bin" => ".snt")
#                 HFobj2 = hf_main(nucs,tsntf,hw,emax;return_obj=true)
#                 Es2 = [HFobj2.E0, HFobj2.EMP2, HFobj2.EMP3]    
#                 @test ((HFobj1.E0-HFobj2.E0)^2 + (HFobj1.EMP2-HFobj2.EMP2)^2 + (HFobj1.EMP3-HFobj2.EMP3)^2) < 1.e-6
#             end
#         end
#         Eref = -4.05225276 
#         @testset "IMSRG results under bare EM500,hw20,e2,nmesh50" begin
#             IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,return_obj=true)
#             Es = IMSRGobj.H.zerobody[1]
#             @test abs(Eref-Es[1]) < 1.e-6
#         end
#         @testset "VSIMSRG results under bare EM500,hw20,e2,nmesh50" begin
#             IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,corenuc=core,ref="nuc",valencespace=vspace,return_obj=true)
#             Es = IMSRGobj.H.zerobody[1]
#             @test abs(Eref - Es[1]) < 1.e-6
#         end
#     end
    
#     include("./ShellModel_test.jl")

#     @testset  "Testing EC on shell model" begin
#         target_nuc = "O18"
#         num_ev = 3
#         targetJ = 0
#         sntpath = "random_snts/"
#         spath = "wavsamples/"
#         if !isdir("wavsamples")
#             mkdir("wavsamples")
#         end
#         if !isdir("appwavs")
#             mkdir("appwavs")
#         end
    
#         println("sampling...")
#         mode = "sample"
#         num_ECsample = 1    
#         Hs = [ sntpath*"tmp_$i"*".snt" for i = 0:4 ]
#         @test prepEC(Hs,target_nuc,num_ev,num_ECsample,targetJ,mode;path_to_samplewav=spath,save_wav=true) == nothing
    
#         println("make TDmat...")
#         num_ECsamples = length(Hs)
#         mode = "TD"
#         @test prepEC(Hs,target_nuc,num_ev,num_ECsample,targetJ,mode;path_to_samplewav=spath) == nothing
    
#         println("solving EC")
#         write_appwav = true
#         @test solveEC(["usdb.snt"],target_nuc,[[targetJ,num_ev]];write_appwav=true,wpath="./wavsamples",tdmatpath="./tdmat/") == nothing       
       
#     end

#     @testset "making msnt file for Quantum Calculations" begin
#         @test main_trans_msnt("vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt","Li6") == nothing
#     end

#     @testset "testing HFMBPT/IMSRG for Rp2 with BetaCM=1.0" begin
#         nucs = ["He4"]
#         sntf = "tbme_emn500n4lo_2n3n_srg10.0hw20emax2.snt.bin"
#         hw = 20
#         emax = 2
#         IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,Operators=["Rp2"],return_obj=true)
#         #checking Rp2 value from IMSRG calculation
#         Rp2_ref = 1.559825^2
#         @test abs(IMSRGobj.ExpectationValues["Rp2"] - Rp2_ref)^2 < 1.e-6
#     end

#     @testset "genuine3NF" begin
#         test3NF()
#     end
#     @testset "2n3n calibration runnable?" begin
#         @test make_chiEFTint(;nucs=["He4"],itnum=3,optimizer="MCMC")
#         @test make_chiEFTint(;nucs=["He4"],itnum=3,optimizer="LHS")
#         @test make_chiEFTint(;nucs=["He4"],itnum=3,optimizer="BayesOpt")
#         #@test make_chiEFTint(;nucs=["He4"],itnum=6,optimizer="MCMC",MPIcomm=true)
#     end
    
#     @testset "Reading 3BME file from NuHamil code (by T. Miyagi)" begin
#         sntf   = "tbme_em500n3lo_barehw20emax2.snt.bin"
#         fn_3nf = "small_3BME_N3LOlnl_4_8_4.me3j.gz"
#         hw = 16
#         emax = 2; e3max = 4; e1max_file=4; e2max_file=8; e3max_file=4
#         HFobj = hf_main(["He6"],sntf,hw,emax;return_obj=true,
#                         e1max_file=e1max_file,e2max_file=e2max_file,e3max_file=e3max_file,e3max=e3max,fn_3nf=fn_3nf)
#         Es = [HFobj.E0, HFobj.EMP2, HFobj.EMP3]
#         Eref = [ 8.160870, -9.58767, -0.82059]
#         @test (Es[1] - Eref[1])^2 + (Es[2] - Eref[2])^2 + (Es[3] - Eref[3])^2 < 1.e-7
#     end
#     rm("tbme_em500n3lo_barehw20emax2.snt")
#     rm("tbme_em500n3lo_barehw20emax2.snt.bin")
#     rm("vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt")
#     rm("vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.msnt")
#     rm("He6_ckpot_j0.wav")
#     rm("He6_ckpot_j4.wav")
#     rm("flowOmega",recursive=true)
#     rm("A3files",recursive=true)
#     rm("tdmat",recursive=true)
#     rm("wavsamples",recursive=true)
#     rm("appwavs",recursive=true)
# end

