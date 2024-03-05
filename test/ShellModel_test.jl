@testset "shell model calculation" begin
    @testset "shell model results with VS-IMSRG interaction" begin        
        Eref = [ -10.720, -8.410]
        vs_sntf = "vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt";  n_eigen=2;targetJ=[]
        Ens = main_sm(vs_sntf,"Be8",n_eigen,targetJ)
        @test ((Eref[1]-Ens[1])^2 + (Eref[2] - Ens[2])^2) < 1.e-6
    end
    @testset "shell model results with CKpot" begin
        vs_sntf = "interaction_file/ckpot.snt"
        Eref = [-31.119,-27.300, -19.162, -18.249, -16.722, -14.925, -14.517, -14.017, -13.951, -13.478]
        Ens = main_sm(vs_sntf,"Be8",10,[];q=2,is_block=true)
        for i = 1:10
            @test (Eref[i] - Ens[i])^2 < 1.e-6 
        end
    end
    @testset "specifying truncation for occupation numbers" begin
        tdict = Dict("n1s1" => [0,1], "n0d5" => [1])
        main_sm("interaction_file/usdb.snt","O18", 1, [0];truncation_scheme="jocc", truncated_jocc=tdict)       
    end
    @testset "pnsystem with pn-pair ansatz" begin 
        vs_sntf = "interaction_file/ckpot.snt"
        main_sm("interaction_file/ckpot.snt","Be8", 1, [0];truncation_scheme="pn-pair")
    end

    @testset "transition check" begin
        vs_sntf = "interaction_file/ckpot.snt"
        main_sm(vs_sntf,"Li6",1,[0]; calc_entropy = true, visualize_occ = true)       
        main_sm(vs_sntf,"He6",1,[0];save_wav=true)
        main_sm(vs_sntf,"He6",1,[4];calc_moment=true, save_wav=true)
        jl2 = 0; jr2 = 4
        in_wfs = ["He6_ckpot_j0.wav","He6_ckpot_j4.wav"]
        transit_main(vs_sntf,"He6",jl2,jr2,in_wfs; calc_EM=true)
    end
    @testset "making msnt file for Quantum Calculations" begin
        vs_sntf = "interaction_file/ckpot.snt"
        @test main_trans_msnt(vs_sntf,"Li6") == nothing
    end    
end

@testset  "Testing EC on shell model" begin
    target_nuc = "O18"
    num_ev = 3
    targetJ = 0
    sntpath = "interaction_file/random_snts/"
    spath = "wavsamples/"
    if !isdir("wavsamples")
        mkdir("wavsamples")
    end
    if !isdir("appwavs")
        mkdir("appwavs")
    end

    println("sampling...")
    mode = "sample"
    num_ECsample = 1    
    Hs = [ sntpath*"tmp_$i"*".snt" for i = 0:4 ]
    @test prepEC(Hs,target_nuc,num_ev,num_ECsample,targetJ,mode;path_to_samplewav=spath,save_wav=true) == nothing

    println("make TDmat...")
    num_ECsample = length(Hs)
    mode = "TD"
    @test prepEC(Hs,target_nuc,num_ev,num_ECsample,targetJ,mode;path_to_samplewav=spath) == nothing

    println("solving EC")
    write_appwav = true
    @test solveEC(["interaction_file/usdb.snt"],target_nuc,[[targetJ,num_ev]];write_appwav=true,wpath="./wavsamples",tdmatpath="./tdmat/") == nothing       
    
end
