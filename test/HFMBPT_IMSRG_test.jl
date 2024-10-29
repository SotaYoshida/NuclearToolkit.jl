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
            HFobj2 = hf_main(nucs,tsntf,hw,emax;return_obj=true,fn_params="parameters/optional_parameters.jl")
            Es2 = [HFobj2.E0, HFobj2.EMP2, HFobj2.EMP3]    
            @test ((HFobj1.E0-HFobj2.E0)^2 + (HFobj1.EMP2-HFobj2.EMP2)^2 + (HFobj1.EMP3-HFobj2.EMP3)^2) < 1.e-6
        end
    end
    Eref = -4.05225276 
    @testset "IMSRG results under bare EM500,hw20,e2,nmesh50" begin
        IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,return_obj=true,fn_params="parameters/optional_parameters.jl")
        Es = IMSRGobj.H.zerobody[1]
        @test abs(Eref-Es[1]) < 1.e-6
    end
    @testset "VSIMSRG results under bare EM500,hw20,e2,nmesh50" begin
        IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,corenuc=core,ref="nuc",valencespace=vspace,return_obj=true,fn_params="parameters/optional_parameters.jl")
        Es = IMSRGobj.H.zerobody[1]
        @test abs(Eref - Es[1]) < 1.e-6
    end

    @testset "testing HFMBPT/IMSRG for Rp2 with BetaCM=1.0" begin
        nucs = ["He4"]
        sntf = "tbme_emn500n4lo_2n3n_srg10.0hw20emax2.snt.bin"
        hw = 20
        emax = 2
        IMSRGobj = hf_main(nucs,sntf,hw,emax;doIMSRG=true,Operators=["Rp2"],return_obj=true,fn_params="parameters/optional_parameters.jl")
        Rp2_ref = 1.559825^2
        @test abs(IMSRGobj.ExpectationValues["Rp2"] - Rp2_ref)^2 < 1.e-6
    end

    @testset "testing DMD" begin      
        nuc = "He4"; emax = 2; hw = 20
        s_pred = [30.0, 50.0]
        ds = 0.5; smin = 15.0; smax=20.0
        # generating Omega 
        sntf = "tbme_emn500n4lo_2n3n_srg10.0hw20emax2.snt.bin"
        pid = getpid()
        hf_main([nuc],sntf,hw,emax;doIMSRG=true,Hsample=1,fn_params="parameters/optional_parameters_forDMD.jl")
        fn_exact = [ "flowOmega/omega_vec_$(pid)$(nuc)_s$(strip(@sprintf("%6.2f",s))).h5" for s in s_pred]
        fns = [ "flowOmega/omega_vec_$(pid)$(nuc)_s$(strip(@sprintf("%6.2f",s))).h5" for s = smin:ds:smax]
        # generating DMD vectors
        trank = 10
        dmd_main(emax, nuc, fns, trank, smin, smax, ds; s_pred=s_pred, fn_exact=fn_exact)
        # restart from DMD
        for s in s_pred
            s_str = strip(@sprintf("%6.2f",s))
            println("s = $s")
            println("DMD:")
            fns_conv = ["flowOmega/omega_dmdvec_e$(emax)_$(nuc)_s$(s_str).h5"]
            hf_main([nuc],sntf,hw,emax;doIMSRG=true,
                    restart_from_files=[fns_conv,[]],                    
                    fn_params="parameters/optional_parameters_forDMD.jl")
        end
    end

end
