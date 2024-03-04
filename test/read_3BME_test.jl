@testset "Reading 3BME file from NuHamil code (by T. Miyagi)" begin
    sntf   = "tbme_em500n3lo_barehw20emax2.snt.bin"
    fn_3nf = "interaction_file/small_3BME_N3LOlnl_4_8_4.me3j.gz"
    hw = 16
    emax = 2; e3max = 4; e1max_file=4; e2max_file=8; e3max_file=4
    HFobj = hf_main(["He6"],sntf,hw,emax;return_obj=true,
                    e1max_file=e1max_file,e2max_file=e2max_file,e3max_file=e3max_file,e3max=e3max,fn_3nf=fn_3nf)
    Es = [HFobj.E0, HFobj.EMP2, HFobj.EMP3]
    Eref = [ 8.160870, -9.58767, -0.82059]
    @test (Es[1] - Eref[1])^2 + (Es[2] - Eref[2])^2 + (Es[3] - Eref[3])^2 < 1.e-7
end
