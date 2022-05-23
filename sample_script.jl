using NuclearToolkit

function run()
    ### generate NN potential
    make_chiEFTint()

    ### HFMBPT & VS-IMSRG calculation 
    hw = 20; emax=4
    nuc = "O16"; core = "O16"; vspace="sd-shell"
    sntf = "tbme_em500n3lo_srg2.0hw"*string(hw)*"emax"*string(emax)*".snt.bin"
    hf_main([nuc],sntf,hw,emax;verbose=false,doIMSRG=true,corenuc=core,ref="nuc",valencespace=vspace)

    ## shell model calculation
    vs_sntf = "vsimsrg_sd-shell_core"*core*"ref"*nuc*"_hw"*string(hw)*"e"*string(emax)*"_Delta0.0.snt"
    n_eigen=10;targetJ=[]
    main_sm(vs_sntf,"Mg24",n_eigen,targetJ)
    return nothing
end

run()
