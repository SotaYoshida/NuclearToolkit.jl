function updateLECs!(chiEFTobj,org_dLECs)
    LECvals = chiEFTobj.LECs.vals
    idxLECs = chiEFTobj.LECs.idxs
    dLECs = chiEFTobj.LECs.dLECs
    for tkey in keys(dLECs)
        idx = idxLECs[tkey] 
        oLEC = org_dLECs[tkey]
        nLEC = oLEC + 0.01 * randn()
        dLECs[tkey] = nLEC
        LECvals[idx] = nLEC
    end 
    return nothing
end

function nn_IMSRG_sampling(;is_show=false,itnum=1,nucs=[],corenuc="",ref="nucl",fn_params="sample_params.jl")
    to = TimerOutput()    
    optimizer = ""
    MPIcomm = false
    io = stdout
    chiEFTobj,OPTobj,d9j,HOBs = construct_chiEFTobj(false,itnum,optimizer,MPIcomm,io,to;fn_params)
    org_dLECs = copy(chiEFTobj.LECs.dLECs)
    nuc = "He4"
    for it = 1:itnum
        calcualte_NNpot_in_momentumspace(chiEFTobj,to)
        V12mom = chiEFTobj.V12mom
        V12mom_2n3n = chiEFTobj.V12mom_2n3n
        SRG(chiEFTobj,to)
        TMtrans(chiEFTobj,HOBs,to)
        sntf = chiEFTobj.params.fn_tbme
        hw = chiEFTobj.params.hw
        emax = chiEFTobj.params.emax
        hf_main([nuc],sntf,hw,emax;is_show=is_show,doIMSRG=true,corenuc="",ref="nuc",fn_params=fn_params,debugmode=2)
        for ch = 1:length(V12mom)
            V12mom[ch] .= 0.0; V12mom_2n3n[ch] .= 0.0
        end
        updateLECs!(chiEFTobj,org_dLECs)

        pid = getpid()

        fname = "flowOmega/HF_$pid"*nuc*"_1.bin"
        run(`mv $fname flowOmega/HF_$(it).snt.bin`)
        fname = "flowOmega/Hs_$pid"*nuc*"_1.bin"
        run(`mv $fname flowOmega/Hs_$(it).snt.bin`)
    end
    show_TimerOutput_results(to;tf=is_show)
end