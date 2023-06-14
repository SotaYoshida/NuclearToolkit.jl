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

function nn_IMSRG_sampling(nucs;is_show=false,itnum=1,corenuc="",ref="nucl",fn_params="ann_sample_params.jl",valencespace="")
    to = TimerOutput()    
    optimizer = ""
    MPIcomm = false
    io = stdout
    chiEFTobj,OPTobj,dWS = construct_chiEFTobj(false,itnum,optimizer,MPIcomm,io,to;fn_params)
    org_dLECs = copy(chiEFTobj.LECs.dLECs)
    
    for it = 1:itnum
        calcualte_NNpot_in_momentumspace(chiEFTobj,to)
        V12mom = chiEFTobj.V12mom
        V12mom_2n3n = chiEFTobj.V12mom_2n3n
        SRG(chiEFTobj,to)
        calc_vmom_3nf(chiEFTobj,it,to)
        TMtrans(chiEFTobj,dWS,to)
        sntf = chiEFTobj.params.fn_tbme
        hw = chiEFTobj.params.hw
        emax = chiEFTobj.params.emax
        if valencespace == ""
            hf_main(nucs,sntf,hw,emax;is_show=is_show,doIMSRG=true,corenuc="",ref="nuc",fn_params=fn_params,Hsample=true,Operators=["Rp2"])
        else
            hf_main(nucs,sntf,hw,emax;is_show=is_show,doIMSRG=true,valencespace=valencespace,corenuc=corenuc,ref=ref,fn_params=fn_params,Hsample=true,Operators=["Rp2"])
        end
        for ch = 1:length(V12mom)
            V12mom[ch] .= 0.0
            V12mom_2n3n[ch] .= 0.0
        end
        updateLECs!(chiEFTobj,org_dLECs)

     end
    show_TimerOutput_results(to;tf=is_show)
end
