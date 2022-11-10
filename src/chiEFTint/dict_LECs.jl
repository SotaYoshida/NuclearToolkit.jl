"""
    prep_valsidxs_dLECs!(vals,idxs,dLECs)
overwrites `vals` and `idxs`, which are to be used in e.g., MCMC algorithms, by dict for LECs `dLECs`.
"""
function prep_valsidxs_dLECs!(vals,idxs,dLECs)
    for (idx,tkey) in enumerate(keys(dLECs))
        tLEC = dLECs[tkey]
        push!(vals,tLEC)
        idxs[tkey] = idx
    end
    return nothing
end

"""
    dict_em500n3lo

returns vector&dict. for EM500N3LO: Entem-Machleidt interaction upto N3LO with 500 MeV cutoff 
"""
function dict_em500n3lo()
    vals = Float64[ ]
    idxs=Dict{String,Int64}()
    dLECs=Dict{String,Float64}() 
    ## LO con  10^4 GeV^-2
    dLECs["C0_1S0"]  = -0.147167
    dLECs["C0_3S1"]  = -0.118972496
    dLECs["C_CSB"] = 0.00049950
    dLECs["C_CIB"] = 0.00069075
    ## NLO con 10^4 GeV^-4
    dLECs["C2_3S1"]  = 0.76
    dLECs["C2_3P0"]  = 1.487
    dLECs["C2_1P1"]  = 0.656
    dLECs["C2_3P1"]  = -0.630
    dLECs["C2_1S0"]  = 2.380
    dLECs["C2_3SD1"] = 0.826
    dLECs["C2_3P2"]  = -0.538
    ## N3LO con 10^4 GeV^-6
    dLECs["hD_1S0"] = -2.545
    dLECs["D_1S0"] = -16.0
    dLECs["D_1P1"] = 5.25
    dLECs["D_3P0"] = 0.245
    dLECs["D_3P1"] = 2.35
    dLECs["D_3P2"] = 2.295
    dLECs["hD_3S1"] = 7.0
    dLECs["D_3S1"] = 6.55
    dLECs["hD_3SD1"] = 2.25
    dLECs["D_3SD1"] = 6.61
    dLECs["D_3D1"] = -2.80
    dLECs["D_1D2"] = -1.77
    dLECs["D_3D2"] = -1.46
    dLECs["D_3PF2"] = -0.465
    dLECs["D_3D3"] = 5.66
    ## For pion exchange contribution
    ## NNLO GeV^-1
    dLECs["c1_NNLO"] = -0.81
    dLECs["c2_NNLO"] =  2.8
    dLECs["c3_NNLO"] = -3.2
    dLECs["c4_NNLO"] =  5.4
    ## for density-dependent NN, c1,c3,c4 can be different from genuine NN part
    dLECs["ct1_NNLO"] = -0.804  
    dLECs["ct3_NNLO"] = -2.828   
    dLECs["ct4_NNLO"] =  4.811   
    dLECs["cD"] =  0.545
    dLECs["cE"] = -0.044 
    ## N3LO GeV^-2
    dLECs["d12"] = 3.06
    dLECs["d3"] = -3.27
    dLECs["d5"] = 0.45
    dLECs["d145"] = -5.65
    ## N4LO (unused)
    dLECs["e14"] = 0.0
    dLECs["e17"] = 0.0
    ### NLOvs  (usually not used)
    dLECs["c_vs_1"] =  0.19162
    dLECs["c_vs_2"] = -0.28374
    dLECs["c_vs_3"] =  0.02680
    dLECs["c_vs_4"] = -0.34499
    dLECs["c_vs_5"] = -0.16235

    prep_valsidxs_dLECs!(vals,idxs,dLECs)
    return LECs(vals,idxs,dLECs)
end

"""
    dict_emn500n3lo

returns vector&dict. for EMN500N3LO: Entem-Machleidt-Nosyk interaction upto N3LO with 500 MeV cutoff 
"""
function dict_emn500n3lo()
    vals = Float64[ ]
    idxs=Dict{String,Int64}()
    dLECs=Dict{String,Float64}() 
    ## LO con  10^4 GeV^-2
    dLECs["C0_1S0"]  = -0.13956293
    dLECs["C0_3S1"]  = -0.159635365
    dLECs["C_CSB"] = 0.00052685
    dLECs["C_CIB"] = 0.00075054
    ## NLO con 10^4 GeV^-4
    dLECs["C2_3S1"]  =  0.823265695
    dLECs["C2_3P0"]  =  1.184990880
    dLECs["C2_1P1"]  =  0.187002997
    dLECs["C2_3P1"]  = -0.819122188
    dLECs["C2_1S0"]  =  2.417201616
    dLECs["C2_3SD1"] =  0.502604767
    dLECs["C2_3P2"]  = -0.757935632
    ## N3LO con 10^4 GeV^-6
    dLECs["hD_1S0"]=  -2.332142565
    dLECs["D_1S0"] = -16.737482775
    dLECs["D_1P1"] =   9.976039933
    dLECs["D_3P0"] =   4.989911864
    dLECs["D_3P1"] =   4.812595281
    dLECs["D_3P2"] =   6.023854940
    dLECs["hD_3S1"]=  -4.319199809
    dLECs["D_3S1"] = -19.171052868
    dLECs["hD_3SD1"]=  1.162495137
    dLECs["D_3SD1"]= 1.759338786
    dLECs["D_3D1"] = -5.59034819
    dLECs["D_1D2"] = -1.946384037
    dLECs["D_3D2"] = -3.200942439
    dLECs["D_3PF2"]= 0.010519022
    dLECs["D_3D3"] = -1.3336602
    ## For pion exchange contribution
    ## NNLO GeV^-1
    dLECs["c1_NNLO"] = -1.07
    dLECs["c2_NNLO"] =  3.2
    dLECs["c3_NNLO"] = -5.32
    dLECs["c4_NNLO"] =  3.56
    ## for density-dependent NN, c1,c3,c4 can be different from genuine NN part
    dLECs["ct1_NNLO"] = dLECs["c1_NNLO"]
    dLECs["ct3_NNLO"] = dLECs["c3_NNLO"]
    dLECs["ct4_NNLO"] = dLECs["c4_NNLO"]
    dLECs["cD"] = 0.0
    dLECs["cE"] = 0.0

    ## N3LO GeV^-2
    dLECs["d12"] = 1.04
    dLECs["d3"] =  -0.48
    dLECs["d5"] = 0.14
    dLECs["d145"] = -1.90
    ## N4LO GeV^-3
    dLECs["e14"] = 0.0
    dLECs["e17"] = 0.0
    ### NLOvs  (usually not used)
    dLECs["c_vs_1"] = 0.0
    dLECs["c_vs_2"] = 0.0
    dLECs["c_vs_3"] = 0.0
    dLECs["c_vs_4"] = 0.0
    dLECs["c_vs_5"] = 0.0
    prep_valsidxs_dLECs!(vals,idxs,dLECs)
    return LECs(vals,idxs,dLECs)
end


"""
    dict_emn500n4lo()

returns vector&dict. for EMN500N4LO: Entem-Machleidt-Nosyk interaction upto N4LO with 500 MeV cutoff 

"""
function dict_emn500n4lo()
    vals = Float64[ ]
    idxs=Dict{String,Int64}()
    dLECs=Dict{String,Float64}() 
    ## LO con  10^4 GeV^-2
    dLECs["C0_1S0"]  = -0.142358500
    dLECs["C0_3S1"]  = -0.151581533
    dLECs["C_CSB"] = 0.00047905
    dLECs["C_CIB"] = 0.000682975
    ## NLO con 10^4 GeV^-4
    dLECs["C2_3S1"]  = 0.867156546
    dLECs["C2_3P0"]  = 1.214966189
    dLECs["C2_1P1"]  = 0.232092261
    dLECs["C2_3P1"]  = -0.859718490
    dLECs["C2_1S0"]  = 2.494633548
    dLECs["C2_3SD1"] = 0.550117472
    dLECs["C2_3P2"]  = -0.763178505
    ## N3LO con 10^4 GeV^-6
    dLECs["hD_1S0"] = -1.779896956
    dLECs["D_1S0"] = -13.385692976
    dLECs["D_1P1"] =  10.956500634
    dLECs["D_3P0"] =   4.571860915
    dLECs["D_3P1"] =   4.133947333
    dLECs["D_3P2"] =   5.342585336
    dLECs["hD_3S1"] = -2.949089421
    dLECs["D_3S1"] = -20.793199632
    dLECs["hD_3SD1"] = 1.3545478412
    dLECs["D_3SD1"] = 2.176852098
    dLECs["D_3D1"] = -6.01826561
    dLECs["D_1D2"] = -1.545851484
    dLECs["D_3D2"] = -3.383001593
    dLECs["D_3PF2"]=  0.076818254
    dLECs["D_3D3"] = -0.664832879
    ## N4LO con 10^4 GeV^-8
    dLECs["E_3F2"] =  1.222149783
    dLECs["E_1F3"] =  1.316452134
    dLECs["E_3F4"] = -0.321096808
    ## For pion exchange contribution
    ## NNLO GeV^-1
    dLECs["c1_NNLO"] = -1.10
    dLECs["c2_NNLO"] =  3.57
    dLECs["c3_NNLO"] = -5.54
    dLECs["c4_NNLO"] =  4.17
    # for density-dependent NN, c1,c3,c4 can be different from genuine NN part
    ## EMN effective c's 
    dLECs["ct1_NNLO"] = -0.73
    dLECs["ct3_NNLO"] = -3.38
    dLECs["ct4_NNLO"] = 1.69
    ## NN sector
    dLECs["ct1_NNLO"] = dLECs["c1_NNLO"]
    dLECs["ct3_NNLO"] = dLECs["c3_NNLO"]
    dLECs["ct4_NNLO"] = dLECs["c4_NNLO"]

    dLECs["cD"] =  0.0
    dLECs["cE"] =  0.0

    # ## old MAP
    # dLECs["ct1_NNLO"] = -1.10
    # dLECs["ct3_NNLO"] = -4.753
    # dLECs["ct4_NNLO"] = 4.089
    # dLECs["cD"] = -0.277
    # dLECs["cE"] = 0.820

    ## N3LO GeV^-2
    dLECs["d12"] = 6.18
    dLECs["d3"] = -8.91
    dLECs["d5"] =  0.86
    dLECs["d145"] = -12.18
    ## N4LO GeV^-3
    dLECs["e14"] = 1.18
    dLECs["e17"] = -0.18
    ### NLOvs  (usually not used)
    dLECs["c_vs_1"] = 0.0
    dLECs["c_vs_2"] = 0.0
    dLECs["c_vs_3"] = 0.0
    dLECs["c_vs_4"] = 0.0
    dLECs["c_vs_5"] = 0.0
    prep_valsidxs_dLECs!(vals,idxs,dLECs)
    return LECs(vals,idxs,dLECs)
end
