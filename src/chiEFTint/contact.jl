"""
    LO(chiEFTobj,to)    

Calculate Leading Order (LO), `Q^0`` contact term.
"""
function LO(chiEFTobj,to)    
    xr = chiEFTobj.xr
    LECs = chiEFTobj.LECs.dLECs
    V12mom = chiEFTobj.V12mom
    dict_pwch = chiEFTobj.dict_pwch
    C_1S0 = LECs["C0_1S0"]
    C_3S1 = LECs["C0_3S1"]
    C_CIB = LECs["C_CIB"]
    C_CSB = LECs["C_CSB"]
    fac_q = (hc/(2*pi))^3 * 1.e-2
    pfunc = f1
    n_reg = 3    
    ress = [2*C_CIB+C_CSB,0.0,2*C_CIB-C_CSB]
    for pnrank = 1:3
        res = ress[pnrank]
        ##1S0 (LO)
        l=0;lp=0;S=0;J=0
        LEC = (C_1S0+res) * fac_q
        calc_Vmom!(chiEFTobj.params,pnrank,V12mom,dict_pwch[pnrank],xr,LEC,LEC,l,lp,S,J,pfunc,n_reg,to)
        ## 3S1 (LO)
        l=0;lp=0;S=1;J=1
        if pnrank%2==1;continue;end
        LEC = (C_3S1) * fac_q 
        calc_Vmom!(chiEFTobj.params,pnrank,V12mom,dict_pwch[pnrank],xr,LEC,LEC,l,lp,S,J,pfunc,n_reg,to)
    end
    return nothing
end

"""
    NLO(chiEFTobj,to;n_regulator=2)

Calculate Next to Leading Order (NLO), ``Q^2`` contact term.
LECs=>C2_3S1,C2_3P0,C2_1P1,C2_3P1,C2_1S0,C2_3SD1,C2_3DS1,C2_3P2

The power for exponential regulator `n_reg` is 2 (default),
For 3P1 channel, `n_reg=3` is used as in the R.Machleidt's fortran code.
I don't know the reason, but the LECs in EMN potential) may have been determined with `n_reg=3` potential.
"""
function NLO(chiEFTobj,to;n_regulator=2)
    xr = chiEFTobj.xr
    LECs = chiEFTobj.LECs.dLECs
    V12mom = chiEFTobj.V12mom
    dict_pwch = chiEFTobj.dict_pwch
    fac = hc^3  *  1.e-8  / (2*pi)^3  
    tLECs = [ LECs["C2_3S1"],LECs["C2_3P0"],LECs["C2_1P1"],
              LECs["C2_3P1"],LECs["C2_1S0"],
              LECs["C2_3SD1"],LECs["C2_3SD1"],LECs["C2_3P2"]]
    funcs = [fxxpyy,fxy,fxy,fxy,fxxpyy,fx2,fy2,fxy] #fx2
    llpSJ_s = [ [0,0,1,1],[1,1,1,0],[1,1,0,1],[1,1,1,1],
                [0,0,0,0],[0,2,1,1],[2,0,1,1],[1,1,1,2]] #
    for pnrank = 1:3
        for n=1:8            
            LEC = tLECs[n] * fac
            pfunc = funcs[n]
            l,lp,S,J = llpSJ_s[n]
            if pnrank%2 == 1 && (l+S+1) % 2 != 1;continue;end
            if occursin("emn500",chiEFTobj.params.pottype) && n==4
                n_reg = 3
            else
                n_reg = n_regulator
            end
            calc_Vmom!(chiEFTobj.params,pnrank,V12mom,dict_pwch[pnrank],xr,LEC,LEC,l,lp,S,J,pfunc,n_reg,to)
        end
    end
    return nothing
end

"""
    N3LO(chiEFTobj,to)

Calculate Next-to-next-to-next-to Leading Order (N3LO), ``Q^4`` contact term.
"""
function N3LO(chiEFTobj,to) 
    xr = chiEFTobj.xr
    LECs = chiEFTobj.LECs.dLECs
    V12mom = chiEFTobj.V12mom
    dict_pwch = chiEFTobj.dict_pwch
    fac = hc^3  *  1.e-14  / (2*pi)^3
    LECs_N3LO = [[LECs["hD_1S0"],LECs["D_1S0"]],[LECs["D_3P0"]],
                 [LECs["D_1P1"]],[LECs["D_3P1"]],
                 [LECs["hD_3S1"],LECs["D_3S1"]],[LECs["D_3D1"]],
                 [LECs["hD_3SD1"],LECs["D_3SD1"]],[LECs["hD_3SD1"],LECs["D_3SD1"]],
                 [LECs["D_1D2"]],[LECs["D_3D2"]],[LECs["D_3P2"]],
                 [LECs["D_3PF2"]],[LECs["D_3PF2"]],[LECs["D_3D3"]]]
    funcs = [ f_442,f_31,f_31,f_31,
              f_442,f_x2y2,f_x42,f_y42,
              f_x2y2,f_x2y2,f_31,f_x3y,f_xy3,f_x2y2]
    llpSJ_s = [ [0,0,0,0],[1,1,1,0],[1,1,0,1],[1,1,1,1],
                [0,0,1,1],[2,2,1,1],
                [0,2,1,1],[2,0,1,1],
                [2,2,0,2],[2,2,1,2],[1,1,1,2],[1,3,1,2],[3,1,1,2],[2,2,1,3]]
    nregs = [2,3,2,4, 2,2,2,2, 4,2,2,4,4,-1]
    if occursin("emn500",chiEFTobj.params.pottype) # power in the regulator 
        nregs[4] = nregs[9] = nregs[10] = 3; nregs[14] = 2
    end
    for pnrank = 1:3
        for (n,tLECs) in enumerate(LECs_N3LO)
            l,lp,S,J = llpSJ_s[n]            
            if pnrank%2 == 1 && (l+S+1) % 2 != 1;continue;end
            LEC = 0.0; LEC2 = 0.0
            if length(tLECs)== 1;
                LEC =LEC2 = tLECs[1]
            else
                LEC = tLECs[1]; LEC2 = tLECs[2]
            end                
            pfunc = funcs[n]
            n_reg = nregs[n]
            LEC *= fac; LEC2 *= fac
            calc_Vmom!(chiEFTobj.params,pnrank,V12mom,dict_pwch[pnrank],xr,LEC,LEC2,l,lp,S,J,pfunc,n_reg,to)
        end
    end
    return nothing
end

"""
    N4LO(chiEFTobj,to;n_reg=2)

Calculate Next-to-next-to-next-to-next-to Leading Order (N4LO) contact term.
"""
function N4LO(chiEFTobj,to;n_reg=2) 
    xr = chiEFTobj.xr
    LECs = chiEFTobj.LECs.dLECs
    V12mom = chiEFTobj.V12mom
    dict_pwch = chiEFTobj.dict_pwch
    fac = hc^3  *  1.e-20  / (2*pi)^3
    LECs_N4LO = [LECs["E_3F2"],LECs["E_1F3"],LECs["E_3F4"]]
    llpSJ_s = [[3,3,1,2],[3,3,0,3],[3,3,1,4]]
    pfunc = f_x3y3
    for pnrank = 1:3
        for (n,LEC) in enumerate(LECs_N4LO)
            l,lp,S,J = llpSJ_s[n]            
            if pnrank%2 == 1 && (l+S+1) % 2 != 1;continue;end
            LEC *= fac; LEC2 = LEC
            calc_Vmom!(chiEFTobj.params,pnrank,V12mom,dict_pwch[pnrank],xr,LEC,LEC2,l,lp,S,J,pfunc,n_reg,to)
        end
    end
    return nothing
end