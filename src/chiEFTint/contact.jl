# ** Leading Order (LO)** Q^0 contact term
function LO(chiEFTobj,xr,LECs,V12mom,dict_numst,to)
    C_1S0 = LECs["C0_1S0"]
    C_3S1 = LECs["C0_3S1"]
    C_CIB = LECs["C_CIB"]
    C_CSB = LECs["C_CSB"]
    fac_q = (hc/(2*pi))^3 * 1.e-2
    pfunc = f1
    n_reg = 3    
    ress = [2*C_CIB+C_CSB,0.0,2*C_CIB-C_CSB]
    for pnrank = 1:3
        TF = true
        res = ress[pnrank]
        ##1S0 (LO)
        l=0;lp=0;S=0;J=0
        LEC = (C_1S0+res) * fac_q
        calc_Vmom!(chiEFTobj,pnrank,V12mom,dict_numst[pnrank],xr,LEC,LEC,l,lp,S,J,pfunc,n_reg,to)
        
        ## 3S1 (LO)
        l=0;lp=0;S=1;J=1
        if pnrank%2==1;continue;end
        LEC = (C_3S1) * fac_q 
        calc_Vmom!(chiEFTobj,pnrank,V12mom,dict_numst[pnrank],xr,LEC,LEC,l,lp,S,J,pfunc,n_reg,to)
    end
    return nothing
end

"""
    NLO(chiEFTobj,xr,LECs,V12mom,dict_numst,to;n_regulator=2)

Calculate Next to Leading Order (NLO), ``Q^2`` contact term.
LECs=>C2_3S1,C2_3P0,C2_1P1,C2_3P1,C2_1S0,C2_3SD1,C2_3DS1,C2_3P2

The power for exponential regulator `n_reg` is 2 (default),
For 3P1 channel, `n_reg=3` is used as in the R.Machleidt's fortran code.
I don't know the reason, but LECs may have been determined with `n_reg=3` potential.
"""
function NLO(chiEFTobj,xr,LECs,V12mom,dict_numst,to;n_regulator=2)
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
            if occursin("emn500",chiEFTobj.pottype) && n==4
                n_reg = 3
            else
                n_reg = n_regulator
            end
            calc_Vmom!(chiEFTobj,pnrank,V12mom,dict_numst[pnrank],xr,LEC,LEC,l,lp,S,J,pfunc,n_reg,to)
        end
    end
    return nothing
end

function N3LO(chiEFTobj,xr,LECs,V12mom,dict_numst,to)
    fac = hc^3  *  1.e-14  / (2*pi)^3

    LECs_N3LO = [[LECs["hD_1S0"],LECs["D_1S0"]],
                 [LECs["D_3P0"]],
                 [LECs["D_1P1"]],
                 [LECs["D_3P1"]],
                 [LECs["hD_3S1"],LECs["D_3S1"]],
                 [LECs["D_3D1"]],
                 [LECs["hD_3SD1"],LECs["D_3SD1"]],
                 [LECs["hD_3SD1"],LECs["D_3SD1"]],
                 [LECs["D_1D2"]],
                 [LECs["D_3D2"]],
                 [LECs["D_3P2"]],
                 [LECs["D_3PF2"]],
                 [LECs["D_3PF2"]],
                 [LECs["D_3D3"]]]
    #hD_1S0,D_1S0,D_3P0,D_1P1,D_3P1,hD_3S1,D_3S1,D_3D1,hD_3SD1,D_3SD1,hD_3SD1,D_3SD1,D_1D2,D_3D2,D_3P2,D_3PF2,D_3PF,D_3D3
    #LECs_N3LO [[-2.545, -16.0], [0.245], [5.25], [2.35], [7.0, 6.55], [-2.8], [2.25, 6.61], [2.25, 6.61], [-1.77], [-1.46], [2.295], [-0.465], [-0.465], [5.66]]

    funcs = [ f_442,f_31,f_31,f_31,
              f_442,f_x2y2,f_x42,f_y42,
              f_x2y2,f_x2y2,f_31,f_x3y,f_xy3,f_x2y2]
    llpSJ_s = [ [0,0,0,0],[1,1,1,0],[1,1,0,1],[1,1,1,1],
                [0,0,1,1],[2,2,1,1],
                [0,2,1,1],[2,0,1,1],
                [2,2,0,2],[2,2,1,2],[1,1,1,2],[1,3,1,2],[3,1,1,2],[2,2,1,3]]
    nregs = [2,3,2,4, 2,2,2,2, 4,2,2,4,4,-1]
    if occursin("emn500",chiEFTobj.pottype)
        nregs[4] = nregs[9] = nregs[10] = 3
        nregs[14] = 2
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
            calc_Vmom!(chiEFTobj,pnrank,V12mom,dict_numst[pnrank],xr,LEC,LEC2,l,lp,S,J,pfunc,n_reg,to)
        end
    end
    return nothing
end

