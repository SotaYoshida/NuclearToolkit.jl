const K = 6144.0
const lam_gav = 1.2701

"""
    eval_betadecay_from_kshell_log(fns_sum_parent::Vector{String},fns_sum_daughter::Vector{String},fns_GT::Vector{String},fns_FF::Vector{String},parentJpi::String;pnuc="",dnuc="",verbose=false)

# Arguments
- `fns_sum_parent::Vector{String}` vector of paths to KSHELL-summary files for the parent nucleus
- `fns_sum_daughter::Vector{String}` vector of paths to KSHELL-summary files for the daughter nucleus
- `fns_GT::Vector{String}` vector of paths to KSHELL-log files for Gamow-Teller transitions.
- `fns_FF::Vector{String}` vector of paths to KSHELL-log files for First-Forbidden transitions.

# Optional Arguments
- `pnuc::String` to specify the parent nuclei explicitly
- `dnuc::String` to specify the daughter nuclei explicitly
- `verbose::Bool` option to show GT/FF transitions for each state
"""
function eval_betadecay_from_kshell_log(fns_sum_parent::Vector{String},fns_sum_daughter::Vector{String},fns_GT::Vector{String},fns_FF::Vector{String},parentJpi::String;pnuc="",dnuc="",verbose=false)
    qfactors = get_quenching_factors()
    parent = read_kshell_summary(fns_sum_parent;nuc=pnuc,targetJpi=parentJpi)
    daughter = read_kshell_summary(fns_sum_daughter;nuc=dnuc)
    qvals = get_qval(parent,daughter,parentJpi)
    eval_logft_f(parent,daughter,qvals,fns_GT,fns_FF,qfactors;verbose=verbose)
end

mutable struct quenching_factors
    qGT::Float64
    qMS0::Float64
    qMT0::Float64
    qx::Float64
    qu::Float64
    qz::Float64
end

"""
    Ecoulomb_SM(Z,N)

Additional Coulomb energy (MeV)
```math
E_C(Z,N) = 0.7 \\frac{Z(Z-1) -0.76 [Z(Z-1)]^{2/3}}{R_c} \\\\
R_c = e^{1.5/A} A^{1/3} \\left( 0.946 -0.573 \\left( \\frac{2T}{A} \\right)^2  \\right) \\\\
T = |Z-N|/2 
```
used in references.
* S.Yoshida et al., Phys. Rev. C 97, 054321 (2018).
* J. Duflo and A. P. Zuker, Phys. Rev. C 52, R23 (1995).
* E. Caurier et al., Phys. Rev. C 59, 2033 (1999).
"""
function Ecoulomb_SM(Z,N)
    A = Z+N
    T = abs(Z-N)/2
    Rc = exp(1.5/A) * A^(1/3) * (0.946-0.573* (2*T/A)^2)
    return 0.7 * (Z*(Z-1) -0.76 * (Z*(Z-1))^(2/3)) / Rc
end

function get_quenching_factors(qtype="")
    if qtype == "" || qtype == "SY" 
        # S. Yoshida et al., PRC 97, 054321(2018).
        return quenching_factors(0.74, 1.0, 1.7, 1.0, 1.0, 0.55)
    elseif qtype == "W" || qtype == "Warburton" 
        # E. Warburton, J. Becker, B. Brown, and D. Millener, Ann. Phys. 187, 471 (1988).
        return quenching_factors(1.0, 1.1, 1.5, 1.0, 1.0, 0.51)
    else
        @warn "qtype=$qtype is not supported! SY value will be used"
        return quenching_factors(0.74, 1.0, 1.7, 1.0, 1.0, 0.55) 
    end
end

function get_qval(parent::kshell_nuc,daughter::kshell_nuc,parentJpi::String)
    pZ = parent.nuc.Z; dZ = daughter.nuc.Z
    pN = parent.nuc.N; dN = daughter.nuc.N
    deltam = (Mn - Mp - Me)
    # get theoretical Qval
    qval_tho = parent.Egs_target + Ecoulomb_SM(pZ,pN) - (daughter.Egs + Ecoulomb_SM(dZ,dN)) + deltam
    # get experimental Qval
    qval_exp = 0.0 
    qvalTF= ""
    try
        qval_exp = ame2020data[parent.nuc.cnuc].Qval
        qvalTF = ifelse(occursin("#",qval_exp),"#","")
        qval_exp = parse(Float64,replace(qval_exp,"#"=>"")) / 1000.0
    catch
        println("experimental qval is not available from amedata.jl, so Qval is set 0.0")
    end
    parEx = parent.Egs_target - parent.Egs
    gsTF = ifelse(parEx > 0.0,",Ex.="*strip(@sprintf("%9.2f",parEx))*" MeV","")
    println(parent.nuc.cnuc,"($(parentJpi)$(gsTF))=>",daughter.nuc.cnuc, " qval_exp$(qvalTF) $qval_exp qval_tho $qval_tho")    
    return [qval_exp,qval_tho]
end

function eval_logft_f(parent::kshell_nuc,daughter::kshell_nuc,Qvals,fns_GT::Vector{String},fns_FF::Vector{String},qfactors::quenching_factors;verbose=false)
    GTdata = eval_bgt_files(fns_GT,qfactors,parent,daughter,Qvals)
    FFkeys,FFdata = eval_bFF_files(fns_FF,qfactors,parent,daughter,Qvals;verbose=verbose)
    calc_halflife(daughter,GTdata,FFkeys,FFdata;verbose=verbose)
end

function radius_nuc_for_FermiIntegarl(A;formula=1)
    if formula == 1
        return (1.15 + 1.8*A^(-2/3) -1.2*A^(-4/3)) * A^(1/3)
    elseif formula == 2
        return 1.123A^(1/3) - 0.941A^(-1/3)
    end
end

function calc_halflife(daughter::kshell_nuc,GTdata,FFkeys,FFdata;verbose=false)
    hl_GT_expQ = hl_FF_expQ = hl_total_expQ = 0.0
    hl_GT_thoQ = hl_FF_thoQ = hl_total_thoQ = 0.0
    Qw_hl_GT_expQ = Qw_hl_FF_expQ = Qw_hl_total_expQ = 0.0
    Qw_hl_GT_thoQ = Qw_hl_FF_thoQ = Qw_hl_total_thoQ = 0.0 
    Sn = NaN
    try
        Sn = parse(Float64,replace(ame2020data[daughter.nuc.cnuc].Sn,"#"=>"")) / 1000.0
    catch   
        @warn "experimental Sn for $(daughter.nuc.cnuc) is not available from amedata.jl"
    end
    for GTstate in GTdata
        if verbose #&& (GTstate.hl_expQ != Inf)
            println("Energy ", @sprintf("%12.4f",GTstate.Energy), " J2prty",@sprintf("%4i",GTstate.J2)*GTstate.prty,
                    " Ex. ", @sprintf("%9.3f",GTstate.Ex),
                    " BGT ", @sprintf("%13.5f", GTstate.BGT),
                    " log10ft(GT) ", @sprintf("%9.2f", GTstate.logft),
                    " hl(expQ) ", @sprintf("%12.4e",GTstate.hl_expQ),
                    " hl(thoQ) ", @sprintf("%12.4e",GTstate.hl_thoQ))
        end
        if GTstate.hl_expQ != 0.0
            hl_GT_expQ += 1.0 / GTstate.hl_expQ
            if GTstate.Ex >= Sn
                Qw_hl_GT_expQ += 1.0 / GTstate.hl_expQ
            end
        end
        if GTstate.hl_thoQ != 0.0
            hl_GT_thoQ += 1.0 / GTstate.hl_thoQ
            if GTstate.Ex >= Sn
                Qw_hl_GT_thoQ += 1.0 / GTstate.hl_thoQ
            end
        end
    end
    for tkey in FFkeys
        FFstate = FFdata[tkey]        
        if FFstate.hl_expQ != 0.0
            hl_FF_expQ += 1.0 / FFstate.hl_expQ
            if FFstate.Ex >= Sn
               Qw_hl_FF_expQ += 1.0 / FFstate.hl_expQ
            end
        end
        if FFstate.hl_thoQ != 0.0
            hl_FF_thoQ += 1.0 / FFstate.hl_thoQ
            if FFstate.Ex >= Sn
                Qw_hl_FF_thoQ += 1.0 / FFstate.hl_thoQ
            end
        end
        if verbose && (FFstate.f0_thoQ != Inf)
            logft_expQ = log10(FFstate.f0_expQ*FFstate.hl_expQ)
            logft_thoQ = log10(FFstate.f0_thoQ*FFstate.hl_thoQ)
            println("Energy ", @sprintf("%12.4f",tkey), " J2prty",@sprintf("%4i",FFstate.J2)*FFstate.prty, " Ex. ", @sprintf("%9.3f",FFstate.Ex),
                    " log10ft(expQ) ", @sprintf("%9.2f", logft_expQ)," log10ft(thoQ) ", @sprintf("%9.2f", logft_thoQ),
                    " f0 ", @sprintf("%9.2e",FFstate.f0_expQ),@sprintf("%9.2e",FFstate.f0_thoQ))
        end
    end
    
    hl_total_expQ = 1.0 / (ifelse(hl_GT_expQ==Inf,0.0,hl_GT_expQ) + ifelse(hl_FF_expQ==Inf,0.0,hl_FF_expQ))
    hl_total_thoQ = 1.0 / (ifelse(hl_GT_thoQ==Inf,0.0,hl_GT_thoQ) + ifelse(hl_FF_thoQ==Inf,0.0,hl_FF_thoQ))
    Qw_hl_total_expQ = 1.0 / (ifelse(Qw_hl_GT_expQ==Inf,0.0,Qw_hl_GT_expQ) + ifelse(Qw_hl_FF_expQ==Inf,0.0,Qw_hl_FF_expQ))
    Qw_hl_total_thoQ = 1.0 / (ifelse(Qw_hl_GT_thoQ==Inf,0.0,Qw_hl_GT_thoQ) + ifelse(Qw_hl_FF_thoQ==Inf,0.0,Qw_hl_FF_thoQ))
    Pn_GT_expQ = ifelse(isnan(Sn),NaN,100 * Qw_hl_GT_expQ / hl_GT_expQ)
    Pn_GT_thoQ = ifelse(isnan(Sn),NaN,100 * Qw_hl_GT_thoQ / hl_GT_thoQ)
    Pn_full_expQ = ifelse(isnan(Sn),NaN,100 / (Qw_hl_total_expQ / hl_total_expQ))
    Pn_full_thoQ = ifelse(isnan(Sn),NaN,100 / (Qw_hl_total_thoQ / hl_total_thoQ))
    hl_GT_expQ = 1.0 / hl_GT_expQ; hl_FF_expQ = 1.0 / hl_FF_expQ
    hl_GT_thoQ = 1.0 / hl_GT_thoQ; hl_FF_thoQ = 1.0 / hl_FF_thoQ
    approp_hl_total_expQ = get_appropunit(hl_total_expQ)
    approp_hl_total_thoQ = get_appropunit(hl_total_thoQ)
    println("ExpQ: "," hl(GT-only) ",@sprintf("%12.4e",hl_GT_expQ), get_appropunit(hl_GT_expQ),
            " hl(GT+FF)[sec] ",@sprintf("%12.4e",hl_total_expQ), approp_hl_total_expQ,
            " Pn_GT", @sprintf("%9.2f",Pn_GT_expQ)," Pn_tot", @sprintf("%9.2f",Pn_full_expQ) )
    println("ThoQ: "," hl(GT-only) ",@sprintf("%12.4e",hl_GT_thoQ), get_appropunit(hl_GT_thoQ),
            " hl(GT+FF)[sec] ",@sprintf("%12.4e",hl_total_thoQ), approp_hl_total_thoQ,
            " Pn_GT", @sprintf("%9.2f",Pn_GT_thoQ)," Pn_tot", @sprintf("%9.2f",Pn_full_thoQ))
    return nothing
end

function get_appropunit(hl_in_sec)
    hl = hl_in_sec
    degit = log10(hl_in_sec)
    aunit = "sec"
    if hl >= 0.8 * 31556926
        aunit = "y  "; hl /= 31556926
    elseif 60 <= hl #< 2 *3600
        aunit = "min"; hl /= 60.0
    elseif -3 <= degit < 0
        aunit = "ms "; hl /= 10^(-3)
    elseif -6 <= degit < -3
        aunit = "μs "; hl /= 10^(-6)
    end
    return  @sprintf("%12.3f",hl) * " " * aunit 
end

"""
    Fermi_function(Z,We,R)

```math 
F(Z,W) = 4 (2p_e R)^{-2(1-\\gamma_1)} \\exp \\left( \\pi y \\frac{|\\Gamma(\\gamma_1 + iy)|^2}{ |\\Gamma(2\\gamma_1 + 1)|^2}  \\right) \\\\
\\gamma_k = \\sqrt{k^2 - \\alpha^2 Z^2} \\\\
y = \\alpha Z W / p_e
```
"""
function Fermi_function(Z,We,R)
    pe =  sqrt(We^2 - 1.0)     
    gamma1 = sqrt(1.0 - (fsalpha *Z)^2)  # invalid for Z > 137
    y = fsalpha * Z * We / pe
    gam_nume = abs(gamma(gamma1 + y*im))^2
    gam_deno = abs(gamma(2*gamma1 + 1.0))^2
    gamfac = gam_nume / gam_deno
    return 4 * (2*pe*R)^(-2*(1-gamma1)) * exp( pi*y) * gamfac
end

"""
    Fermi_integral(Qval,Ex,pZ,A,nmesh=40)

Fermi integrals are evaluated with Eqs. in Chapter 5 of "Weak Interactions and Nuclear Beta Decay" by H. F. Schopper.]
```math
f_0 = \\int^{W_0}_{1} F(Z,W) \\sqrt{W^2-1} W(W_0-W)^2 dW, \\\\
W_0 = Q(\\beta^-) + 1 - E_\\mathrm{ex.}
```
Note that Ws are in natural unit, divided by ``m_ec^2``.
"""
function Fermi_integral(qval,Z,R,p=0;nmesh=40)
    W0 = qval / Me + 1.0
    if qval < 0.0; return 0.0; end
    f0 = 0.0
    xs,ws = Gauss_Legendre(1.0,W0,nmesh)
    for i in eachindex(xs)
        We = xs[i]
        weight = ws[i]
        f0 +=  weight * We^(p+1) * (W0-We)^2 * sqrt(We^2-1) * Fermi_function(Z,We,R)
    end
    return f0
end

##################
## Gamow-Teller ##
##################
struct bgtdata
    n::Int
    Energy::Float64
    Ex::Float64
    J2::Int
    prty::String
    S::Float64
    BGT::Float64
    logft::Float64
    f0_expQ::Float64
    f0_thoQ::Float64
    hl_expQ::Float64
    hl_thoQ::Float64
end

function eval_bgt_files(fns,qfactors::quenching_factors,parent,daughter,Qvals;using_strength_funciton_method=true)
    data = bgtdata[ ]
    for fn in fns
        if using_strength_funciton_method 
            read_bgtstrength_file!(fn,qfactors,parent,daughter,Qvals,data)
        else
            @error "only the log files using strength function method is supported now"
            #read_bgt_file!(fn,qfactors,parent,daughter,Qvals,data)
        end
    end
    return data
end

"""

"""
function read_bgtstrength_file!(fn::String,qfactors::quenching_factors,parent::kshell_nuc,daughter::kshell_nuc,Qvals,data)
    dZ = parent.nuc.Z
    A = parent.nuc.A
    R = radius_nuc_for_FermiIntegarl(A)
    gAV = - lam_gav
    qGT = qfactors.qGT
    f = open(fn,"r")
    lines = readlines(f) 
    close(f)
    maxit = count = n_eigen = 0
    J2 = -1
    prty = ""
    for line in lines
        if occursin("MTOT",line)        
            J2 = parse(Int, replace(split(split(line, "MTOT")[end],",")[1], "="=>""))
            #J2 = parse(Int,split(split(line)[end],",")[1])
        end
        if occursin("parity =",line)
            prty = replace(split(split(line, "parity")[end],",")[1], "="=>"")
            #prty = strip(split(line,"parity =")[end])
        end
        if occursin("N_EIGEN",line)
            n_eigen = parse(Int, replace(split(split(line, "N_EIGEN")[end],",")[1], "="=>""))
            #n_eigen = parse(Int,split(split(line,"=")[end],",")[1])
        end
        if !occursin("strength function ",line); continue; end
        tl = split(line)
        itnum = parse(Int,tl[3])
        if itnum > maxit; maxit = itnum; end
    end
    for line in lines
        if !occursin("strength function ",line); continue; end
        if split(line)[3] != string(maxit); continue; end
        tl = split(line)
        n = parse(Int,tl[4])
        Energy = parse(Float64,tl[5])
        S = parse(Float64,tl[6])
        BGT = qGT^2 * S
        logft = log10( K / (gAV^2 * BGT) )
        Ex = Energy-daughter.Egs
        @assert Ex >= 0.0 "Energy $Energy must be higher than daughter Egs $(daughter.Egs)"        
        f0_expQ = Fermi_integral(Qvals[1]-Ex,dZ,R)
        f0_thoQ = Fermi_integral(Qvals[2]-Ex,dZ,R)
        hl_expQ =  K / (gAV^2 * BGT) / f0_expQ
        hl_thoQ =  K / (gAV^2 * BGT) / f0_thoQ
        push!(data,bgtdata(n,Energy,Ex,J2,prty,S,BGT,logft,f0_expQ,f0_thoQ,hl_expQ,hl_thoQ))
        count += 1
        if count == maxit
            break
        end
    end
    if n_eigen != maxit 
        @warn "n_eigen $n_eigen != maxit $maxit  The bgt file $fn may be not completed, or the number of possible states with the given total J is smaller than the specified n_eigen $n_eigen."
    end
    return data
end

##################
## First Forbidden
##################

mutable struct bFFdata
    J2::Int64
    prty::String
    Ex::Float64
    Qd_exp::Float64
    Qd_tho::Float64
    M0rs::Float64
    M0sp::Float64
    M0rsd::Float64
    M1r::Float64
    M1rs::Float64
    M1p::Float64
    M1rd::Float64
    M1rsd::Float64
    M2rs::Float64
    f0_expQ::Float64
    f0_thoQ::Float64
    hl_expQ::Float64
    hl_thoQ::Float64
end

function eval_bFF_files(fns,qfactors::quenching_factors,parent,daughter,Qvals;verbose=false,nonrela=true)
    dZ = daughter.nuc.Z
    @assert dZ <= 137 "Z = $dZ > 137 is not assumed in Fermi's beta-decay theory"
    A = daughter.nuc.A
    dprty = ifelse(parent.states[1].prty==1,"-","+")
    lambda_Ce = hc / Me
    mu1 = lambda2 = 1.0
    R_in_nu = radius_nuc_for_FermiIntegarl(A) / lambda_Ce
    xi = fsalpha*dZ/(2*R_in_nu)
    gamma1 = sqrt(1.0-(fsalpha*dZ)^2)
    data = Dict{Float64,bFFdata}()
    for fn in fns
        read_bff_file!(fn,qfactors,parent,daughter,Qvals,data, verbose)
    end
    fis = zeros(Float64,4)
    keylist = [ tmp for tmp in keys(data)]
    sorted_keys =  [ keylist[i] for i in sortperm(keylist) ]
    Ecoul = Ecoulomb_SM(parent.nuc.Z,parent.nuc.N) - Ecoulomb_SM(daughter.nuc.Z,daughter.nuc.N)
    deltam = (Mn - Mp - Me)
    for tkey in sorted_keys
        ffobj = data[tkey]
        ffobj.prty = dprty
        for use_expQ in [true,false]
            qval = ifelse(use_expQ,ffobj.Qd_exp,ffobj.Qd_tho)
            J2 = ffobj.J2
            if qval < 0.0; continue;end
            Egamma = qval - Ecoul - deltam
            W0 = qval / Me + 1.0
            fis[1] = Fermi_integral(qval,dZ,R_in_nu,-1)
            fis[2] = Fermi_integral(qval,dZ,R_in_nu,0)
            fis[3] = Fermi_integral(qval,dZ,R_in_nu,1)
            fis[4] = Fermi_integral(qval,dZ,R_in_nu,2)
            # Rank 0:
            zeta0 = ffobj.M0sp + xi * ffobj.M0rsd + 1/3 * ffobj.M0rs * W0
            K0_0 = zeta0^2 + 1/9 * ffobj.M0rs^2
            K0_m1 = -2/3 * mu1 * gamma1 * zeta0 * ffobj.M0rs
            f_0 = K0_m1 * fis[1] + K0_0 * fis[2]
            # Rank 1:
            x = ffobj.M1r
            xd = ffobj.M1rd
            u = ffobj.M1rs
            ud = ffobj.M1rsd
            xidy = ifelse(nonrela,Egamma/Me*x,ffobj.M1p)
            Y = xidy - xi * (ud+xd)
            zeta1 = Y + 1/3 * (u-x) * W0
            K1_0 = zeta1^2 + 1/9 * (x + u)^2 - 4/9 * mu1*gamma1*u*(x + u)
            K1_0 += 1/18 * W0^2 * (2*x + u)^2 - 1/18 * lambda2 * (2*x - u)^2 
            K1_1 = -4/3*u*Y - 1/9 * W0 * (4*x^2 + 5*u^2)
            K1_m1= 2/3 * mu1 * gamma1 * zeta1 * (x+u)
            K1_2 = 1/18 * (8*u^2 + (2*x+u)^2 + lambda2*(2*x-u)^2)
            f_1 =  K1_m1 * fis[1] + K1_0 * fis[2] + K1_1 * fis[3] + K1_2 * fis[4]
            # Rank 2:
            z = ffobj.M2rs
            K2_0  =  1/12 * z^2 * (W0^2 - lambda2)
            K2_1  = -1/6 * z^2 * W0
            K2_2  =  1/12 * z^2 * (1.0 + lambda2)
            f_2 = K2_0 * fis[2] + K2_1 * fis[3]  + K2_2 * fis[4]
            fFF = f_0 + f_1 + f_2
            thalf = K/fFF
            logft = log10(fis[2]*thalf) 
            if verbose                
                println("Energy $tkey Ex.",@sprintf("%8.3f",tkey-daughter.Egs)," J2",@sprintf("%4i",J2),ffobj.prty)
                println("M0s,M0t,M0s'", @sprintf("%12.4e", ffobj.M0rs),@sprintf("%12.4e", ffobj.M0sp),@sprintf("%12.4e", ffobj.M0rsd))
                println("x, u, xi'y  ", @sprintf("%12.4e", x),@sprintf("%12.4e", u),@sprintf("%12.4e", xidy))
                println("x', u' z    ", @sprintf("%12.4e", xd),@sprintf("%12.4e", ud),@sprintf("%12.4e", z) )
                println(@sprintf("zeta0 %12.4e breakdown %12.4e %12.4e %12.4e", zeta0, ffobj.M0sp, xi * ffobj.M0rsd, 1/3 * ffobj.M0rs * W0))
                println("K0_01,K0_0 ",@sprintf("%12.4e", K0_m1),@sprintf("%12.4e",K0_0))
                #println("K2_0,K2_1,K2_2 ",@sprintf("%12.4e", K2_0),@sprintf("%12.4e",K2_1),@sprintf("%12.4e", K2_2) )
                println("  W0 ",@sprintf("%12.4e",W0), " qval $qval Egamma $Egamma Ecoul $Ecoul")
                println("f0, f1, f2  ",@sprintf("%12.4e", f_0),@sprintf("%12.4e", f_1),@sprintf("%12.4e", f_2) )
                println("log10(f0t) ", @sprintf("%6.2f",logft) ,@sprintf("%10.4f",logft))
                println("thalf ",@sprintf("%10.4f",thalf))
                println("")
            end
            if use_expQ
                ffobj.f0_expQ = fis[2]
                ffobj.hl_expQ = thalf 
            else
                ffobj.f0_thoQ = fis[2]
                ffobj.hl_thoQ = thalf 
            end
        end
    end
    return sorted_keys,data
end

# """
#     get_wf_if_order(fn,parent,daughter)
    
# Since the logfile by KSHELL is given as log_NucA_NucB_sntname_tr_JpiA_JpiB.txt,
# """
# function get_wf_if_order(fn,parent,daughter)
#     pnuc = parent.nuc.cnuc
#     dnuc = daughter.nuc.cnuc
#     if occursin(pnuc*"_"*dnuc, fn)
#         return true
#     elseif occursin(dnuc*"_"*pnuc,fn)
#         return false
#     else
#         #@error "filename:$fn does not include $(pnuc)_$(dnuc) nor $(dnuc)_$(pnuc)"
#         return nothing
#     end
# end

"""

Note: 
When you evaluate FF transition matrix elements, you should have prepared using KSHELL log_*_tr* files under the following order: Left = Daughter & Right = Parent.
In the f0-integrals for FF transitions the reduced matrix elements such as, M2rs=<f|| [rC1 x sigma]^(2) t^- ||i>.
This has a nontrivial sign (with respect to flip w.r.t. left and right), which originates from the conventions in the KSHELL code.
"""
function read_bff_file!(fn,qfactors::quenching_factors,parent,daughter,Qvals,data,verbose)
    dZ = daughter.nuc.Z
    @assert dZ <= 137 "Z = $dZ > 137 is not assumed in Fermi's beta-decay theory"
    A = daughter.nuc.A
    dN = A - dZ
    lambda_Ce = hc / Me    
    mass_n_nu = 0.5*(Mp + Mn) / Me
    if !isfile(fn)
        @warn "file $fn is not found"
        return nothing
    end
    f = open(fn,"r"); lines = readlines(f); close(f)
    hit = 0; rank = -1
    FFchannel = ""
    tar_text1 = "First Forbidden"
    tar_text2 = "First forbidden"
    order = false
    for line in lines
        if occursin("FN_LOAD_WAVE_L",line)
            nuc_l = strip(split(split(replace(split(line,"=")[end],"\""=>""),"/")[end],"_")[1])
            order = ifelse(nuc_l == parent.nuc.cnuc,true,false)
            #@assert order == false "When you evaluate FF transitions, wavefunctions must be ordered as <daughter|O|parent>."
        end
        if occursin(tar_text1,line) || occursin(tar_text2,line) 
            hit += 1
            rank = parse(Int,split(split(line,",")[2])[2])
            FFchannel = split(split(split(line,",")[end],"=")[1])[end]
            continue
        end
        if hit == 0; continue;end
        if occursin("E_gam  D_Ec",line) || occursin("total elapsed",line);break;end
        if !( occursin("(",line) && occursin(")",line)); continue; end
        if occursin("B(EM )->",line); continue;end  
        line = replace(line, "("=>" ")
        line = replace(line, ")"=>" ")        
        tl = split(line)
        J2f, NJf, Ef, J2i, NJi, Ei, Ex, Mred, B1, B2, Mom = tl
        Nth_parent = parse(Int,ifelse(order, NJf, NJi))
        if Nth_parent != 1; continue; end # If the target state is not the lowest (within the same Jpi), skip
        Edaughter = parse(Float64,ifelse(order,Ei,Ef))
        J2f = parse(Int,J2f)
        J2i = parse(Int,J2i)
        J2 = ifelse(order,J2i,J2f)
        Mred = parse(Float64,Mred)
        C_J = 1.0 / sqrt( ifelse(order,J2f,J2i) + 1)
        Bval = (C_J * Mred)^2
        Ex_calc = Edaughter - daughter.Egs
        @assert Ex_calc >= 0.0 " Ex_calc $Ex_calc must be >= 0 "
        if order 
            B1 = parse(Float64,B1)    
            @assert abs(B1-Bval) < 1.e-3 "B1 $B1 B(FF) $Bval"
        else
            B2 = parse(Float64,B2)
            @assert abs(B2-Bval) < 1.e-3 "B2 $B2 B(FF) $Bval"
        end
        Qd_exp = Qvals[1] - Ex_calc
        Qd_tho = Qvals[2] - Ex_calc
        # if verbose
        #     println("rank $rank channel $FFchannel Qexp $(@sprintf("%6.2f", Qd_exp)) Qtho $(@sprintf("%6.2f", Qd_tho)) Bval $(@sprintf("%7.3f", Bval))  Mred $(@sprintf("%7.3f", Mred))",
        #             "Edaughter $(@sprintf("%9.3f", Edaughter)) Ex_calc $(@sprintf("%9.3f", Ex_calc)) order $order")
        # end
        if !haskey(data, Edaughter)
            data[Edaughter] = bFFdata(J2,"",Ex_calc,Qd_exp,Qd_tho,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0)
        end
        ME_FF_func!(FFchannel,qfactors,C_J,Mred,lambda_Ce,mass_n_nu,data[Edaughter])
    end
    return nothing
end

"""
    ME_FF_func!(chan,qfactors,C_J,Mred,lambda_Ce,mass_n_nu,dict)

Function to calc FF matrix elements.
Since M0rs,M1r,M1rs,M2rs (M0sp,M1p) are evaluated in fm (1/fm) in KSHELL, it is needed to multiply 1/lambda_Ce (lambda_Ce) to get transition matrix element in natural unit.
"""
function ME_FF_func!(chan,qfactors,C_J,Mred,lambda_Ce,mass_n_nu,dict)
    if chan == "M0rs" #M0rs=<f||[rC1 x sigma]^(0)t^-||i>
        dict.M0rs = lam_gav * sqrt(3.0) * C_J * Mred / lambda_Ce * qfactors.qMS0
    elseif chan == "M0sp"  #M0sp=<f||[sigma x nabla]^(0)t^-||i> 
        dict.M0sp = -lam_gav * sqrt(3.0) * C_J / mass_n_nu * Mred * lambda_Ce * qfactors.qMT0
    elseif chan == "M0rsd" # M0rsd=<f||2/3I[rC1 x sigma]^(0)t^-||i>
        dict.M0rsd = lam_gav * sqrt(3.0) * C_J * Mred / lambda_Ce * qfactors.qMS0
    elseif chan == "M1r" #xi'y M1r =<f||[rC1 t^-||i> 
        dict.M1r = - C_J * Mred / lambda_Ce * qfactors.qx
    elseif chan == "M1rs" #u M1rs=<f||[rC1 x sigma]^(1)t^-||i> 
        dict.M1rs = lam_gav * sqrt(2.0)* C_J * Mred / lambda_Ce * qfactors.qu
    elseif chan == "M1p" #  M1p =<f||nabla t^-||i> 
        dict.M1p = - C_J / mass_n_nu * Mred * lambda_Ce  * qfactors.qx
    elseif chan == "M1rd" # M1rd =<f||2/3I[rC1 t^-||i> " 
        dict.M1rd = - C_J * Mred / lambda_Ce  * qfactors.qx
    elseif chan == "M1rsd" # M1rsd=<f||2/3I[rC1 x sigma]^(1)t^-||i>
        dict.M1rsd = lam_gav * sqrt(2.0) * C_J * Mred / lambda_Ce * qfactors.qu
    elseif chan == "M2rs" #M2rs=<f|| [rC1 x sigma]^(2) t^- ||i>
        dict.M2rs = -2.0 * lam_gav * C_J * Mred / lambda_Ce * qfactors.qz
    else
        @error "FFchannel $chan is not supported"
    end
end 
