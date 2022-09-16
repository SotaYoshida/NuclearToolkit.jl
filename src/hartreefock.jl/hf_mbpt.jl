"""
    HF_MBPT2(binfo,modelspace,fp,fn,e1b_p,e1b_n,Chan2b,Gamma)
Calculate 2nd order correction to HF energy
```math
E^{(2)} = \\frac{1}{4}\\sum_{abij} \\frac{\\bar{H}^{[2]}_{abij} \\bar{H}^{[2]}_{ijab}}{\\epsilon^{ab}_{ij}} = \\frac{1}{4} \\sum_{\\tilde{a}\\tilde{b}\\tilde{i}\\tilde{j}}\\sum_{\\{m\\}}\\sum_{JJ'MM'}
\\frac{{}^J\\bar{H}^{[2]}_{\\tilde{a}\\tilde{b}\\tilde{i}\\tilde{j}} {}^J\\bar{H}^{[2]}_{\\tilde{i}\\tilde{j}\\tilde{a}\\tilde{b}}}{\\epsilon^{ab}_{ij}}
(j_a j_b m_a m_b|J M)
(j_a j_b m_a m_b|J' M')
(j_i j_j m_i m_j|J M)
(j_i j_j m_i m_j|J' M')
```
```math
=
\\frac{1}{4} \\sum_{\\tilde{a}\\tilde{b}\\tilde{i}\\tilde{j}}\\sum_{JJ'MM'}
\\frac{{}^J\\bar{H}^{[2]}_{\\tilde{a}\\tilde{b}\\tilde{i}\\tilde{j}} {}^{J'}\\bar{H}^{[2]}_{\\tilde{i}\\tilde{j}\\tilde{a}\\tilde{b}}}{\\epsilon^{ab}_{ij}}
\\delta_{JJ'} \\delta_{MM'} = \\frac{1}{4} \\sum_{\\tilde{a}\\tilde{b}\\tilde{i}\\tilde{j}}\\sum_{J}(2J+1)
\\frac{{}^J\\bar{H}^{[2]}_{\\tilde{a}\\tilde{b}\\tilde{i}\\tilde{j}} {}^J\\bar{H}^{[2]}_{\\tilde{i}\\tilde{j}\\tilde{a}\\tilde{b}}}{\\epsilon^{ab}_{ij}}
```
"""
function HF_MBPT2(binfo,modelspace,fp,fn,e1b_p,e1b_n,Chan2b,Gamma;verbose=false,io=stdout)
    p_sps = modelspace.p_sps
    n_sps = modelspace.n_sps
    sps = modelspace.sps
    holes = modelspace.holes
    particles = modelspace.particles
    nchan = length(Chan2b)
    EMP2 = EMP2_1b = EMP2_pp = EMP2_pn = EMP2_nn = 0.0    
    ### One-body correction (basically zero)
    for pn = 1:2     
        f = ifelse(pn==1,fp,fn)
        tsps = ifelse(pn==1,p_sps,n_sps)   
        for α in holes[pn]
            iα = div(α,2) + ifelse(pn==1,1,0)
            eα = f[iα,iα]
            Nα = tsps[iα].j + 1
            for i in particles[pn]
                ii = div(i,2) + ifelse(pn==1,1,0)
                ei = f[ii,ii]
                EMP2 += Nα * f[ii,iα]^2 / (eα - ei)
                EMP2_1b += Nα * f[ii,iα]^2 / (eα - ei)
            end
        end
    end
    ## Two-body correction <hh(αα_)| pp(ββ_)>
    for ch = 1:nchan
        Gam = Gamma[ch]
        tmp = Chan2b[ch]
        Tz = tmp.Tz; prty=tmp.prty; J=tmp.J; kets = tmp.kets
        npq = length(kets)
        for ib = 1:npq
            α, α_ = kets[ib]
            oα = sps[α]; oα_= sps[α_]
            nafac = oα.occ * oα_.occ
            if nafac == 0.0; continue;end
            if (oα.tz + oα_.tz != Tz);continue;end
            iα = div(α,2) + α%2 
            iα_= div(α_,2) + α_%2 
            e1b_α  = ifelse(α%2==1,e1b_p,e1b_n)
            e1b_α_ = ifelse(α_%2==1,e1b_p,e1b_n)
            for ik = 1:npq
                β, β_ = kets[ik]
                if sps[β].occ + sps[β_].occ !=0.0;continue;end 
                if (sps[β].tz + sps[β_].tz != Tz);continue;end
                iβ = div(β,2) + β%2 
                iβ_= div(β_,2) + β_%2 
                e1b_β  = ifelse(β%2==1,e1b_p,e1b_n)
                e1b_β_ = ifelse(β_%2==1,e1b_p,e1b_n)
                nume = Gam[ib,ik]^2
                deno = e1b_α[iα] + e1b_α_[iα_] - e1b_β[iβ] - e1b_β_[iβ_]
                EMP2 += (2*J+1) * nume/deno * nafac
                if Tz == -2; EMP2_pp +=  (2*J+1) * nume/deno;end
                if Tz ==  0; EMP2_pn +=  (2*J+1) * nume/deno;end
                if Tz ==  2; EMP2_nn +=  (2*J+1) * nume/deno;end
            end
        end
    end
    if verbose
        println(io,"EMP2 ",@sprintf("%9.3f",EMP2)," 1b ",@sprintf("%9.3f",EMP2_1b),
                " pp ",@sprintf("%9.3f",EMP2_pp)," pn ",@sprintf("%9.3f",EMP2_pn)," nn ",@sprintf("%9.3f",EMP2_nn))
    end
    return EMP2
end

"""
    HF_MBPT3(binfo,modelspace,e1b_p,e1b_n,Chan2b,dict_2b_ch,dict6j,Gamma,to;io=stdout)

Calculate 2nd order correction to HF energy
```math
E^{(3)}=
\\frac{1}{8} 
\\sum_{\\tilde{a}\\tilde{b}\\tilde{c}\\tilde{i}\\tilde{j}\\tilde{k}} \\sum_{J}(2J+1)
\\frac{
{}^{J}\\bar{H}^{[2]}_{\\tilde{a}\\tilde{b}\\tilde{i}\\tilde{j}} 
{}^{J}\\bar{H}^{[2]}_{\\tilde{i}\\tilde{j}\\tilde{c}\\tilde{d}}
{}^{J}\\bar{H}^{[2]}_{\\tilde{c}\\tilde{d}\\tilde{a}\\tilde{b}} 
}
{\\epsilon^{\\tilde{a}\\tilde{b}}_{\\tilde{i}\\tilde{j}}\\epsilon^{\\tilde{c}\\tilde{d}}_{\\tilde{i}\\tilde{j}}}
+
\\frac{1}{8} 
\\sum_{\\tilde{a}\\tilde{b}\\tilde{i}\\tilde{j}\\tilde{k}\\tilde{l}} \\sum_{J}(2J+1)
\\frac{
{}^{J}\\bar{H}^{[2]}_{\\tilde{i}\\tilde{j}\\tilde{a}\\tilde{b}} 
{}^{J}\\bar{H}^{[2]}_{\\tilde{a}\\tilde{b}\\tilde{k}\\tilde{l}}
{}^{J}\\bar{H}^{[2]}_{\\tilde{k}\\tilde{l}\\tilde{i}\\tilde{j}} 
}
{\\epsilon^{\\tilde{a}\\tilde{b}}_{\\tilde{i}\\tilde{j}}\\epsilon^{\\tilde{c}\\tilde{d}}_{\\tilde{i}\\tilde{j}}}
-\\sum_{\\tilde{a}\\tilde{b}\\tilde{c}\\tilde{i}\\tilde{j}\\tilde{k}} \\sum_{J}(2J+1)
\\frac{
{}^JH^{XC}_{\\tilde{a}\\tilde{i}\\tilde{j}\\tilde{b}}
{}^JH^{XC}_{\\tilde{j}\\tilde{b}\\tilde{k}\\tilde{c}}
{}^JH^{XC}_{\\tilde{k}\\tilde{c}\\tilde{a}\\tilde{i}}
}{
\\epsilon^{\\tilde{a}\\tilde{b}} \\epsilon_{\\tilde{i}\\tilde{j}}
\\epsilon^{\\tilde{a}\\tilde{c}} \\epsilon_{\\tilde{k}\\tilde{j}}
}
```

Ref. Many-Body Methods in Chemistry and Physics by Isaiah Shavitt and Rodney J. Bartlett (2009, Cambridge Molecular Science).
More details can be found in e.g. Dr. thesis by A.Tichai (2017, TU Darmstadt).
"""
function HF_MBPT3(binfo,modelspace,e1b_p,e1b_n,Chan2b,dict_2b_ch,dict6j,Gamma,to;verbose=false,io=io)
    p_sps = modelspace.p_sps
    n_sps = modelspace.n_sps
    sps = modelspace.sps
    holes = modelspace.holes
    particles = modelspace.particles
    if particles[1] == [ ] && particles[2] == []
        println("modelspace(emax) is too small for HFMBPT(3)")
        return 0.0
    end
    e1b = Float64[ ] 
    for i in eachindex(e1b_p)
        push!(e1b,e1b_p[i]); push!(e1b,e1b_n[i])
    end
    EMP3 = EMP3_pp = EMP3_ph = EMP3_hh = 0.0
    nchan = length(Chan2b)
    """
    E3_pp: h=>(α,α') p=>(ββ'γγ') <αα'|H|ββ'><ββ'|H|γγ'><γγ'|H|αα'> /(eα+eα'-eβ-eβ')(eα+eα'-eγ-eγ')
    E3_hh: h=>(α,α'γγ') p=>(ββ') <αα'|H|ββ'><ββ'|H|γγ'><γγ'|H|αα'> /(eα+eα'-eβ-eβ')(eγ+eγ'-eβ-eβ')
    """
    for ch = 1:nchan
        Gam = Gamma[ch]
        tbc = Chan2b[ch]; J=tbc.J; kets = tbc.kets
        npq = length(kets)
        for i = 1:npq
            α, α_ = kets[i]
            nafac = sps[α].occ * sps[α_].occ
            if nafac == 0; continue;end
            for j = 1:npq                       
                β, β_ = kets[j]                
                if sps[β].occ + sps[β_].occ == 0
                    v1 = Gam[i,j]
                    for k = 1:npq
                        γ, γ_ = kets[k]
                        if sps[γ].occ + sps[γ_].occ != 0.0; continue;end
                        v2 = Gam[j,k]
                        v3 = Gam[k,i]
                        nume = v1 * v2 * v3
                        deno = (e1b[α] + e1b[α_] - e1b[β] - e1b[β_]) * (e1b[α] + e1b[α_] - e1b[γ] - e1b[γ_])
                        EMP3_pp += (2*J+1) * nume/deno # * nafac
                    end
                end
                nbfac = sps[β].occ * sps[β_].occ
                if nbfac != 0
                    v1 = Gam[i,j]
                    for k = 1:npq
                        γ, γ_ = kets[k]
                        if sps[γ].occ + sps[γ_].occ != 0; continue;end
                        v2 = Gam[j,k]
                        v3 = Gam[k,i]
                        nume = v1 * v2 * v3
                        deno = (e1b[α] + e1b[α_] - e1b[γ] - e1b[γ_]) * (e1b[β] + e1b[β_] - e1b[γ] - e1b[γ_])
                        EMP3_hh += (2*J+1) * nume/deno  #* nbfac
                    end 
                end                    
            end
        end
    end   
    allhs = vcat(holes[1],holes[2])
    allps = vcat(particles[1],particles[2])
    nthre = nthreads()
    keychs = [ zeros(Int64,3) for i=1:nthre]
    Ethreads = zeros(Float64,nthre)      
    @threads for idxa in eachindex(allps)
        a = allps[idxa]
        threid = threadid()
        keych = keychs[threid]
        oa = sps[a]; la = oa.l; ja = oa.j; tz_a = oa.tz
        Etmp = 0.0 
        for i in allhs
            oi = sps[i]; li = oi.l; ji = oi.j; tz_i = oi.tz
            for j in allhs
                oj = sps[j]; lj = oj.l; jj = oj.j; tz_j = oj.tz
                tz_ij = tz_i + tz_j
                Jmin = div(abs(ja-jj),2)
                Jmax = div(ja+jj,2)
                for totJ = Jmin:Jmax                           
                    if tri_check(ja,jj,totJ*2)==false;continue;end
                    Jfac = (2.0*totJ+1.0)
                    tdict6j = dict6j[totJ+1]  
                    ehole = e1b[i]+e1b[j]
                    prty_ij = (-1)^(li+lj)
                    for b in allps                        
                        ob = sps[b]; lb = ob.l; jb = ob.j; tz_b = ob.tz
                        tz_ab = tz_a + tz_b                                                    
                        if tz_ab != tz_ij;continue;end
                        if (-1)^(la+lb) != prty_ij;continue;end                            
                        if tri_check(ji,jb,totJ*2)==false;continue;end
                        keych[1] = tz_a + tz_b; keych[2] = (-1)^(la+lb)
                        v1 = vPandya(a,b,i,j,ja,jb,ji,jj,totJ,dict_2b_ch,tdict6j,Gamma,keych) 
                        if v1 == 0.0;continue;end
                        v1 = v1 / (ehole - e1b[a] - e1b[b])
                        for k in allhs                                 
                            ok = sps[k]; lk = ok.l; jk = ok.j; tz_k = ok.tz
                            prty_kb = (-1)^(lk+lb)
                            for c in allps 
                                oc = sps[c]
                                tz_c = oc.tz
                                if tz_i + tz_c != tz_k + tz_b;continue;end
                                if tz_k + tz_j != tz_a + tz_c;continue;end  
                                lc = oc.l
                                if (-1)^(li+lc) != prty_kb;continue;end                                
                                jc = oc.j
                                if tri_check(jk,jc,totJ*2)==false;continue;end
                                keych[1] = tz_i + tz_c; keych[2] = prty_kb # prty_ic 
                                v2 = vPandya(i,c,k,b,ji,jc,jk,jb,totJ,dict_2b_ch,tdict6j,Gamma,keych)
                                if v2==0.0;continue;end
                                keych[1] = tz_k + tz_j; keych[2] = (-1)^(lk+lj)  
                                v3 = vPandya(k,j,a,c,jk,jj,ja,jc,totJ,dict_2b_ch,tdict6j,Gamma,keych)
                                v3 = v3 / (e1b[k] + e1b[j] -e1b[a] -e1b[c])
                                Etmp += Jfac * v1 * v2 * v3                                   
                            end
                        end
                    end
                end
            end
        end
        Ethreads[threid] += Etmp
    end
    EMP3_ph = sum(Ethreads)
    EMP3 = EMP3_pp + EMP3_hh + EMP3_ph 
    if verbose
        println(io,"pp ",@sprintf("%9.3f",EMP3_pp)," hh ",@sprintf("%9.3f",EMP3_hh),
                " ph ",@sprintf("%9.3f",EMP3_ph)," EMP3 ",@sprintf("%9.3f",EMP3))
    end 
    return EMP3
end

function get_intkey_2(a,b;ofst=1000)
    return ofst*a + b
end

"""
    vPandya(a,b,c,d,ja,jb,jc,jd,totJ,dict_2b_ch,tdict6j,Gamma,keych,key6j,keyab;verbose=false) 

returns generalized Pandya transformed matrix element:
```math
\\tilde{V}^J_{ajib} = -\\sum_{J'} [J'] 
\\begin{Bmatrix} 
j_a & j_j & J \\\\ j_i & j_d & J'
\\end{Bmatrix}
V^{J'}_{abij}
```
"""
function vPandya(a,b,c,d,ja,jb,jc,jd,totJ,dict_2b_ch,tdict6j,Gamma,keych;verbose=false) 
    Jmin = div(max(abs(ja-jb),abs(jc-jd)),2)
    Jmax = div(min(ja+jb,jc+jd),2)
    nkeyab = get_intkey_2(a,b); nkeycd = get_intkey_2(c,d)
    if a > b # this can happen only when Tz=0 now
        nkeyab = get_intkey_2(b,a)
    end
    if c > d # this can happen only when Tz=0 now
        nkeycd = get_intkey_2(d,c)
    end
    v = 0.0
    @inbounds for Jp = Jmin:Jmax
        if Jp % 2 == 1 && (a==b || c==d);continue;end          
        nkey = get_nkey_from_key6j(ja,jd,jc,jb,Jp); t6j = tdict6j[nkey]
        if t6j == 0.0; continue;end
        keych[3] = Jp
        vch = dict_2b_ch[keych].Vch; vdict = dict_2b_ch[keych].Vdict
        norfac = ifelse(a == b,sqrt(2.0),1.0) *  ifelse(c == d,sqrt(2.0),1.0)        
        if a > b; norfac *= (-1)^(div(ja+jb,2)+Jp+1); end 
        if c > d; norfac *= (-1)^(div(jc+jd,2)+Jp+1); end
        ib = vdict[nkeyab]; ik = vdict[nkeycd]   
        v -= t6j * (2*Jp+1) * norfac * Gamma[vch][ib,ik]      
    end
    return v
end
