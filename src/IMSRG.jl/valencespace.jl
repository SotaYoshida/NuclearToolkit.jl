
"""
    check_valence_space(HFobj,valencespace)
check validity of specified valence space
"""
function check_valence_space(HFobj,valencespace)
    v = Int64[ ]
    if length(valencespace)==0
        return v
    elseif typeof(valencespace)== String
        check_major_valencespace(valencespace,HFobj,v)
    else
        check_str_valencespace(valencespace,HFobj,v)
    end
    return v
end

"""
    check_major_valencespace(str::String,HFobj,v)

Function to check valence space and overwrite `v` and `q` fields of SingleParticleState
The valencespace is specified by argument `str` (e.g. "p-shell") 
"""
function check_major_valencespace(str::String,HFobj,v)
    sps = HFobj.modelspace.sps
    p_sps = HFobj.modelspace.p_sps
    n_sps = HFobj.modelspace.n_sps
    target_nljtz = Vector{Int64}[ ]
    if str == "p-shell" 
        target_nljtz = [[0,1,1,-1],[0,1,3,-1], [0,1,1,1],[0,1,3,1]]
    elseif str == "sd-shell"
        target_nljtz = [[1,0,1,-1],[0,2,3,-1],[0,2,5,-1],[1,0,1,1],[0,2,3,1],[0,2,5,1]]
    elseif str == "pf-shell"
        target_nljtz = [[1,1,1,-1],[1,1,3,-1],[0,3,5,-1],[0,3,7,-1],[1,1,1,1],[1,1,3,1],[0,3,5,1],[0,3,7,1]]
    elseif str == "psd-shell"
        target_nljtz = [[0,1,1,-1],[0,1,3,-1], [1,0,1,-1],[0,2,3,-1],[0,2,5,-1],[0,1,1,1],[0,1,3,1],[1,0,1,1],[0,2,3,1],[0,2,5,1]]
    elseif str == "sdpf-shell"
        target_nljtz = [[1,0,1,-1],[0,2,3,-1],[0,2,5,-1],[1,1,1,-1],[1,1,3,-1],[0,3,5,-1],[0,3,7,-1],[1,0,1,1],[0,2,3,1],[0,2,5,1],[1,1,1,1],[1,1,3,1],[0,3,5,1],[0,3,7,1]]
    elseif str == "pfsdg-shell"
        target_nljtz = [[1,1,1,-1],[1,1,3,-1],[0,3,5,-1],[0,3,7,-1],[2,0,1,-1],[1,2,3,-1],[1,2,5,-1],[0,4,7,-1],[0,4,9,-1],
                        [1,1,1,1],[1,1,3,1],[0,3,5,1],[0,3,7,1],[2,0,1,1],[1,2,3,1],[1,2,5,1],[0,4,7,1],[0,4,9,1]]
    else
        println("modelaspace:$str is not supported!")
        println("If you wannna add $str, please edit function `check_major_valencespace` in `src/IMSRG.jl/valencespace.jl`")
        exit()
    end
    for (i,o) in enumerate(sps)
        idx = div(i,2) + i%2 
        tsps = ifelse(o.tz==-1,p_sps,n_sps)[idx]
        if [o.n,o.l,o.j,o.tz] in target_nljtz
            push!(v,i)
            o.v = tsps.v = true
            o.c = tsps.c = false
        else
            if o.occ == 0.0
                o.q = tsps.q = true
            end
        end
    end
    return nothing
end

"""
    check_str_valencespace(valencespace::Vector{Vector{Int64}},HFobj,v)

check valence space and overwrtie SingleParticleState.v/q

specified by or Vector{Int} (e.g., [[0,1,1,-1],[0,1,3,-1], [0,1,1,1],[0,1,3,1]])
"""
function check_str_valencespace(valencespace::Vector{Vector{Int64}},HFobj,v)
    sps = HFobj.modelspace.sps
    for target in valencespace
        for (i,o) in enumerate(sps)
            if [o.n,o.l,o.j,o.tz] in target
                push!(v,i)
            end
        end
    end
    return nothing
end 
   
"""
    update_vsspace_chs!(HFobj,valencespace,Chan2b)   

overwrite cc/vc/qc/vv/qv/qq channals
"""
function update_vsspace_chs!(HFobj,valencespace,Chan2b)    
    sps = HFobj.modelspace.sps
    spaces = HFobj.modelspace.spaces
    cc = spaces.cc; vc = spaces.vc; qc = spaces.qc
    vv = spaces.vv; qv = spaces.qv; qq = spaces.qq
    for ch in eachindex(Chan2b)
        tbc = Chan2b[ch]
        kets = tbc.kets
        for (ik,ket) in enumerate(kets)
            i,j = ket
            v_i = ifelse(sps[i].v,1,0)
            v_j = ifelse(sps[j].v,1,0)            
            ## cc 
            if sps[i].c && sps[j].c
                add_ch_ket!(ch,ik,cc) 
            end
            ## vc or qc case
            if (sps[i].c || sps[j].c) && !(sps[i].c && sps[j].c)
                if v_i+v_j == 1
                    add_ch_ket!(ch,ik,vc)
                else 
                    add_ch_ket!(ch,ik,qc)
                end
            end
            ## vv,qv,qq
            if v_i + v_j == 2 
                # if you need additional truncations, specify here
                add_ch_ket!(ch,ik,vv)
            elseif v_i + v_j == 1 && (sps[i].q || sps[j].q)
                add_ch_ket!(ch,ik,qv) 
            elseif sps[i].q && sps[j].q
                add_ch_ket!(ch,ik,qq) 
            end
        end       
    end
    return nothing
end

"""
    calc_Eta_smatan!(HFobj,IMSRGobj,Chan2b,dictMono,norms)

``\\eta(s)`` with shell-model atan generator to decouple the specified valence space.
"""
function calc_Eta_smatan!(HFobj::HamiltonianNormalOrdered,IMSRGobj,Chan2b,dictMono,norms)
    MS = HFobj.modelspace
    sps = MS.sps; p_sps = MS.p_sps; n_sps = MS.n_sps    
    key = zeros(Int64,2)
    ## one-body piece a->h i->v/p
    Hs = IMSRGobj.H; f = Hs.onebody;  Gamma = Hs.twobody
    Eta1b = IMSRGobj.eta.onebody
    Eta2b = IMSRGobj.eta.twobody
    Delta = IMSRGobj.denominatorDelta    
    for pn = 1:2
        eta1b = Eta1b[pn]; eta1b .*= 0.0
        tf = f[pn]
        pnrank = ifelse(pn==1,1,3)
        dMono = dictMono[pnrank]
        tz_pn = ifelse(pn==1,-1,1)
        for (a,oa) in enumerate(sps)
            if oa.tz != tz_pn;continue;end
            if oa.q == true;continue;end 
            idx_a = div(a,2) +a%2
            for (i,oi) in enumerate(sps)
                if i==a;continue;end
                if oi.tz != tz_pn; continue;end
                if oi.c; continue;end
                idx_i = div(i,2) + i%2
                nume = 2 * tf[idx_a,idx_i]
                key[1] = a; key[2] = i
                if a > i; key[1] = i; key[2]= a; end
                mono_ai = dMono[key].monopole[1] 
                deno = tf[idx_i,idx_i] -tf[idx_a,idx_a] + (oi.occ-oa.occ)*mono_ai + Delta                
                tmp = 0.5 * atan( nume / deno)
                eta1b[idx_i,idx_a] = tmp
                eta1b[idx_a,idx_i] = -tmp
            end
        end 
    end 

    ## two-body piece, decoupling of core 
    for ch in eachindex(Chan2b)
        Gam = Gamma[ch]; tbc = Chan2b[ch]; Tz = tbc.Tz; kets = tbc.kets
        eta2b = Eta2b[ch]; eta2b .*= 0.0
        pnrank = div(Tz,2) + 2
        ### decoupling of core - valence/qspace
        for ik in eachindex(kets) # to be cc/vc (cq is not considered)
            i,j = kets[ik]  
            oi =sps[i]; oj = sps[j]      
            if !oi.c && !oj.c;continue;end
            ni = oi.occ; nj = oj.occ
            #if ni + nj == 1 && (oi.q || oj.q); continue;end
            if (ni*nj==0.0 && ni+nj!=0.0) && (oi.q || oj.q); continue;end
            tz_ij = oi.tz + oj.tz
            if tz_ij != Tz;continue;end
            for ib in eachindex(kets) # to be vv/qv/qq
                a,b = kets[ib]
                oa = sps[a]; ob = sps[b]
                if oa.c || ob.c; continue;end
                na = oa.occ; nb = ob.occ
                tz_ab = oa.tz + ob.tz
                if tz_ab != Tz;continue;end
                nume = 2 * Gam[ib,ik]
                deno = Get2bDenominator(ch,pnrank,a,b,i,j,na,nb,ni,nj,f,Delta,dictMono,key)
                tmp = 0.5 * atan(nume / deno)                
                eta2b[ib,ik] = tmp
                if ib != ik; eta2b[ik,ib] = -tmp;end
            end
        end          
        ### decoupling of valence - qspace
        for ik in eachindex(kets) # to be vv
            i,j = kets[ik]
            if !sps[i].v || !sps[j].v; continue;end
            ni = sps[i].occ; nj = sps[j].occ 
            for ib in eachindex(kets) # to be qv/qq
                a,b = kets[ib]
                if sps[a].c || sps[b].c; continue;end
                if !sps[a].q && !sps[b].q; continue;end
                na = sps[a].occ; nb = sps[b].occ
                nume = 2 * Gam[ib,ik]
                deno = Get2bDenominator(ch,pnrank,a,b,i,j,na,nb,ni,nj,f,Delta,dictMono,key)
                tmp = 0.5 * atan(nume / deno)                
                eta2b[ib,ik] = tmp
                if ib != ik; eta2b[ik,ib] = -tmp;end
            end
        end 
    end 
    p_sps = HFobj.modelspace.p_sps
    n_sps = HFobj.modelspace.n_sps
    norms[3] = getNorm1b(IMSRGobj.eta.onebody,p_sps,n_sps)
    norms[4] = getNorm2b(IMSRGobj.eta.twobody,Chan2b)    
    return nothing
end

"""
    write_vs_snt(binfo,HFobj,IMSRGobj,Operators,effOps,Chan1b,Chan2bD,vspace)

Function to write out valence space effective interaction in snt (KSHELL/ShellModel.jl) format.
"""
function write_vs_snt(binfo,HFobj,IMSRGobj,Operators,effOps,Chan1b,Chan2bD,vspace)
    hw = binfo.hw
    emax = binfo.emax
    dict_ms2snt = Chan1b.ms2snt
    ref = binfo.ref
    MS = HFobj.modelspace; sps = MS.sps
    cDelta = strip(@sprintf("%5.1f",IMSRGobj.denominatorDelta))
    for i=1:1+length(Operators)
        Op = IMSRGobj.H
        cOp = "Hamiltonian"
        refnuc = ifelse(ref=="core",binfo.nuc.corenuc,binfo.nuc.cnuc)
        fname = "vsimsrg_" * vspace * "_core"*binfo.nuc.corenuc* "ref"*refnuc* "_" * binfo.nuc.cnuc *"_"*"hw"*string(hw)*"e"*string(emax)* "_Delta$cDelta.snt"
        if i !=1
            Op = effOps[i-1]; cOp = Operators[i-1]
            fname = "vsimsrg_" * vspace * "_core" *binfo.nuc.corenuc* "ref"*refnuc * "_" * binfo.nuc.cnuc * "_"*"hw"*string(hw)*"e"*string(emax)*"_Delta"*cDelta*"_"*cOp*".snt"
        end
        io = open(fname,"w")

        ## one-body piece a->h i->v/p
        println(io,"!input interaction => ",binfo.sntf)
        println(io,"!Op:$cOp, zerobody: ",@sprintf("%15.8f",Op.zerobody[1]))
        f = Op.onebody    
        vp_sps = Int[ ]
        vn_sps = Int[ ]
        for (a,oa) in enumerate(sps)
            if !oa.v;continue;end
            pn = ifelse(oa.tz==-1,1,2)
            if pn == 1
                push!(vp_sps,a)
            else
                push!(vn_sps,a)
            end
        end 
        v_sps = Dict{Int64,Int64}() 
        for (i,tmp) in enumerate(vcat(vp_sps,vn_sps))
            v_sps[tmp] = i
        end
        # header #p_sps,#n_sps,cp,cn
        txt = @sprintf("%4i",length(vp_sps)) * @sprintf("%4i",length(vn_sps))  
        txt *=  @sprintf("%4i",binfo.nuc.cZ) * @sprintf("%4i",binfo.nuc.cN)
        println(io, txt)

        for isps in vcat(vp_sps,vn_sps)
            ii = v_sps[isps]
            oi = sps[isps]; n = oi.n; l = oi.l; j = oi.j; tz = oi.tz    
            txt = @sprintf("%3i",ii) * @sprintf("%5i",n) * @sprintf("%4i",l) * @sprintf("%4i",j) * @sprintf("%4i",tz)
            println(io,txt)
        end

        txt = @sprintf("%6i",length(vp_sps)+length(vn_sps))
        txt *=  @sprintf("%3i",0) * @sprintf("%8i", Int(binfo.hw))
        println(io, txt)

        for (a,oa) in enumerate(sps)
            if !oa.v;continue;end
            idx_a = div(a,2) +a%2   
            pn = ifelse(oa.tz==-1,1,2)
            spe = f[pn][idx_a,idx_a]
            ii = v_sps[a]
            txt = @sprintf("%3i",ii) * @sprintf("%3i",ii) * @sprintf("%15.6f",spe) 
            println(io,txt)
        end

        ## two-body piece, decoupling of core 
        Chan2b = Chan2bD.Chan2b
        vv = MS.spaces.vv
        txt = ""
        ntbme = 0
        ntbmes = zeros(Int64,3)
        for ch in keys(vv)
            tbc = Chan2b[ch]
            H2 = Op.twobody[ch]
            tkets = tbc.kets
            idxs_vv = vv[ch]
            J = tbc.J 
            abcd = zeros(Int64,4)
            pnrank = 2 + div(tbc.Tz,2)
            for ib in idxs_vv
                a,b = tkets[ib]
                snt_a = v_sps[a]
                snt_b = v_sps[b]
                flip_ab = ifelse(snt_a>snt_b,true,false)
                phase_ab = 1.0
                if flip_ab
                    ja = sps[a].j; jb = sps[b].j
                    phase_ab *= (-1)^(div(ja+jb,2)+J+1)
                    abcd[1] = b; abcd[2] = a
                else
                    abcd[1] = a; abcd[2] = b
                end
                #sq_ab = ifelse(a==b,sqrt(2.0),1.0)* phase_ab 
                for ik in idxs_vv
                    if ib > ik;continue;end
                    c,d = tkets[ik]
                    snt_c = v_sps[c]
                    snt_d = v_sps[d]
                    flip_cd = ifelse(snt_c>snt_d,true,false)
                    phase_cd = 1.0
                    if flip_cd
                        jc = sps[c].j; jd = sps[d].j
                        phase_cd *= (-1)^(div(jc+jd,2)+J+1)                
                    end
                    abcd[1] = snt_a; abcd[2] =  snt_b; abcd[3] = snt_c; abcd[4] = snt_d
                    if flip_ab
                        abcd[1] = snt_b; abcd[2] = snt_a
                    end
                    if flip_cd
                        abcd[3] = snt_d; abcd[4] = snt_c
                    end
                    if abcd[1] > abcd[3]
                        ta = abcd[3]; tb = abcd[4]
                        abcd[3] = abcd[1]; abcd[4] = abcd[2]
                        abcd[1] = ta; abcd[2] = tb
                    end
                    tbme = H2[ib,ik] * phase_ab * phase_cd 

                    txt *=  @sprintf("%3i",abcd[1]) * @sprintf("%3i",abcd[2]) 
                    txt *= @sprintf("%3i",abcd[3]) * @sprintf("%3i",abcd[4])
                    txt *= @sprintf("%4i",J) * @sprintf("%12.6f",tbme)
                    txt *= "\n"
                    ntbmes[pnrank] += 1
                    ntbme += 1
                end
            end
        end          
        txt2b = @sprintf("%10i",ntbme) * @sprintf("%5i",0) * @sprintf("%12.5f",binfo.hw)
        println(io,txt2b)
        println(io,rstrip(txt))
        close(io)
    end
    return nothing
end
