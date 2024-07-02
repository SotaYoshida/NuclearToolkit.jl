"""
    get_ZNref(ref,Z,N,corenuc)
get ``Z`` and ``N`` of the target reference
"""
function get_ZNref(ref::String,Z::Int,N::Int,corenuc::String)
    Zref = Z; Nref = N
    cZ = cN = 0
    if ref =="core"
        Zr,Nr = cZN_from_corenuc(Z,N,corenuc)
        cZ = Zref; cN = Nref
    elseif ref == "nucl"
        nothing
    else
        println("ref must be 'core' or 'nucl'");exit()
    end
    return Zr,Nr
end
"""
    def_nuc(nuc::Vector{Int},ref,corenuc)
constructor of `nucleus` strict from given `Z`,`N`,`ref`,`corenuc`
"""
function def_nuc(nuc::Vector{Int},ref::String,corenuc::String)
    Z,N = nuc
    el = nuclist[Z]
    A = Z+N
    cnuc = el * string(A)
    Zref = Z; Nref = N    
    cZ,cN = cZN_from_corenuc(Zref,Nref,corenuc)
    if ref == "core"
        Zref = cZ; Nref = cN
    end    
    Aref = Zref + Nref
    return nucleus(Zref,Nref,A,Aref,el,cnuc,cZ,cN,corenuc)
end 
"""
    def_nuc(Z,N,ref,corenuc)
constructor of `nucleus` strict from given `cnuc`,`ref`,`corenuc`
"""
function def_nuc(cnuc::String,ref::String,corenuc::String)
    reg = r"[0-9]+"
    A = match(reg,cnuc).match 
    el = replace(cnuc, A=>"")
    Z = 0
    for tZ=1:length(nuclist)
        if el == nuclist[tZ]
            Z = tZ;break
        end
    end 
    A = parse(Int64,A); N = A-Z
    Zref = Z; Nref = N    
    cZ,cN = cZN_from_corenuc(Zref,Nref,corenuc)
    if ref == "core"
        Zref = cZ; Nref = cN
    end    
    Aref = Zref + Nref
    return nucleus(Zref,Nref,A,Aref,el,cnuc,cZ,cN,corenuc)
end 

"""
    cZN_from_corenuc(rZ,rN,corenuc)
get ``Z`` and ``N`` of the core nucleus
"""
function cZN_from_corenuc(rZ,rN,corenuc::String)
    if corenuc == ""
        return rZ,rN
    end
    reg = r"[0-9]+"
    A = match(reg,corenuc).match 
    el = replace(corenuc, A=>"")
    Z = 0
    for tZ=1:length(nuclist)
        if el == nuclist[tZ]
            Z = tZ;break
        end
    end 
    A = parse(Int64,A); N = A-Z   
    return Z,N
end

"""
    rm_comment(lines)
remove fortran like comment from input (snt fmt) strings
"""
function rm_comment(lines)
    nlines = []
    for line in lines
        line = strip(line)
        if length(line) > 0
            if startswith(line,"!")||startswith(line,"#")
                continue
            end
        end
        push!(nlines,line)
    end
    return nlines
end

function rm_nan(array)
    na = []
    for tmp in array
        if tmp != "";push!(na,tmp); end
    end
    return na
end

""" 
    readsnt(sntf,binfo,to)
Function to read snt file. Note that it is slightly different from `readsnt()` in ShellModel.jl.
"""
function readsnt(sntf,binfo,to) 
    Anum=binfo.nuc.Aref;hw=binfo.hw;emax_calc = binfo.emax;emax_calc = binfo.emax
    f = open(sntf,"r");tlines = readlines(f);close(f)
    lines = rm_comment(tlines)
    line = lines[1]
    lp,ln,cp,cn = map(x->parse(Int,x),rm_nan(split(line," ")))
    @assert lp ==ln "lp&ln must be the same!"
    p_sps = SingleParticleState[ ]; n_sps = SingleParticleState[ ]
    @inbounds for i = 1:lp
        ith,n,l,j,tz = map(x->parse(Int,x),rm_nan(split(lines[1+i]," "))[1:5])
        if (2*n+l <= emax_calc)
            push!(p_sps,SingleParticleState(n,l,j,tz,[0.0],[false],[false],[false]))
        end
    end
    @inbounds for i = 1:lp
        ith, n,l,j,tz = map(x->parse(Int,x),rm_nan(split(lines[1+i+ln]," "))[1:5])
        if (2*n+l <= emax_calc)
            push!(n_sps,SingleParticleState(n,l,j,tz,[0.0],[false],[false],[false]))
        end
    end
    lpn_calc = get_lpln_from_emax(emax_calc)
    idxofst = ifelse(lp > lpn_calc,lp-lpn_calc,0)
    sps,dicts1b = make_sps_and_dict_isnt2ims(p_sps,n_sps,emax_calc)
    dict_snt2ms = dicts1b.snt2ms

    nsp,zero = map(x->parse(Int,x),split(lines[1+ln+lp+1])[1:2])
    ntbme  = parse(Int,rm_nan(split(lines[1+ln+lp+1+nsp+1]," "))[1])
    dicts=[ Dict{Int64,Vector{Vector{Float64}}}() for pnrank=1:3]
    ofst = ln+lp+nsp+3 
    tls = @view lines[ofst+1:end]
    ci=cj=ck=cl=cJ="-1"
    cVjj=cVjj_2n3n=cVpp="0.0"
    @inbounds for ith = 1:ntbme
        tkey = zeros(Int64,4)
        tl = tls[ith]
        sptl = split(tl) 
        nsptl = length(sptl)
        if nsptl == 6
            ci,cj,ck,cl,cJ,cVjj = sptl
        elseif nsptl == 7
            ci,cj,ck,cl,cJ,cVjj,cVpp=sptl
        elseif nsptl == 8 
            ci,cj,ck,cl,cJ,cVjj,cVjj_2n3n,cVpp=sptl
        end
        totJ = parse(Float64,cJ); Vjj = parse(Float64,cVjj); Vjj_2n3n = parse(Float64,cVjj_2n3n); Vpp = parse(Float64,cVpp)
        tkey[1] = parse(Int64,ci)
        tkey[2] = parse(Int64,cj)
        tkey[3] = parse(Int64,ck)
        tkey[4] = parse(Int64,cl)
        if !check_truncated_abcd(tkey,lp,lpn_calc,idxofst,dict_snt2ms,to); continue;end
        nth = 2
        if tkey[1] % 2 == 1  && tkey[2] % 2 == 1; nth = 1;
        elseif tkey[3] % 2 == 0 && tkey[4] %2 == 0; nth=3;end
        phase = 1.0        
        #if nth == 2
            i,j,k,l = tkey
            oi = sps[i];oj = sps[j];ok = sps[k];ol = sps[l]
            ji = oi.j; jj = oj.j; jk = ok.j; jl = ol.j    
            phaseij = (-1)^(div(ji+jj,2)+totJ+1)
            phasekl = (-1)^(div(jk+jl,2)+totJ+1)
            flipij = ifelse(i>j,true,false)
            flipkl = ifelse(k>l,true,false)
            if flipij; tkey[1] = j; tkey[2] = i; phase *= phaseij;end
            if flipkl; tkey[3] = l; tkey[4] = k; phase *= phasekl;end
            if tkey[1] > tkey[3] || (tkey[1] == tkey[3] && tkey[2] > tkey[4])
                a,b,c,d = tkey
                tkey[1] = c; tkey[2] = d; tkey[3] = a; tkey[4] = b
            end
        #end
        tdict = dicts[nth]   
        Vjj *= phase; Vjj_2n3n *= phase; Vpp *= phase    
        V2b = Vjj + Vjj_2n3n + Vpp*hw/Anum
        nkey = get_nkey_from_abcdarr(tkey)

        Vcm = 0.0
        if !haskey(tdict,nkey)
            tdict[nkey] = [ [totJ,V2b,Vjj,Vjj_2n3n,Vpp*hw,Vcm] ]
        else
            push!(tdict[nkey],[totJ,V2b,Vjj,Vjj_2n3n,Vpp*hw,Vcm])
        end
    end
    return sps,dicts1b,dicts
end

"""
    make_sps_and_dict_isnt2ims(p_sps,n_sps,lp)

make dicts1b, snt-idx(sntidx) = 1-lp (proton) & lp+1~lp+ln (neutron), modelspace-idx(msidx) = odd(1,3,...)-> proton, even(2,4,...) -> neutron

returns:
- `dict_snt2ms`: from sntidx to msidx 
- `dict_ms2snt`: from msidx to sntidx
"""
function make_sps_and_dict_isnt2ims(p_sps,n_sps,emax_calc)
    lpn_calc = get_lpln_from_emax(emax_calc)
    dict_snt2ms = Dict{Int64,Int64}()
    dict_ms2snt = Dict{Int64,Int64}()
    sps = SingleParticleState[ ]
    msidx = 0
    @inbounds for i = 1:lpn_calc
        msidx +=1
        push!(sps,p_sps[i])
        dict_snt2ms[i] = msidx
        dict_ms2snt[msidx] = i
        msidx +=1
        push!(sps,n_sps[i])
        dict_snt2ms[i+lpn_calc] = msidx
        dict_ms2snt[msidx] = i+lpn_calc
    end
    return sps,Dict1b(dict_snt2ms,dict_ms2snt)
end
"""
    make_sps_from_pnsps(p_sps,n_sps,Chan1b)

construct `sps` from `p_sps` and `n_sps`
"""
function make_sps_from_pnsps(p_sps,n_sps,Chan1b)
    dict_snt2ms = Chan1b.snt2ms; dict_ms2snt = Chan1b.ms2snt
    lp = length(p_sps)
    sps = SingleParticleState[ ]
    hit = 0
    @inbounds for i = 1:lp
        push!(sps,p_sps[i]); hit +=1; dict_snt2ms[i] = hit;   dict_ms2snt[hit] = i
        push!(sps,n_sps[i]); hit +=1; dict_snt2ms[i+lp] = hit;dict_ms2snt[hit] = i+lp
    end
    return sps
end

function get_lpln_from_emax(emax)
    lpn = 0
    for e = 0:emax
        for n = 0:div(e,2)
            l = e-2*n
            lpn += ifelse(l==0,1,2)
        end
    end
    return lpn
end

""" 
    readsnt(sntf,binfo,to)
Function to read snt.bin file.
"""
function readsnt_bin(sntf,binfo,to;use_Float64=false) 
    Anum=binfo.nuc.Aref;hw=binfo.hw;emax_calc = binfo.emax
    f = open(sntf,"r")
    lp = read(f,Int); ln = read(f,Int)
    if lp != ln; pringln("lp&ln must be the same! err in readsnt_bin");exit();end
    cp = read(f,Int); cn = read(f,Int)
    p_sps = SingleParticleState[ ]; n_sps = SingleParticleState[ ]
    lpn_calc = get_lpln_from_emax(emax_calc)
    idxofst = ifelse(lp > lpn_calc,lp-lpn_calc,0)
    @inbounds for i = 1:lp
        ith = read(f,Int); n = read(f,Int); l = read(f,Int)
        j = read(f,Int); tz = read(f,Int)
        if 2*n + l <= emax_calc;
            push!(p_sps,SingleParticleState(n,l,j,tz,[0.0],[false],[false],[false]))
        end
    end
    @inbounds for i = 1:ln
        ith = read(f,Int); n = read(f,Int); l = read(f,Int)
        j = read(f,Int); tz = read(f,Int)
        if 2*n + l <= emax_calc;
            push!(n_sps,SingleParticleState(n,l,j,tz,[0.0],[false],[false],[false]))
        end
    end
    sps,dicts1b = make_sps_and_dict_isnt2ims(p_sps,n_sps,emax_calc)
    dict_snt2ms = dicts1b.snt2ms

    nsp = read(f,Int); zero = read(f,Int); thw = read(f,Float64)
    spes = [ [read(f,Int),read(f,Int),read(f,Float64)] for i=1:nsp]
    ntbme = read(f,Int); massop = read(f,Int); thw = read(f,Float64)
    dicts=[ Dict{Int64,Vector{Vector{Float64}}}() for pnrank=1:3]

    for n = 1:ntbme
        tkey = zeros(Int64,4)        
        org_ijkl = [read(f,Int16) for k=1:4];tkey .= org_ijkl                
        totJ = read(f,Int16)
        if use_Float64
            Vjj = read(f,Float64); Vjj_2n3n = read(f,Float64); Vpp = read(f,Float64)
        else
            Vjj = Float64(read(f,Float32)); Vjj_2n3n = Float64(read(f,Float32)); Vpp = Float64(read(f,Float32))        
        end
        if !check_truncated_abcd(tkey,lp,lpn_calc,idxofst,dict_snt2ms,to); continue;end
        nth = 2
        if tkey[1] % 2 == 1  && tkey[2] % 2 == 1; nth = 1;
        elseif tkey[3] % 2 == 0 && tkey[4] %2 == 0; nth=3;end
        phase = 1.0
        if nth == 2
            i,j,k,l = tkey
            oi = sps[i];oj = sps[j];ok = sps[k];ol = sps[l]
            ji = oi.j; jj = oj.j; jk = ok.j; jl = ol.j
            phaseij = (-1)^(div(ji+jj,2)+totJ+1); flipij = ifelse(i>j,true,false)
            phasekl = (-1)^(div(jk+jl,2)+totJ+1); flipkl = ifelse(k>l,true,false)
            if flipij; tkey[1] = j; tkey[2] = i; phase *= phaseij;end
            if flipkl; tkey[3] = l; tkey[4] = k; phase *= phasekl;end
            if tkey[1] > tkey[3] || (tkey[1]==tkey[3]&&tkey[2]>tkey[4])
                a,b,c,d = tkey
                tkey[1] = c; tkey[2] = d; tkey[3] = a; tkey[4] = b
            end
        end
        tdict = dicts[nth]
        Vjj *= phase
        Vjj_2n3n *= phase
        Vpp *= phase
        V2b = Vjj + Vjj_2n3n + Vpp*hw/Anum
        nkey = get_nkey_from_abcdarr(tkey)
        Vcm= 0.0
        tJ = totJ * 1.0
        if !haskey(tdict,nkey)
            tdict[nkey] = [ [tJ,V2b,Vjj,Vjj_2n3n,Vpp*hw,Vcm] ]
        else
            push!(tdict[nkey],[tJ,V2b,Vjj,Vjj_2n3n,Vpp*hw,Vcm])
        end
    end
    close(f)
    return sps,dicts1b,dicts
end

function check_truncated_abcd(tkey,lp,lpn_calc,idxofst,dict_snt2ms,to)
    tf = true
    for k in eachindex(tkey)
        org_sntidx = tkey[k]
        if org_sntidx > lp # neutron sntidx
            sntidx = org_sntidx - idxofst
            dicget = get(dict_snt2ms,sntidx,0)
            if dicget == 0
                return false
            end
            tkey[k] = dicget
        else #proton sntidx
            if org_sntidx <= lpn_calc
                tkey[k] = dict_snt2ms[tkey[k]]
            else
                return false
            end
        end
    end
    return tf
end

"""
    get_abcdarr_from_intkey!(nkey,abcdarr;ofst=1000)

A destructive function to get original array from integer array.
"""
function get_abcdarr_from_intkey!(nkey,abcdarr;ofst=1000)
    intkey = nkey
    q = div(intkey,ofst^3); abcdarr[4] = q;intkey -= q *ofst^3
    q = div(intkey,ofst^2); abcdarr[3] = q;intkey -= q *ofst^2
    q = div(intkey,ofst^1); abcdarr[2] = q;intkey -= q *ofst^1
    abcdarr[1] = intkey 
    return nothing
end

"""
    store_1b2b(sps,dicts1b,dicts,binfo)
"""
function store_1b2b(sps,dicts1b::Dict1b,dicts,binfo)
    Anum = binfo.nuc.Aref; hw = binfo.hw
    dim1b = div(length(sps),2)
    ### store one-body part
    p1b = zeros(Float64,dim1b,dim1b)
    n1b = zeros(Float64,dim1b,dim1b)
    massfac = hw * (Anum-1)/Anum
    for pn = 1:2
        v1b = ifelse(pn==1,p1b,n1b)
        for ith = 1:dim1b 
            idx = pn + 2*(ith-1)        
            orbit = sps[idx]
            n = orbit.n; l = orbit.l; tj = orbit.j
            temax = 2*n+l
            v1b[ith,ith] = 0.5 * hw * (temax+3/2) * (Anum-1) /Anum
            nlj_i = [n,l,tj]
            for jth = ith:dim1b
                jdx = pn + 2*(jth-1)
                if jth==ith;continue;end
                nlj_j = [sps[jdx].n,sps[jdx].l,sps[jdx].j]
                tmp = kinetic_ob(nlj_i,nlj_j) * massfac
                v1b[ith,jth] = v1b[jth,ith] = tmp 
            end
        end
    end
    ### store two-body part  dictMonopole is used in IMSRG
    dictTBMEs = [ Dict{Vector{Int64},Vector{Float64}}( ) for pnrank=1:3]
    dictMonopole = [ Dict{Vector{Int64},valDictMonopole}( ) for pnrank=1:3]
    for pnrank=1:3
        tdictl = dicts[pnrank]                
        tdict = dictTBMEs[pnrank]
        tmdict = dictMonopole[pnrank]
        for intkey in keys(tdictl)    
            tkey = zeros(Int64,4)                
            get_abcdarr_from_intkey!(intkey,tkey)
            ja = sps[tkey[1]].j; jb = sps[tkey[2]].j
            jc = sps[tkey[3]].j; jd = sps[tkey[4]].j
            if div(ja+jb,2) != div(jc+jd,2);continue;end
            ts = tdictl[intkey]
            sqfac = sqrt( (1+delta(tkey[1],tkey[2])) * (1+delta(tkey[3],tkey[4])))
            vmono,vmono2n3n = calc_vmono_for_given_abcd(ts,Anum,1.0,false)
            tdict[tkey] = [0.0,0.0]
            tdict[tkey][1] = vmono *sqfac
            tdict[tkey][2] = vmono2n3n *sqfac
            exkey = [tkey[3],tkey[4],tkey[1],tkey[2]]
            tdict[exkey] = [0.0,0.0]
            tdict[exkey][1] = vmono *sqfac
            tdict[exkey][2] = vmono2n3n *sqfac
            if tkey[1] == tkey[3] && tkey[2] == tkey[4]  
                vmval = vmono * sqfac / ( (ja+1)*(jb+1) )              
                tmdict[ [tkey[1],tkey[2]]] = valDictMonopole([vmval,0.0],Vector{Int64}[ ])
            end
            #exchange term 
            #i<->j is nontrivial for monopole components
            if tkey[1] != tkey[2] 
                exkey = copy(tkey)
                exkey[1] = tkey[2]; exkey[2] = tkey[1]
                ja = sps[tkey[1]].j; jb = sps[tkey[2]].j
                phasefac = (-1)^(div(ja+jb,2)+1)
                vmono,vmono2n3n = calc_vmono_for_given_abcd(ts,Anum,phasefac,true)
                tdict[exkey] = [0.0,0.0]; tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]] = [0.0,0.0]
                tdict[exkey][1] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][1] = vmono *sqfac 
                tdict[exkey][2] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][2] = vmono2n3n *sqfac 
            end
           
            if tkey[3] != tkey[4]
                exkey = copy(tkey)
                exkey[3] = tkey[4]; exkey[4] = tkey[3]
                jc = sps[tkey[3]].j; jd = sps[tkey[4]].j
                phasefac = (-1)^(div(jc+jd,2)+1)
                vmono,vmono2n3n = calc_vmono_for_given_abcd(ts,Anum,phasefac,true)
                tdict[exkey] = [0.0,0.0]                
                tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]] = [0.0,0.0]
                tdict[exkey][1] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][1] = vmono *sqfac 
                tdict[exkey][2] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][2] = vmono2n3n *sqfac 
            end

            if ((tkey[1] != tkey[2]) && (tkey[3] != tkey[4]))
                exkey = copy(tkey)
                exkey[1] = tkey[2]; exkey[2] = tkey[1]
                exkey[3] = tkey[4]; exkey[4] = tkey[3]    
                ja = sps[tkey[1]].j; jb = sps[tkey[2]].j
                jc = sps[tkey[3]].j; jd = sps[tkey[4]].j
                phase = (-1)^(div(jc+jd,2)+div(jc+jd,2))
                vmono,vmono2n3n = calc_vmono_for_given_abcd(ts,Anum,phase,false)
                tdict[exkey] = [0.0,0.0]
                tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]] = [0.0,0.0]
                tdict[exkey][1] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][1] = vmono *sqfac 
                tdict[exkey][2] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][2] = vmono2n3n *sqfac 
            end
        end
    end
    Chan1b = def_chan1b(dim1b,sps,dicts1b)
    Chan2bD,Gamma,maxnpq,V2 = def_chan2b(binfo,dicts,sps)
    dictsnt = dictSnt(dictTBMEs,dictMonopole)
    Hamil = Operator([0.0],[p1b,n1b],V2,[true],[false])
    return Hamil,dictsnt,Chan1b,Chan2bD,Gamma,maxnpq
end

function calc_vmono_for_given_abcd(ts,Anum,phasefac,phaseJ::Bool;calcT=true,calcCM=true)    
    vmono = vmono2n3n = 0.0            
    for t in ts
        J = t[1]
        v = 0.0
        if calcT || calcCM  
            v = t[3] + t[4] + ifelse(calcT,t[5] /Anum,0.0) + ifelse(calcCM,t[6] * Anum,0.0)
        else
            v = t[3] + t[4]
        end
        v2n3n = t[4]
        vmono += (2*J+1) * v * phasefac * ifelse(phaseJ,(-1)^J,1.0)
        vmono2n3n += (2*J+1) * v2n3n * phasefac * ifelse(phaseJ,(-1)^J,1.0)
    end
    return vmono,vmono2n3n
end

function get_pn_sps(sps)
    dim1b = div(length(sps),2)
    p_sps = sps[1:2:2*dim1b-1]
    n_sps = sps[2:2:2*dim1b]
    return p_sps,n_sps
end 

"""
    def_chan1b(dim1b,sps,dicts1b)

define Chan1b: dict. to get 1b-channels to be coupled to a given channel
Chan1b = [ dict_for_proton_sps, dict_for_neutron_sps] dicts1b
Basically ToBeCoupled will be used, but redundant one is needed in some cases
"""
function def_chan1b(dim1b,sps,dicts1b)
    dict_snt2ms = dicts1b.snt2ms
    dict_ms2snt = dicts1b.ms2snt 
    ToBeCoupled = [Dict{Int64,Vector{Int64}}() for pn = 1:2]
    AllTBC = [Dict{Int64,Vector{Int64}}() for pn = 1:2]
    for pn = 1:2
        target = ToBeCoupled[pn]
        target2 = AllTBC[pn]
        for ii = 1:dim1b
            idx = pn + 2*(ii-1)
            for jj = ii:dim1b
                jdx = pn + 2*(jj-1)
                ljtz1 = [ sps[idx].l, sps[idx].j, sps[idx].tz]
                ljtz2 = [ sps[jdx].l, sps[jdx].j, sps[jdx].tz]
                if ljtz1 == ljtz2
                    t = get(target,idx,-100)
                    if t == -100
                        target[idx] = [jdx]
                    else
                        push!(target[idx],jdx)
                    end
                end
            end
            for jj = 1:dim1b
                jdx = pn + 2*(jj-1)
                ljtz1 = [ sps[idx].l, sps[idx].j, sps[idx].tz]
                ljtz2 = [ sps[jdx].l, sps[jdx].j, sps[jdx].tz]
                if ljtz1 == ljtz2
                    t = get(target2,idx,-100)
                    if t == -100
                        target2[idx] = [jdx]
                    else
                        push!(target2[idx],jdx)
                    end
                end
            end
        end
    end
    return chan1b(ToBeCoupled,AllTBC,dict_snt2ms,dict_ms2snt)
end

"""
    def_chan2b(binfo,dicts,sps)
define two-body utils and returns them as `Chan2bD` struct 
"""
function def_chan2b(binfo,dicts,sps)
    Anum = binfo.nuc.Aref; emax = binfo.emax
    Jmax = 2*emax+1
    dim1b = div(length(sps),2)
    nchan = 0
    Chan2b = chan2b[ ]
    maxnpq = 0
    #dict_2b_ch = Dict{Vector{Int64},VdictCh}()
    dict_2b_ch = Dict{UInt64,VdictCh}()
    Gamma = Matrix{Float64}[ ]
    V2b = Matrix{Float64}[ ]
    tkey = zeros(Int64,4); dkey = zeros(Int64,3)
    for (tidx,Tz) in enumerate(collect(-2:2:2)) ## pp/pn/nn
        tdict = dicts[tidx]
        pn1 = ifelse(Tz==2,2,1) 
        pn2 = ifelse(Tz==-2,1,2)
        for prty = 1:-2:-1        
            for J = 0:Jmax
                kets = Vector{Int64}[ ]
                vdict = Dict{Int64,Int64}()
                for idx_a = 1:dim1b
                    a = pn1 + 2*(idx_a-1)
                    la = sps[a].l; ja = sps[a].j
                    bmin = ifelse(Tz==0,1,idx_a)
                    for idx_b=1:dim1b
                        b = pn2 + 2*(idx_b-1)
                        lb = sps[b].l; jb = sps[b].j
                        if tri_check(ja//2,jb//2,J)==false;continue;end
                        if abs(Tz)==2 && a == b && J % 2 == 1;continue;end                        
                        tprty = (-1)^(la+lb)
                        if prty != tprty; continue;end
                        ta = a; tb = b
                        if a > b;ta=b;tb=a;end
                        if !([ta,tb] in kets)
                            push!(kets,[ta,tb])
                            nkey_ab = get_nkey2(ta,tb)
                            vdict[nkey_ab] = length(kets)
                        end
                    end
                end
                nchan += 1
                kets = sort(kets)
                for (ik,ket) in enumerate(kets)
                    intkey = get_nkey2(ket[1],ket[2])
                    vdict[intkey] = ik
                end

                tuplekets = [ (ket[1],ket[2]) for ket in kets]
                push!(Chan2b, chan2b(Tz,prty,J,tuplekets,length(kets)))
                dim = length(kets)
       
                # org
                #push!(Chan2b, chan2b(Tz,prty,J,kets,length(kets)))
                #dim = length(kets)
                push!(Gamma,zeros(Float64,dim,dim))

                dkey[1] = Tz; dkey[2]=prty; dkey[3]=J
                maxnpq = ifelse(dim>maxnpq,dim,maxnpq)
                #dict_2b_ch[copy(dkey)] = VdictCh(nchan,vdict)
                dict_2b_ch[get_nkey3_JPT(copy(dkey))] = VdictCh(nchan,vdict)
                
                vmat = zeros(Float64,dim,dim)
                for i = 1:dim
                    a,b = kets[i]
                    tkey[1] = a; tkey[2] = b
                    for j=i:dim
                        c,d = kets[j]
                        tkey[1] = a; tkey[2] = b
                        tkey[3] = c; tkey[4] = d
                        if tkey[1] > tkey[3] 
                            tkey[1]=c; tkey[2]=d
                            tkey[3]=a; tkey[4]=b
                        end
                        intkey = get_nkey_from_abcdarr(tkey)
                        #println("tkey $tkey intkey $intkey")
                        for JV in tdict[intkey]
                            tJ = JV[1]
                            v = JV[3] + JV[4] + JV[5] /Anum  + JV[6] * Anum                           
                            if Int(tJ) != J;continue;end                            
                            vmat[i,j] = vmat[j,i] = v
                        end
                    end
                end
                push!(V2b, vmat)
            end
        end
    end
    nch = length(Chan2b)
    dict_ch_idx_from_ket = [ Dict{UInt64,NTuple{2,Int64}}( ) for Tz = -2:2:2] # Tz vector
    dict_idx_from_chket = [Dict{Vector{Int64},Int64}( ) for _ = 1:nch]
    for ch = 1:nch
        tbc = Chan2b[ch]; tkets = tbc.kets; Tz = tbc.Tz; J = tbc.J
        target = dict_ch_idx_from_ket[div(Tz+2,2)+1]
        target2 = dict_idx_from_chket[ch]
        for (idx,ket) in enumerate(tkets)
            nkey = get_nkey3_ketJ(ket[1],ket[2],J)
            target[nkey] = (ch,idx)
            target2[[ket[1],ket[2]]] = idx
        end
    end
    return chan2bD(Chan2b,dict_2b_ch,dict_ch_idx_from_ket,dict_idx_from_chket),Gamma,maxnpq,V2b
end

"""
    update_1b!(binfo,sps,Hamil)
Update one-body(1b) part of Hamiltonian for different target nuclei
"""
function update_1b!(binfo,sps,Hamil::Operator)
    p1b = Hamil.onebody[1]; n1b = Hamil.onebody[2]
    Anum = binfo.nuc.Aref; hw = binfo.hw
    ### store one-body part
    p_sps,n_sps = get_pn_sps(sps)
    lp =length(p_sps); ln =length(n_sps)       
    massfac = hw * (Anum-1)/Anum
    for pn = 1:2
        tsps = ifelse(pn==1,p_sps,n_sps)
        v1b = ifelse(pn==1,p1b,n1b)
        lmax = ifelse(pn==1,lp,ln)
        v1b .= 0.0
        for i = 1:lmax
            n = tsps[i].n; l = tsps[i].l; j = tsps[i].j
            temax = 2*n+l
            v1b[i,i] = 0.5 * hw * (temax+3/2) * (Anum-1) /Anum
            nlj_i = [n,l,j]
            for jj = i:lmax
                if jj==i;continue;end
                nlj_j = [ tsps[jj].n, tsps[jj].l, tsps[jj].j]
                tmp = kinetic_ob(nlj_i,nlj_j) * massfac
                v1b[i,jj] = v1b[jj,i] = tmp 
            end
        end
    end
    return nothing
end

"""
    update_2b!(binfo,sps,Hamil,dictTBMEs,Chan2bD,dicts)

Update two-body(2b) kinetic part for different target nuclei
"""
function update_2b!(binfo,sps,Hamil,dictTBMEs,Chan2bD,dicts;Hcm=nothing)
    emax = binfo.emax; A = binfo.nuc.A
    V2 = Hamil.twobody
    Chan2b = Chan2bD.Chan2b
    tkey = zeros(Float64,4)
    for pnrank=1:3
        tdictl = dicts[pnrank]        
        tdict = dictTBMEs[pnrank]
        for intkey in keys(tdictl)
            tkey = zeros(Int64,4)
            get_abcdarr_from_intkey!(intkey,tkey)
            ja = sps[tkey[1]].j
            jb = sps[tkey[2]].j
            jc = sps[tkey[3]].j
            jd = sps[tkey[4]].j
            if div(ja+jb,2) != div(jc+jd,2);continue;end
            ts = tdictl[intkey]
            sqfac = sqrt( (1+delta(tkey[1],tkey[2])) * (1+delta(tkey[3],tkey[4])))
            vmono,vmono2n3n = calc_vmono_for_given_abcd(ts,A,1.0,false)
            tdict[tkey][1] = sqfac * vmono
            tdict[tkey][2] = sqfac * vmono2n3n
            exkey = [tkey[3],tkey[4],tkey[1],tkey[2]]
            tdict[exkey][1] = vmono *sqfac
            tdict[exkey][2] = vmono2n3n *sqfac
            if tkey[1] != tkey[2]
                exkey = copy(tkey)
                exkey .= tkey                
                exkey[1] = tkey[2]; exkey[2] = tkey[1]
                ja = sps[tkey[1]].j ; jb = sps[tkey[2]].j
                phasefac = (-1)^(div(ja+jb,2)+1)
                vmono,vmono2n3n = calc_vmono_for_given_abcd(ts,A,phasefac,true)
                tdict[exkey][1] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][1] = vmono *sqfac 
                tdict[exkey][2] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][2] = vmono2n3n *sqfac 
            end
            if tkey[3] != tkey[4]
                exkey = copy(tkey)
                exkey[3] = tkey[4]; exkey[4] = tkey[3]
                jc = sps[tkey[3]].j; jd = sps[tkey[4]].j
                phasefac = (-1)^(div(jc+jd,2)+1)
                vmono,vmono2n3n = calc_vmono_for_given_abcd(ts,A,phasefac,true)
                tdict[exkey][1] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][1] = vmono *sqfac
                tdict[exkey][2] = tdict[[exkey[3],exkey[4],exkey[1],exkey[2]]][2] = vmono2n3n *sqfac 
            end
        end
    end
    
    Jmax = 2*emax+1
    nchan = 0
    tkey = [0,0,0,0]
    dkey = [0,0,0]
    for (tidx,Tz) in enumerate(collect(-2:2:2)) ## pp/pn/nn
        tdict = dicts[tidx]
        for prty = 1:-2:-1        
            for J = 0:Jmax
                nchan += 1
                kets = Chan2b[nchan].kets
                dim = length(kets)
                dkey[1] = Tz; dkey[2]=prty; dkey[3]=J
                vmat = V2[nchan]
                for i = 1:dim
                    a,b = kets[i]                
                    for j=i:dim
                        c,d = kets[j]
                        tkey[1] = a; tkey[2] = b
                        tkey[3] = c; tkey[4] = d
                        if tkey[1] > tkey[3]
                            tkey[1] = c; tkey[2]=d
                            tkey[3] = a; tkey[4]=b
                        end
                        intkey = get_nkey_from_abcdarr(tkey)
                        for JV in tdict[intkey]
                            tJ = JV[1]
                            v = JV[3] + JV[4] + JV[5] /A  + JV[6] * A
                            if Int(tJ) != J;continue;end
                            vmat[i,j] = vmat[j,i] = v
                        end
                    end
                end         
            end
        end
    end
    return nothing 
end

## for debug
function print_rho(rho_p,rho_n)
    println("rho_p & rho_n")
    for i=1:size(rho_p)[1]; print_vec("", @views rho_p[i,:]);end
    for i=1:size(rho_n)[1]; print_vec("", @views rho_n[i,:]);end
    return nothing
end
function print_F(h_p,h_n)
    println("Fock Mat.")
    lp = size(h_p)[1];ln = size(h_n)[1]
    for i=1:lp;print_vec("", @views h_p[i,:]);end
    for i=1:ln;print_vec("", @views h_n[i,:]);end
    return nothing
end
function print_V2b(h_p,p1b,h_n,n1b)
    println("V2b p/n")
    lp = size(h_p)[1];ln = size(h_n)[1]
    for i=1:lp;print_vec("", @views h_p[i,:] - @views p1b[i,:]);end
    for i=1:ln;print_vec("", @views h_n[i,:] - @views n1b[i,:]);end
    return nothing
end
