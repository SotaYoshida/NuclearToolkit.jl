const element = ["H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
                 "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
                 "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
                 "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
                 "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
                 "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

const reg = r"[0-9]+"

struct bit2b
    a::Int64
    b::Int64
    c::Int64
    d::Int64
end

struct ifph
    i::Int64
    f::Int64
    phase::Bool
end

struct Jpninfo
    fac::Float64
    pjump::Array{ifph,1}
    njump::Array{ifph,1}
end

"""
Struct to store 2b-jumps (T=1).

This includes the following information:
`f::Int64` the final index `f` for a given index `i`/
`coef::Float64` the coefficient of the 2b-jump.
"""
struct T1info
    f::Int64
    coef::Float64
end

struct MiMf
    Mi::Int64
    Mf::Int64
    fac::Float64
end

"""
main_sm(sntf,target_nuc,num_ev,target_J;
        save_wav=false,q=1,is_block=false,is_show=false,num_history=3,lm=300,ls=20,tol=1.e-8,
        in_wf="",mdimmode=false,calc_moment=false, visualize_occ=false, gfactors=[1.0,0.0,5.586,-3.826],effcharge=[1.5,0.5])

Digonalize the model-space Hamiltonian 

# Arguments: 

- `sntf`:       path to input interaction file (.snt fmt)
- `target_nuc`: target nucleus
- `num_ev`:     number of eigenstates to be evaluated
- `target_J`:   target total J (specified by e.g. [0]). Set to [] if you want lowest states with any J. 
   Note that J should be doubled (J=0=>[0], J=1/2=>[1], J=1=>[2],...) 

# Optional arguments:
- `q=1`              block size for Block-Lanczos methods 
- `is_block=false`   whether or not to use Block algorithm 
- `save_wav=false`   whether or not to save wavefunction file 
- `is_show = true`   to show elapsed time & allocations 
- `lm = 100`         number of Lanczos vectors to store 
- `ls = 20`          number of vectors to be used for Thick-Restart 
- `tol= 1.e-8`       tolerance for convergence check in the Lanczos method 
- `in_wf=""`      path to initial w.f. (for preprocessing) 
- `mdimmode=false`   `true` => calculate only the M-scheme dimension
- `calc_moment=false`  `true` => calculate mu&Q moments 
- `calc_entropy=false` `true` => calculate entropy, which is not optimized yet.
- `visualize_occ=false` `true` => visualize all configurations to be considered
- `gfactors=[1.0,0.0,5.586,-3.826]` angular momentum and spin g-factors 
- `effcgarge=[1.5,0.5]` effective charges 
- `truncation_scheme==""` option to specify a truncation scheme, "jocc" and "pn-pair" is supported and you can use multiple schemes by separating them with camma like "jocc,pn-pair".
- `truncated_jocc=Dict{String,Vector{Int64}}()` option to specify the truncation scheme for each orbit, e.g. Dict("p0p1"=>[1],"n0p3"=>[2,3]) means that the occupation number for proton 0p1 is truncated upto 1, and for neutron 0p3 min=2 and max=3"
- `debugmode`= "" option to specify debug mode. Lanczos/BlockLanczos are supported
"""
function main_sm(sntf,target_nuc,num_ev,target_J;save_wav=false,q=1,is_block=false,is_show=false,num_history=3,lm=300,ls=20,tol=1.e-8,
                 in_wf="",mdimmode=false, print_evec=false, calc_moment = false, calc_entropy = false, visualize_occ = false,
                 gfactors = [1.0,0.0,5.586,-3.826], effcharge=[1.5,0.5], truncation_scheme="",truncated_jocc=Dict{String,Vector{Int64}}(),
                 debugmode="")

    to = TimerOutput()
    Anum = parse(Int64, match(reg,target_nuc).match)
    lp,ln,cp,cn,massop,Aref,pow,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsmsnt(sntf,Anum)
    massformula = ifelse(16<=Anum<= 40,2,1)
    ## massformula=2: J. Blomqvist and A. Molinari, Nucl. Phys. A106, 545 (1968).
    ## we use this for the sd-shell nuclei
    hw, bpar = init_ho_by_mass(Anum,massformula)

    @assert length(target_J) <= 1 "Invalid `target_J` Multiple J is not supported now."
    Mtot = ifelse(Anum%2==0,0,1)
    tJ = -1; eval_jj = -1.0
    if length(target_J) > 0; Mtot = minimum(target_J);tJ=target_J[1]
        eval_jj = 0.5*tJ*(tJ/2+1)
    end
    if Anum % 2 != Mtot % 2; println("invalid targetJ $tJ");exit();end
    target_el = replace(target_nuc, string(Anum)=>"")
    Z,N,vp,vn = getZNA(target_el,Anum,cp,cn)
    tf_pairwise_truncation = ifelse(occursin("-pair",truncation_scheme),true,false)
    if tf_pairwise_truncation
        #if Z!=N; @error "Invalid truncation scheme: pn-pair is not supported for Z != N nuclei"; exit();end
        if is_block; @error "Block Lanczos method is not supported for pn-pair truncation"; exit();end
    end
    msps_p,msps_n,mz_p,mz_n = construct_msps(p_sps,n_sps)
    pbits,nbits,jocc_p,jocc_n,Mps,Mns,tdims = occ(p_sps,msps_p,mz_p,vp,n_sps,msps_n,mz_n,vn,Mtot;truncation_scheme=truncation_scheme,truncated_jocc=truncated_jocc)

    lblock=length(pbits)
    mdim = tdims[end]; if mdim==0; println("Aborted due to the mdim=0");return true;end
    num_ev = min(num_ev,mdim)
    mdim_print(target_nuc,Z,N,cp,cn,vp,vn,mdim,tJ)
    if mdimmode; return nothing;end
    if visualize_occ; visualize_configurations(msps_p,msps_n,pbits,nbits,mdim); end
    if tf_pairwise_truncation && mdim >= 10^4
        @warn "For pairwise truncation, constructing H explicitly, is not recommended for large model space."
    end

    ##Making bit representation of Hamiltonian operators, and storing them
    @timeit to "prep. 1b&2b-jumps" begin
        ## storing two-body jumps for pp/nn 2b interaction
        pp_2bjump,nn_2bjump = prep_Hamil_T1(p_sps,n_sps,msps_p,msps_n,pbits,nbits,labels,TBMEs)        
        bVpn,Vpn,delMs = prep_Hamil_pn(p_sps,n_sps,msps_p,msps_n,labels[3],TBMEs[3])

        ## storing one-body jumps for pn 2b interaction
        bis,bfs,p_NiNfs,n_NiNfs,num_task = prep_pn(lblock,pbits,nbits,Mps,delMs,bVpn,Vpn)
        bVpn=nothing
        ## distribute task
        block_tasks = make_distribute(num_task)
    end
   
    @timeit to "Jcalc." begin
        oPP,oNN,oPNu,oPNd = prep_J(tdims,p_sps,n_sps,msps_p,msps_n,pbits,nbits)
        Js = [ 0.5*Mtot*(0.5*Mtot+1) for i = 1:num_ev]
        Jtasks = zeros(Int64,lblock)
        for i = 1:lblock
            Jtasks[i] = length(pbits[i])*length(nbits[i])
        end
        Jidxs = make_distribute_J(Jtasks)
    end

    @timeit to "Lanczos" begin
        en =[ [1.e+4 for j=1:num_ev] for i = 1:num_history]
        Rvecs = [ zeros(Float64,mdim) for i=1:num_ev]
        Tmat = zeros(Float64,lm,lm)
        vks = Vector{Float64}[]; uks=Vector{Float64}[];itmin = 1; elit=1
        doubleLanczos = false
        if tJ !=-1; doubleLanczos = true;end
            
        if is_block #Thick-Restart double Block Lanczos: TRdBL
            ls_sub = div(ls,q)
            vks = [ zeros(Float64,q,mdim) for i=1:div(lm,q)]
            uks = [ zeros(Float64,q,mdim) for i=1:ls_sub*q]
            Beta_H = zeros(Float64,q,q)
            V = vks[1]
            if in_wf !=""
                try
                    read_appwav(in_wf,mdim,V,q,true,is_block)
                    bl_QR!(V',Beta_H,mdim,q)
                catch
                    println("error @preprocessing: failed to read appwav")
                    initialize_bl_wav(mdim,q,vks[1])
                    bl_QR!(V',Beta_H,mdim,q)
                end
            else
                initialize_bl_wav(mdim,q,vks[1])
                bl_QR!(V',Beta_H,mdim,q)
            end
            elit = TRBL(q,vks,uks,Tmat,Beta_H,pbits,nbits,jocc_p,jocc_n,SPEs,
                         pp_2bjump,nn_2bjump,bis,bfs,block_tasks,
                         p_NiNfs,n_NiNfs,Mps,delMs,Vpn,tdims,
                         eval_jj,oPP,oNN,oPNu,oPNd,Jidxs,
                         num_ev,num_history,lm,ls_sub,en,tol,to;doubleLanczos=doubleLanczos,debugmode=debugmode)
        else #Thick Restart (double) Lanczos: TRL
            vks = [ zeros(Float64,mdim) for i=1:lm]
            uks = [ zeros(Float64,mdim) for i=1:ls]
            if in_wf==""
                initialize_wf(vks[1],"rand",tJ,mdim)
            else
                read_appwav(in_wf,mdim,vks[1],1,true,is_block)
            end
            elit = TRL(vks,uks,Tmat,itmin,
                       pbits,nbits,jocc_p,jocc_n,SPEs,
                       pp_2bjump,nn_2bjump,bis,bfs,block_tasks,
                       p_NiNfs,n_NiNfs,Mps,delMs,Vpn,
                       eval_jj,oPP,oNN,oPNu,oPNd,Jidxs,
                       tdims,num_ev,num_history,lm,ls,en,tol,to;doubleLanczos=doubleLanczos,debugmode=debugmode)
        end
    end

    @timeit to "Rvecs" begin
        vals,vecs = eigen(@views Tmat[1:elit*q,1:elit*q])
        @inbounds for (nth,Rvec) in enumerate(Rvecs)
            if is_block == false
                for k=1:length(vals)
                    Rvec .+= vecs[k,nth] .* vks[k]
                end
            else
                for k=1:length(vals)
                    it = div(k-1,q)
                    b = k - q*it 
                    Rvec .+= @views vecs[k,nth] .* vks[it+1][b,:]
                end
            end
            Rvec .*= 1.0/sqrt(dot(Rvec,Rvec))
        end        
    end

    vt = zeros(Float64,mdim)
    for (nth,Rv) in enumerate(Rvecs)
        vt .= 0.0
        operate_J!(Rv,vt,pbits,nbits,tdims,Jidxs,oPP,oNN,oPNu,oPNd)
        Js[nth] += dot(Rv,vt)
    end

    if print_evec
        print("\n")
        if target_nuc == "Li6"
            svd_li6(Rvecs)
        end
    end
    
    totJs = J_from_JJ1.(Js)
    #println("totJs $totJs")
    tx_mom =""
    if calc_moment 
        tx_mom = eval_moment(Mtot,Rvecs,totJs,p_sps,n_sps,msps_p,msps_n,tdims,
                             jocc_p,jocc_n,pbits,nbits,bpar,gfactors,effcharge)
    end        
    if save_wav
        @timeit to "I/O" begin
            csnt = split(split(sntf,"/")[end],".")[1]
            oupf="./"*target_nuc*"_"*csnt*".wav"
            if tJ != -1;oupf="./"*target_nuc*"_"*csnt*"_j"*string(tJ)*".wav";end
            writeRitzvecs(mdim,Mtot,en[1],totJs,Rvecs,oupf)
        end
    end
    #println("sntf: $sntf")
    if length(target_J) == 0
        println("J $totJs")
        print("En. ");map(x -> @printf("%10.4f ",x), en[1])
        print("\nEx. ")
        map(x -> @printf("%10.4f ",x),[en[1][i]-en[1][1] for i=1:num_ev])
    else
        tJ = Int(2*totJs[1])
        print("2J= $tJ  En.")
        map(x -> @printf("%12.5f ",x), en[1])
    end
    print("\n")
    if tx_mom != ""
        println(tx_mom)
    end

    if tf_pairwise_truncation
        main_pairwise_truncation(truncation_scheme,Rvecs,vals[1:num_ev],tdims,msps_p,msps_n,pbits,nbits,jocc_p,jocc_n,SPEs,pp_2bjump,nn_2bjump,bis,bfs,block_tasks,p_NiNfs,n_NiNfs,Mps,delMs,Vpn,Jidxs,oPP,oNN,oPNu,oPNd)
    end

    if calc_entropy
        entropy(Rvecs,pbits,nbits,tdims,to)
    end

    show_TimerOutput_results(to;tf=is_show)
    return en[1]
end

"""
    readsmsnt(sntf,Anum)

To read interaction file in ".snt" format.
- `sntf`: path to the interaction file
- `Anum`: mass number (used for "scaling" of TBMEs)

!!! note
    The current version supports only "properly ordered",like ``a \\leq b,c \\leq d,a \\leq c`` for ``V(abcd;J)``, snt file.
"""
function readsmsnt(sntf,Anum;ignore_massop=false) 
    f = open(sntf,"r");tlines = readlines(f);close(f)
    lines = rm_comment(tlines)
    line = lines[1]    
    lp,ln,cp,cn = map(x->parse(Int,x),rm_nan(split(line," ")))
    p_sps = [zeros(Int64,4) for i=1:lp]
    n_sps = [zeros(Int64,4) for i=1:ln]
    for i = 1:lp
        ith,n,l,j,tz = map(x->parse(Int64,x),rm_nan(split(lines[1+i]," "))[1:5])
        tp = p_sps[ith]; tp[1] = n; tp[2] = l; tp[3] = j; tp[4] = tz
    end
    for i = 1:ln
        ith,n,l,j,tz = map(x->parse(Int64,x),rm_nan(split(lines[1+i+ln]," "))[1:5])
        tp = n_sps[ith-lp] 
        tp[1] = n; tp[2] = l; tp[3] = j; tp[4] = tz
    end
    nsp,zero = map(x->parse(Int,x),rm_nan(split(lines[1+ln+lp+1]," "))[1:2])
    SPEs = [ [0.0 for i=1:lp],[0.0 for i=1:ln]]
    for i = 1:nsp
        ttxt = rm_nan(split(lines[1+ln+lp+1+i]," "))
        j = parse(Int64,ttxt[1])
        idx = ifelse(j<=lp,1,2)
        if idx ==2; j -= lp; end
        SPEs[idx][j] =parse(Float64,ttxt[3])
    end
    p = 0.0
    tmp = rm_nan(split(lines[1+ln+lp+1+nsp+1]," "))
    ntbme,massop,Aref = tmp[1:3]
    if length(tmp)>3; p = tmp[4]; p = parse(Float64,p); end
    ntbme = parse(Int,ntbme); massop=parse(Int,massop)
    Aref = parse(Float64,string(Aref))
    labels = [ [ [0,0] ] for i=1:3]
    olabels = [ [0,0] ]
    for i=1:3; deleteat!(labels[i],1);end
    deleteat!(olabels,1)
    TBMEs=[ Float64[] for i =1:3] #pp/nn/pn
    oTBMEs= Float64[]
    @inbounds for ith = 1:ntbme
        tmp = rm_nan(split(lines[1+ln+lp+1+nsp+1+ith], " "))
        i = tmp[1]; j = tmp[2]; k = tmp[3]; l =tmp[4]; totJ = tmp[5]; TBME= tmp[6]
        i = parse(Int,i);j = parse(Int,j);k = parse(Int,k);l = parse(Int,l);
        nth = 0
        if i<=lp && j<=lp
            nth = 1
        elseif i>lp && j > lp
            nth = 2
        elseif i<=lp && j>lp
            nth = 3
        else
            println("i $i j $j lp $lp")
            println("err");exit()
        end
        ## snt file must be "ordered"; a<=b & c=d & a<=c
        TBME = parse(Float64,TBME)
        if massop == 1 && !ignore_massop; TBME*= (Anum/Aref)^(p);end
        if unique([i,j]) != unique([k,l])
            push!(labels[nth],[k,l,i,j,parse(Int,totJ),ith])
            push!(TBMEs[nth],TBME)
        end
        push!(labels[nth],[i,j,k,l,parse(Int,totJ),ith])
        push!(TBMEs[nth],TBME)

        push!(olabels,[i,j,k,l,parse(Int,totJ),ith])
        push!(oTBMEs,TBME)
    end
    return lp,ln,cp,cn,massop,Aref,p,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs
end

"""
   construct_msps(p_sps,n_sps)

Function to define the single particle states specified by `[n,l,j,tz,mz,p(n)idx]`.
The last elements `pidx` and `nidx` to represent original index of `p_sps`/`n_sps`, which is specified by `[n,l,j,tz]``.
"""
function construct_msps(p_sps,n_sps)
    msps_p = Vector{Int64}[] ; msps_n = Vector{Int64}[]; mz_p = Int64[]; mz_n = Int64[]
    for (pidx,tsps) in enumerate(p_sps)
        n,l,j,tz = tsps
        for mz = -j:2:j;push!(msps_p,[n,l,j,tz,mz,pidx]);push!(mz_p,mz);end
    end
    for (nidx,tsps) in enumerate(n_sps)
        n,l,j,tz = tsps
        for mz = -j:2:j;push!(msps_n,[n,l,j,tz,mz,nidx]);push!(mz_n,mz);end
    end
    return msps_p, msps_n,mz_p,mz_n
end

"""
    all_perm!(ln::Int64,num_valence::Int64,occs::Array{Array{Bool,1}})

make all possible permutation of 'bits'

Example:
If 2 protons and 1 neutron are in the 0p-shell space,
valence orbits(0p1/2,0p3/2) => -1/2, 1/2, -3/2, -1/2, 1/2, 3/2

configurations are represented like:

proton: 000011, 000101, ..., 110000

neutron: 000001, 000010, ..., 100000
"""
function all_perm!(ln::Int64, num_valence::Int64, occs::Array{Array{Bool,1}})
    for (i,tcomb) in enumerate(collect(combinations(collect(1:ln),num_valence)))
        for nth in tcomb
            occs[i][nth] = 1
        end
    end
    return nothing
end

function Mcount!(ln::Int64,mzs::Array{Int64,1}, occ::Array{Bool,1}, Mret::Array{Int64,1})
    Mret[1] = 0
    @inbounds for i = 1:ln
        if occ[i]
            Mret[1] += mzs[i]
        end
    end
    return nothing
end

"""
Function to check if the specified configuration `tocc_j` is valid within a given truncation scheme.
"""
function check_jocc_truncation(tz,m_sps,tocc_j,truncated_jocc,truncation_scheme="")
    if truncation_scheme == ""
        return true
    end
    if occursin("jocc",truncation_scheme)
        if typeof(truncated_jocc) == Dict{String,Vector{Int64}}
            if length(keys(truncated_jocc)) == 0 
                @error "You should specify `truncated_jocc::Dict{String,Vector{Int64}}`, see the document of `main_sm` (main API of ShellModel) for more details."
                exit()
            end
            return jocc_truncation_dict(tz,m_sps,tocc_j,truncated_jocc)        
        elseif typeof(truncated_jocc) == Vector{Vector{Int64}}
            if length(truncated_jocc) == 0 
                @error "Inappropriate specification of `truncated_jocc::Vector{Vector{Int64}`! You must specify it like [ [n,l,j,tz,max_occ],... ] or [ [n,l,j,tz,min_occ,max_occ],... ]."
                exit()
            end
            return jocc_truncation(tz,m_sps,tocc_j,truncated_jocc)
        else
            @error "Invalid type of truncated_jocc !!! $(typeof(truncated_jocc)) is not supported!!"   
            exit()     
        end
    elseif occursin("pn-pair",truncation_scheme) # it doesn't care pn-pair for now
        return true
    elseif occursin("nn-pair",truncation_scheme) 
        return true
    else
        @error "truncation_scheme $truncation_scheme is not defined now!!"
        exit()
    end 
    return false
end

function jocc_truncation_dict(target_tz,m_sps,tocc_j,truncated_jocc_dict::Dict{String,Vector{Int64}})
    regex = r"(\d+)(\D)(\d+)"  
    for tkey in keys(truncated_jocc_dict)
        if length(tkey) < 4
            @error "KeyError: length(tkey=$tkey) < 4,  key must include proton/neutron, n, l, j such as p0p1, p1s1, n0d5"
            exit()
        end
        if string(tkey[1]) !="p" && string(tkey[1]) != "n"
            @error "KeyError: ($tkey) key must include proton/neutron, n, l, j such as p0p1, p1s1, n0d5, etc."
            exit()
        end       
        tz = ifelse(string(tkey[1])=="p",-1,1)
        if tz != target_tz; continue;end

        n = l = j = -1
        tmp = match(regex, tkey)    
        n = parse(Int64,tmp.captures[1])
        j = parse(Int64,tmp.captures[3])
        lstr = tmp.captures[2]
        l = findall(x->x==lstr, chara_l)[1] - 1 
        if n == -1 || l == -1 || j == -1        
            println("something wrong! key $tkey => n $n l $l j $j")
            exit()
        end
        min_occ = 0; max_occ = truncated_jocc_dict[tkey][end]
        if length(truncated_jocc_dict[tkey]) == 2
            min_occ = truncated_jocc_dict[tkey][1]
        end
        occupation_number = 0
        for (idx, occ_num) in enumerate(tocc_j)
            tn, tl, tj, ttz, tjz, idx_orb = m_sps[idx]
            if tn != n || tl != l || tj != j || ttz != tz
                continue
            end
            if occ_num #bool
                occupation_number += 1
            end
        end
        if occupation_number < min_occ || occupation_number > max_occ
            return false
        end
    end
    return true
end

function jocc_truncation(target_tz,m_sps,tocc_j,truncated_jocc)
    for tkey in truncated_jocc
        ln = length(tkey)
        if ln < 5 || ln > 6
            @error "KeyError: key must be either [n,l,j,tz,occ_max] or [n,l,j,tz,occ_min,occ_max]"
            exit()
        end
        n = l = j = -1; min_occ = max_occ = 0

        if ln == 5 
            n,l,j,tz,max_occ = tkey
        elseif ln == 6
            n,l,j,tz,min_occ,max_occ = tkey
        end
        if tz != target_tz; continue;end

        occupation_number = 0
        for (idx, occ_bool) in enumerate(tocc_j)
            tn, tl, tj, ttz, tjz, idx_orb = m_sps[idx]
            if tn != n || tl != l || tj != j || ttz != tz
                continue
            end
            if occ_bool
                occupation_number += 1
            end
        end
        if occupation_number < min_occ || occupation_number > max_occ
            return false
        end
    end
    return true
end

"""
Function to check occupations in proton and neutron are only in time-reversal pairs.
"""
function check_pnpair_truncation(msps_p,msps_n,ovec_p,ovec_n,truncation_scheme="")
    if !occursin("pn-pair",truncation_scheme); return true; end
    tf = true
    for (idx_p,occ_p) in enumerate(ovec_p)
        n_p,l_p,j_p,tz_p,mz_p = msps_p[idx_p][1:5]
        for (idx_n,occ_n) in enumerate(ovec_n)
            n_n,l_n,j_n,tz_n,mz_n = msps_n[idx_n][1:5]
            if n_p != n_n || l_p != l_n || j_p != j_n || tz_p != -tz_n || mz_p != -mz_n
                continue
            end
            tf *= ifelse(occ_p == occ_n,true,false)
        end
    end
    return tf
end

function possible_mz(nljtz,mstates)
    n,l,j,tz = nljtz
    mzs = Int64[]; midxs=Int64[]
    for mz = -j:2:j
        push!(mzs,mz)
        for k = 1:length(mstates)
            if @views mstates[k][1:5] == [n,l,j,tz,mz]
                push!(midxs,k)
                break
            end
        end
    end
    return j,mzs, midxs
end

function initialize_tvec!(tvec::Array{Bool,1})
    tvec .= false
    return nothing
end

"""
    function prep_Hamil_T1(p_sps::Array{Array{Int64,1}},n_sps::Array{Array{Int64,1}},
                msps_p::Array{Array{Int64,1},1},msps_n::Array{Array{Int64,1},1},
                labels::Array{Array{Array{Int64,1},1},1},TBMEs::Array{Array{Float64,1}})

First, this makes bit representations of T=1 (proton-proton&neutron-neutron) interactions for each {m_z}.

"""
function prep_Hamil_T1(p_sps::Array{Array{Int64,1}},n_sps::Array{Array{Int64,1}},
                 msps_p::Array{Array{Int64,1},1},msps_n::Array{Array{Int64,1},1},
                 pbits,nbits,
                 labels::Array{Array{Array{Int64,1},1},1},TBMEs::Array{Array{Float64,1}})
    lp = length(msps_p)
    ln = length(msps_n)
    bV1 = [ bit2b[] for i=1:2 ]
    V1 = [ Float64[] for i=1:2]
    mstates = [msps_p,msps_n]
    sps = [p_sps,n_sps]
    loffs = [ 0, length(p_sps)]

    for vrank =1:2 #pp:1, nn:2
        loff = loffs[vrank]
        vecs= [ [ [ false for i = 1:lp] for j=1:2],
                [ [ false for i = 1:ln] for j=1:2]]
        blist=bit2b[]
        Vs=Float64[]
        @inbounds for (i,ME) in enumerate(TBMEs[vrank])
            a,b,c,d,totJ,dummy = labels[vrank][i]            
            J2  = 2*totJ
            ja,ma_s,ma_idxs = possible_mz(sps[vrank][a-loff],mstates[vrank])
            jb,mb_s,mb_idxs = possible_mz(sps[vrank][b-loff],mstates[vrank])
            jc,mc_s,mc_idxs = possible_mz(sps[vrank][c-loff],mstates[vrank])
            jd,md_s,md_idxs = possible_mz(sps[vrank][d-loff],mstates[vrank])
            for (ic,mc) in enumerate(mc_s)
                for (id,md) in enumerate(md_s)
                    if c == d && mc >= md; continue;end
                    if abs(mc + md) > J2; continue;end
                    M_ani = mc + md
                    initialize_tvec!(vecs[vrank][1]);vecs[vrank][1][mc_idxs[ic]] = true
                    bit_c = bitarr_to_int(vecs[vrank][1])
                    initialize_tvec!(vecs[vrank][1]); vecs[vrank][1][md_idxs[id]] = true
                    bit_d = bitarr_to_int(vecs[vrank][1])
                    for (ia,ma) in enumerate(ma_s)
                        for (ib,mb) in enumerate(mb_s)
                            if a == b && ma>=mb; continue;end
                            if ma + mb != M_ani;continue;end
                            initialize_tvec!(vecs[vrank][2]);vecs[vrank][2][ma_idxs[ia]] = true
                            bit_a = bitarr_to_int(vecs[vrank][2])
                            initialize_tvec!(vecs[vrank][2]);vecs[vrank][2][mb_idxs[ib]] = true
                            bit_b = bitarr_to_int(vecs[vrank][2])
                            CG1 = clebschgordan(Float64,ja//2,ma//2,jb//2,mb//2, J2//2, M_ani//2)
                            CG2 = clebschgordan(Float64,jc//2,mc//2,jd//2,md//2, J2//2, M_ani//2)
                            tl = bit2b(bit_a,bit_b,bit_c,bit_d)
                            if ( tl in blist) == false
                                push!(blist,tl)
                                push!(Vs, ME * sqrt( (1.0+deltaf(a,b)) *(1.0+deltaf(c,d)) ) * CG1 * CG2)
                                # 2024/07/08
                                if a == c && b == d == 2 && ma == mc == -1  && mb == md == 3
                                    println("Def: ma,mb = $ma, $mb  mc,md = $mc, $md J $totJ CG = $CG1, $CG2 ME = $ME")                                    
                                end
                                #
                                continue
                            end
                            for kk = 1:length(blist)
                                if blist[kk] == tl
                                    Vs[kk] += ME * sqrt( (1.0+deltaf(a,b)) *(1.0+deltaf(c,d)) ) * CG1 * CG2
                                    # 2024/07/08
                                    if  a == c && b == d == 2 && ma == -1  && mb == md == 3
                                        println("Add: ma,mb = $ma, $mb  mc,md = $mc, $md J $totJ CG = $CG1, $CG2 ME = $ME")                                    
                                    end
                                    #
                                    break
                                end
                            end
                        end
                    end
                end
            end
        end
        bV1[vrank] = blist
        V1[vrank] = Vs
    end
    pp_2bjump = prep_2bjumps(msps_p,pbits,bV1[1],V1[1])
    nn_2bjump = prep_2bjumps(msps_n,nbits,bV1[2],V1[2])

    
    return pp_2bjump,nn_2bjump
end

"""
    prep_Hamil_pn(p_sps::Array{Array{Int64,1}},n_sps::Array{Array{Int64,1}},
           msps_p::Array{Array{Int64,1},1},msps_n::Array{Array{Int64,1},1},
           labels::Array{Array{Int64,1}},TBMEs::Array{Float64,1},zeroME=false)

make bit representation of T=0 (proton-neutron) interactions for each {m_z}
"""
function prep_Hamil_pn(p_sps::Array{Array{Int64,1}},
                n_sps::Array{Array{Int64,1}},
                msps_p::Array{Array{Int64,1},1},
                msps_n::Array{Array{Int64,1},1},
                labels::Array{Array{Int64,1}},
                TBMEs::Array{Float64,1},
                zeroME=false)
    lp = length(msps_p); ln = length(msps_n)
    loff = length(p_sps)
    Mzs = Int64[]
    for j = 1:lp
        for i = 1:lp
            push!(Mzs,msps_p[i][5]-msps_p[j][5])
        end
    end
    unique!(Mzs);sort!(Mzs,rev=true)
    lenMz = length(Mzs)

    bVpn= [ [[1,1]] for i=1:lenMz ]
    for i = 1:lenMz;deleteat!(bVpn[i],1);end
    Vpn=[ Float64[ ] for i=1:lenMz]

    ret = [-1,-1,-1]
    vec_ani_p = [false for i = 1:lp];vec_cre_p = [false for i = 1:lp]
    vec_ani_n = [false for i = 1:ln];vec_cre_n = [false for i = 1:ln]

    @inbounds for (i,ME) in enumerate(TBMEs)
        a,b,c,d,totJ,dummy = labels[i]
        J2  = 2*totJ
        ja,ma_s,ma_idxs = possible_mz(p_sps[a],msps_p)
        jc,mc_s,mc_idxs = possible_mz(p_sps[c],msps_p)
        jb,mb_s,mb_idxs = possible_mz(n_sps[b-loff],msps_n)
        jd,md_s,md_idxs = possible_mz(n_sps[d-loff],msps_n)
        ja = ja//2; jb=jb//2;jc = jc//2; jd=jd//2
        for (ic,mc) in enumerate(mc_s)
            initialize_tvec!(vec_ani_p); vec_ani_p[mc_idxs[ic]] = true
            bit_c = bitarr_to_int(vec_ani_p)
            for (ia,ma) in enumerate(ma_s)
                initialize_tvec!(vec_cre_p); vec_cre_p[ma_idxs[ia]] = true
                bit_a = bitarr_to_int(vec_cre_p)
                Mp = ma - mc
                bisearch!(Mzs,Mp,ret);idx = ret[1]
                tV = Vpn[idx];  bV = bVpn[idx]
                for (id,md) in enumerate(md_s)
                    if abs(mc + md) > J2; continue;end
                    initialize_tvec!(vec_ani_n); vec_ani_n[md_idxs[id]] = true
                    bit_d = bitarr_to_int(vec_ani_n)
                    CG1 = clebschgordan(Float64,jc,mc//2,jd,md//2,J2//2,(mc+md)//2)
                    tfac = ME .* CG1
                    for (ib,mb) in enumerate(mb_s)
                        if mb - md + Mp != 0; continue;end
                        initialize_tvec!(vec_cre_n); vec_cre_n[mb_idxs[ib]] = true
                        bit_b = bitarr_to_int(vec_cre_n)
                        fac = tfac .* clebschgordan(Float64,ja,ma//2,jb,mb//2,J2//2,(ma+mb)//2)
                        tl = [bit_a,bit_b,bit_c,bit_d]
                        if (tl in bV) == false
                            push!(bV,copy(tl))
                            push!(tV,fac)
                            continue
                        end
                        for kk = 1:length(bV)
                            if bV[kk] == tl
                                tV[kk] += fac
                                break
                            end
                        end
                    end
                end
            end
        end
    end
    return bVpn,Vpn,Mzs
end

function bitarr_to_int(arr::Array{Bool,1}, val = 0)
    v = 2^(length(arr)-1)
    for i in eachindex(arr)
        val += v*arr[i]
        v >>= 1
    end
    return val
end
function bitarr_to_int!(arr::Array{Bool,1},ret)
    ret[1] = 0
    ret[2] = 2^(length(arr)-1)
    for i in eachindex(arr)
        ret[1] += ret[2]*arr[i]
        ret[2] >>= 1
    end
    return nothing
end

function deltaf(i::Int64,j::Int64)
    ifelse(i==j,1.0,0.0)
end

"""
    occ(p_sps::Array{Array{Int64,1}},msps_p::Array{Array{Int64,1}},mzp::Array{Int64,1},num_vp::Int64,
        n_sps::Array{Array{Int64,1}},msps_n::Array{Array{Int64,1}},mzn::Array{Int64,1},num_vn::Int64,Mtot::Int64;pnpair_truncated=false)

prepare bit representations of proton/neutron Slater determinants => pbits/nbits

jocc_p, jocc_n: corresponding occupation numbers for a "j" state,  
which is used for one-body operation and OBTDs.  

Mps/Mns: total M for proton/neutron "blocks"

For 6Li in the p shell and M=0, Mps = [-3,-1,1,3] & Mns = [3,1,-1,-3]  
blocks => [ (Mp,Mn)=(-3,3),(Mp,Mn)=(-1,1),...]  

tdims: array of cumulative number of M-scheme dimensions for "blocks"  

tdims =[ # of possible configurations of (-3,3),  
         # of possible configurations of (-1,1),...]  
"""
function occ(p_sps::Vector{Vector{Int64}}, msps_p::Vector{Vector{Int64}}, mzp::Array{Int64,1}, num_vp::Int64,
             n_sps::Vector{Vector{Int64}}, msps_n::Vector{Vector{Int64}}, mzn::Array{Int64,1}, num_vn::Int64, Mtot::Int64;
             truncation_scheme="",truncated_jocc=Dict{String,Vector{Int64}}())
    lp = length(mzp); ln = length(mzn)
    all_pDim = binomial(lp,num_vp)
    all_occs_p =[ [false for j=1:lp] for i = 1:all_pDim]
    all_perm!(lp,num_vp,all_occs_p)

    all_nDim = binomial(ln,num_vn)
    all_occs_n =[ [false for j=1:ln] for i = 1:all_nDim]
    all_perm!(ln,num_vn,all_occs_n)

    occs_p = Vector{Bool}[ ]
    occs_n = Vector{Bool}[ ]

    mdim = 0;tdims=Int64[];push!(tdims,0)
    Mret = Int64[0]; Mps = Int64[]; Mns = Int64[]
    for i=1:all_pDim
        tocc = all_occs_p[i]
        Mcount!(lp,mzp,tocc,Mret)
        if !check_jocc_truncation(-1,msps_p,tocc,truncated_jocc,truncation_scheme); continue; end
        push!(occs_p,copy(tocc))
        push!(Mps, Mret[1])
    end
    for i=1:all_nDim
        tocc = all_occs_n[i]
        Mcount!(ln,mzn,tocc,Mret)
        if !check_jocc_truncation(1,msps_n,tocc,truncated_jocc,truncation_scheme); continue;end
        push!(occs_n,copy(tocc))
        push!(Mns, Mret[1])
    end
    pDim = length(occs_p); nDim = length(occs_n)

    Mps = unique(Mps); Mns = unique(Mns)
    sort!(Mps,rev=false); sort!(Mns,rev=true)
    possidxs = Vector{Int64}[]
    for i = 1:length(Mps)
        for j = 1:length(Mns)
            if Mps[i] + Mns[j] == Mtot
                push!(possidxs,[i,j])
            end
        end
    end
    pbits = [ Int64[] for i=1:length(possidxs) ]
    nbits = [ Int64[] for i=1:length(possidxs) ]
    occ_p_j = [ [[0.0,0.0]] for i=1:length(possidxs) ]
    occ_n_j = [ [[0.0,0.0]] for i=1:length(possidxs) ]
    for i =1:length(possidxs)
        deleteat!(occ_p_j[i],1);deleteat!(occ_n_j[i],1)
    end
    tocc_p_j = zeros(Float64,length(p_sps))
    tocc_n_j = zeros(Float64,length(n_sps))
    for ith = 1:length(possidxs)
        ni,nj = possidxs[ith]
        Mp = Mps[ni]; Mn = Mns[nj]
        if Mp + Mn != Mtot;@error "something wrong in occ func" ;end
        for i = 1:pDim
            Mcount!(lp,mzp,occs_p[i],Mret)
            if Mp != Mret[1]; continue;end
            push!(pbits[ith], bitarr_to_int(occs_p[i]))
            count_jocc!(p_sps,msps_p,occs_p[i],tocc_p_j)
            push!(occ_p_j[ith],copy(tocc_p_j))
            for j=1:nDim
                Mcount!(ln,mzn,occs_n[j],Mret)
                if Mp + Mret[1] != Mtot; continue; end   
                mdim += 1
            end
        end
        for j=1:nDim
            Mcount!(ln,mzn,occs_n[j],Mret)
            if Mn != Mret[1]; continue; end   
            nbit = bitarr_to_int(occs_n[j])
            push!(nbits[ith], nbit)
            count_jocc!(n_sps,msps_n,occs_n[j],tocc_n_j)
            push!(occ_n_j[ith],copy(tocc_n_j))
        end
        push!(tdims,mdim)
    end
    return pbits,nbits,occ_p_j,occ_n_j,Mps,Mns,tdims
end

function TF_connectable(Phi::Int64,bVs::bit2b,TF::Array{Bool,1})
    # bVs = [ bit(Ca^+), bit(Cb^+), bit(Cc), bit(Cd) ]
    if bVs.c != (bVs.c & Phi);TF[1]=false;return nothing;end
    if bVs.d != (bVs.d & Phi);TF[1]=false;return nothing;end
    APhi =  (~bVs.d) & ( (~bVs.c) & Phi)
    if bVs.a != (bVs.a & (~APhi));TF[1]=false;return nothing;end
    if bVs.b != (bVs.b & (~APhi));TF[1]=false;return nothing;end
    TF[1] = true
    return nothing
end

function TF_connectable_1(Phi::Int64,a::Int64,c::Int64,TF::Array{Bool,1})
    # a: creation c: anihilation
    if c != c & Phi;TF[1]=false;return nothing;end
    if a != (a & (~((~c) & Phi)));TF[1]=false;return nothing;end
    TF[1] = true
    return nothing
end

function calc_phase!(Phi::Int64,bVs::bit2b,ln::Int64,ret::Array{Int64,1})
    ret[1] = 0; ret[2] = -1 # ret=[phase,nPhi]
    ## anihilation c
    ret[1] += count_ones(Phi& (~(bVs.c-1)))
    ret[2] = (~bVs.c) & Phi #|phi>:= Cc|v>
    ## anihilation d
    ret[1] += count_ones(ret[2]& (~(bVs.d-1)))
    ret[2] = (~bVs.d) & ret[2] #|phi>:=CdCc|v>
    ## creation b
    ret[1] += count_ones(ret[2]& (~(bVs.b-1)))
    ret[2] = bVs.b | ret[2] #|phi>:=C^+bCdCc|v>
    ## creation a
    ret[1] += count_ones(ret[2]& (~(bVs.a-1)))
    ret[2] = bVs.a | ret[2] #|phi>:=C^+aC^+bCdCc|v>
    ret[1] = (-1)^ret[1]
    return nothing
end

function calc_phase_1!(Phi::Int64,cre::Int64,ani::Int64,ret::Array{Int64,1})
    ret[1] = 1; ret[2] = -1
    ## anihilation #|phi>:=Cc|v>
    ret[1] += count_ones(Phi & (~(ani-1)))
    ret[2] = (~ani) & Phi
    ## creation  #|phi>:=C^+aCc|v>
    ret[1] += count_ones(ret[2] & (~(cre-1)))
    ret[2] = cre | ret[2]
    ret[1] = (-1)^ret[1]
    return nothing
end

function count_jocc!(sps,mstates,occ_bit,tocc)
    tocc .= 0.0
    for i = 1:length(occ_bit)
        if occ_bit[i]
            tocc[mstates[i][6]] += 1.0
        end
    end
    return nothing
end
function getZNA(target_el,Anum,cp,cn)
    Z = 0
    for i = 1:length(element)
        if element[i] == target_el;Z = i; break; end
    end
    N = Anum - Z
    vp = Z - cp; vn = N - cn
    return Z,N,vp,vn
end

function ThickRestart_J(vks,uks,Tmat,lm,ls,eval_jj,Jtol)
    vals,vecs = eigen(@views Tmat[1:lm-1,1:lm-1])
    k = argmin(abs.(vals.-eval_jj))
    beta = Tmat[lm-1,lm]
    Tmat .= 0.0
    Tmat[1,1] = vals[k]
    uk=uks[1]
    tmp = beta * vecs[lm-1,k]
    Tmat[ls+1,k] = tmp; Tmat[k,ls+1] = tmp
    uk .= 0.0
    @inbounds for j=1:lm-1
        axpy!(vecs[j,k],vks[j],uk)
    end
    vks[1] .= uk .* (1.0/sqrt(dot(uk,uk)))
    vks[2] .= vks[lm]
    for k = ls+2:length(vks)
        vks[k] .= 0.0
    end
    return nothing
end

function initialize_wf(v1,method,tJ,mdim)
    Random.seed!(1234)
    if method=="rand" || tJ == -1
        v1 .= randn(mdim,)
        v1 ./= sqrt(dot(v1,v1))
    elseif method == "one"
        v1 .= 1.0 / sqrt(mdim)
    end
    return nothing
end

function initialize_bl_wav(mdim,q,vks)
    Random.seed!(1234)
    for i=1:q
        v = @views vks[i,:]
        v .= randn(mdim,)
        if i>1
            for k=i-1:-1:1
                tv = @views vks[k,:]
                v .-= dot(v,tv) .* tv
            end
        end
        tmp = 1.0/sqrt(dot(v,v))
        v.*= tmp
    end
    return nothing
end

function read_wav(inpf,mdim,n;all=false,verbose=false)
    fid =open(inpf)
    neig = read(fid,Int32)
    mtot = read(fid,Int32)
    Es = [read(fid,Float64) for i = 1:neig]
    jj = [read(fid,Int32) for i=1:neig]
    V = [ [read(fid,Float64) for i=1:mdim] for nth=1:neig]
    close(fid)
    if all
        return V,jj
    else
        return V[n],jj[n]
    end
end

function read_appwav(inpf,mdim,V,q,verbose=false,is_block=false)
    fid =open(inpf)
    neig = read(fid,Int32)
    if neig < q; println("#(appwav) must be >= q");exit();end
    mtot = read(fid,Int32)
    Es = [read(fid,Float64) for i = 1:neig]
    jj = [read(fid,Int32) for i=1:neig]
    if verbose
        println("appwav: $inpf")
        print_vec("Es(n=$neig)",Es)
    end
    if is_block
        for j=1:q
            V[j,:] .= [read(fid,Float64) for i=1:mdim]
        end
    else        
        V .= [read(fid,Float64) for i=1:mdim]
    end
        
    close(fid)
    return nothing
end

function ReORTH(it,vtarget,vks)
    for l = it:-1:1
        v = vks[l]
        alpha = dot(v,vtarget)
        axpy!(-alpha,v,vtarget)
    end
    return nothing
end

function Check_Orthogonality(it::Int,vks,en)
    svks = @views vks[1:it+1]
    for i = 1:it+1
        for j = i:it+1
            tdot = dot(svks[i],svks[j])
            if (abs(tdot) > 1.e-10 && i != j) || (i==j && abs(1-tdot) > 1.e-10)
                println("dot(",@sprintf("%3i",i),",",@sprintf("%3i",j),") = ",@sprintf("%15.4e",tdot))
            end
        end
    end
    return nothing
end

function myQR!(Q,R::FA2,d1::Int64,d2::Int64) where{FA2<:Array{Float64,2}}
    R .= 0.0
    @inbounds for j = 1:d2
        for i = 1:d1
            R[j,j] += Q[i,j] .* Q[i,j]
        end
        R[j,j] = sqrt(R[j,j])
        for i = 1:d1
            Q[i,j] /= R[j,j]
        end
        for k = j+1:d2
            for i = 1:d1
                R[j,k] += Q[i,j] .* Q[i,k]
            end
            for i = 1:d1
                Q[i,k] -= Q[i,j] .* R[j,k]
            end
        end
    end
    return nothing
end

function bl_QR!(Q,R::FA2,d1::Int64,d2::Int64) where{FA2<:Array{Float64,2}}
    R .= 0.0
    @inbounds for j = 1:d2
        q = @views Q[:,j]
        t = sqrt(dot(q,q))
        R[j,j] = t
        q .*= 1.0/t
        for k = j+1:d2
            for i = 1:d1
                R[j,k] += q[i] .* Q[i,k]
            end
            for i = 1:d1
                Q[i,k] -= q[i] .* R[j,k]
            end
        end
    end
    return nothing
end

function bisearch!(v::Array{Int64,1}, target::Int64,ret::Array{Int64,1})
    hi = length(v); lo = 1
    ret[2] = 1; ret[3]=length(v)
    @inbounds while ret[2] <= ret[3]
        ret[1] = div(ret[2]+ret[3],2)
        if v[ret[1]] > target
            ret[2] = ret[1]+1
        elseif v[ret[1]] < target
            ret[3] = ret[1]-1
        else
            return nothing
        end
    end
    ret[1] = 0
    return nothing
end

function bisearch_ord!(v::Array{Int64,1}, target::Int64, ret::Array{Int64,1})
    ret[3] = length(v); ret[2] = 1
    @inbounds while ret[2] <= ret[3]
        ret[1] = div(ret[2]+ret[3],2)
        if v[ret[1]] > target
            ret[3] = ret[1]-1
        elseif v[ret[1]] < target
            ret[2] = ret[1]+1
        else
            return nothing
        end
    end
    ret[1] = 0
    return nothing
end

function func_j(j1::I,j2::I,m1::I,m2::I) where{I<:Int64}
    sqrt( (0.5*j1*(0.5*j1+1.0)-0.5*m1*(0.5*m1-1.0))*(0.5*j2*(0.5*j2+1.0)-0.5*m2*(0.5*m2+1.0)) )
end

function J_from_JJ1(JJ,tol=1.e-6)
    for J = 0:100
        hJ = 0.5*J
        if abs(hJ*(hJ+1.0)-JJ) <tol
            return hJ
        end
    end
    return -1.0
end

function prep_2bjumps(msps::Array{Array{Int64,1},1},bits::Array{Array{Int64,1}},bVT1::Array{bit2b,1},VT1::Array{Float64,1}) 
    n_block=length(bits)
    lmsps = length(msps)
    twobody_jump = [ [ T1info[] for j=1:length(bits[i]) ] for i = 1:n_block]

    @inbounds @threads for block_i = 1:n_block # number of proton subblock
        TF=[true]; ret=[1,-1]; ridx=[-1,-1,-1]
        for (i_idx_conf,Phi) in enumerate(bits[block_i])
            for nth = 1:length(VT1)
                if VT1[nth] == 0.0; continue;end
                TF_connectable(Phi,bVT1[nth],TF)
                if TF[1]==false; continue;end
                calc_phase!(Phi,bVT1[nth],lmsps,ret)#ret=>phase,final p/n bit
                bisearch!(bits[block_i],ret[2],ridx)#ridx=Npf
                if i_idx_conf <= ridx[1]
                    coef = VT1[nth]*ret[1]*ifelse(i_idx_conf==ridx[1],0.5,1.0)
                    push!(twobody_jump[block_i][i_idx_conf], T1info(ridx[1],coef))
                end
            end
        end
    end
    return twobody_jump
end

function prep_pn(lblock::Int64,pbits::Array{Array{Int64,1}},nbits::Array{Array{Int64,1}},
                 Mps::IA,delMs::IA,bVpn::Array{Array{Array{Int64,1},1},1},Vpn::Array{Array{Float64,1},1}) where {IA<:Array{Int64,1}}
    #to = TimerOutput()
    maxDeltaM = maximum(delMs)
    bfs = [ Int64[ ] for i = 1:lblock]
    bis = [ Int64[ ] for i = 1:lblock]
    for i = 1:lblock
        for j = i:lblock
            if Mps[j] - Mps[i] in delMs
                push!(bfs[i],j); push!(bis[j],i)
            end
        end
    end

    p_NiNfs = [ [ [ ifph[] ] for j=1:length(bfs[i]) ] for i=1:lblock]
    n_NiNfs = [ [ [ ifph[] ] for j=1:length(bfs[i]) ] for i=1:lblock]
    hits = [ [0,0] for i=1:lblock]
    @threads for bi = 1:lblock
        ret = [0,0,0]
        for j = 1:length(bfs[bi])
            bf = bfs[bi][j]
            deltaM = Mps[bf] - Mps[bi]
            if (deltaM in delMs) == false;continue;end
            bisearch!(delMs,deltaM,ret) # Vidx=ret[1]
            hits[bi] += calc_1b_jumps!(bi,bf,j,pbits,nbits,bVpn[ret[1]],Vpn[ret[1]],p_NiNfs,n_NiNfs)
        end
    end
    return bis,bfs,p_NiNfs,n_NiNfs, hits
end

function calc_1b_jumps!(bi::Int64,bf::Int64,j::Int64,pbits::AAI,nbits::AAI,bVpn::AAI,Vpn::Array{Float64,1},
                        p_NiNfs::Aifph,n_NiNfs::Aifph) where {AAI<:Array{Array{Int64,1},1},Aifph<:Array{Array{Array{Array{ifph,1},1},1},1}}
    TF=[true]; ret=[1,-1] # [phase,possible]
    ridx=[-1,-1,-1]
    plist = [ ifph[ ] for i =1:length(Vpn) ]
    nlist = [ ifph[ ] for i =1:length(Vpn) ]
    hits = [0,0]
    @inbounds for nth = 1:length(Vpn)
        a,b,c,d = bVpn[nth]; if Vpn[nth] == 0.0; continue;end
        for (Npi,pPhi) in enumerate(pbits[bi])# Npi:idx pPhi:bit pSD
            TF_connectable_1(pPhi,a,c,TF)
            if TF[1]==false; continue;end
            calc_phase_1!(pPhi,a,c,ret)# ret<=[phase,bit]
            bisearch!(pbits[bf],ret[2],ridx)# ridx=[Npf]
            if ridx[1] ==0 ;continue;end
            push!(plist[nth],ifph(Npi,ridx[1],ret[1]==-1))
            hits[1] += 1
        end
        for (Nni,nPhi) in enumerate(nbits[bi]) # Nni:idx
            TF_connectable_1(nPhi,b,d,TF);if TF[1]==false; continue;end
            calc_phase_1!(nPhi,b,d,ret)
            bisearch!(nbits[bf],ret[2],ridx) # ridx=[Nnf]
            if ridx[1] == 0;continue;end
            push!(nlist[nth],ifph(Nni,ridx[1],ret[1]==-1))
            hits[2] += 1
        end
    end
    p_NiNfs[bi][j] = plist
    n_NiNfs[bi][j] = nlist
    return hits
end

function add_bl_T!(q::Int64,k::Int64,Tmat::FA2,R::FA2) where{FA2<:Array{Float64,2}}
    Tmat[q*k+1:q*k+q,q*k-q+1:q*k] .= R
    Tmat[q*k-q+1:q*k,q*k+1:q*k+q] .= R'
    return nothing
end

function vec_to_block(vecs,V)
    for (b,v) in enumerate(vecs)        
        V[:,b] .= v
    end
    return nothing
end


function writeRitzvecs(mdim,mtot,vals,totJs,Rvecs,oupf)
    num_ev = length(vals)
    io = open(oupf,"w")
    write(io,Int32(num_ev))
    write(io,Int32(mtot))
    for i = 1:num_ev
        write(io,vals[i])
    end
    js= [ Int32(2*totJs[i]) for i=1:num_ev]
    for i = 1:num_ev
        write(io,Int32(js[i]))
    end
    for nth = 1:num_ev; write(io,Rvecs[nth]); end
    close(io)
end

function make_distribute(num_task)
    lblock = length(num_task)
    n = nthreads()
    r = div(lblock,n)
    tasks = [ num_task[bi][1]*num_task[bi][2] for bi=1:lblock]
    idxs = sortperm(tasks,rev=true)
    tmp = [ Int64[ ] for i=1:n]
    imin = 1;order = true

    while length([(tmp...)...]) < lblock
        imax = minimum([imin+n-1,lblock])
        if order == true
            for i = imin:imax
                push!(tmp[i-imin+1],idxs[i])
            end
            order = false
        else
            #for i = imax:-1:imin
            for i = n:-1:1
                j = imax -(n-i)
                if j < imin
                    push!(tmp[i],0)
                else
                    push!(tmp[i],idxs[j])
                end
            end
            order = true
        end
        imin += n
    end
    task_idxs =[(tmp...)...]
    return task_idxs
end

function make_distribute_J(Jtasks)
    lblock = length(Jtasks)
    n = nthreads()
    r = div(lblock,n)
    idxs = sortperm(Jtasks,rev=true)
    tmp = [ Int64[ ] for i=1:n]
    imin = 1;order = true

    while length([(tmp...)...]) < lblock
        imax = minimum([imin+n-1,lblock])
        if order == true
            for i = imin:imax
                push!(tmp[i-imin+1],idxs[i])
            end
            order = false
        else
            for i = n:-1:1
                j = imax -(n-i)
                if j < imin
                    push!(tmp[i],0)
                else
                    push!(tmp[i],idxs[j])
                end
            end
            order = true
        end
        imin += n
    end
    Jidxs = [(tmp...)...]
    return Jidxs
end

function idx_from_ij(i::I,j::I,Dim::I) where{I<:Int64}
    return Int64((Dim+1)*(i-1) -(i*(i-1))/2+j-i+1)
end

function sparse_mat_plot(Mat::Array{Float64,2},ba="")
    n1,n2 = size(Mat)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for j = 1:n2
        for i=1:n1
            if Mat[i,j] != 0.0
                ax.scatter(i,j,marker=",",color="green", alpha=0.6)
            end
        end
    end
    plt.savefig("./mat_plot_"*ba*".pdf")
    plt.close()
end

function mdim_print(target_nuc,Z,N,cp,cn,vp,vn,mdim,tJ=-1)
    s0 = @sprintf "%6s %14s %10s %10s" target_nuc "Z,N=($Z,$N)" "c($cp,$cn)" "v($vp,$vn)"
    if tJ != -1; s0 * " 2J = " * @sprintf("%3i", tJ);end
    s1 = @sprintf "%s %12i %s %5.2f %s " "  mdim:" mdim "(" log10(mdim) ")"
    print(s0,s1)
end


function JT1(bi,oPP,oNN, p_sps,n_sps, msps_p,msps_n, pbits,nbits, tdims)
    pMat = zeros(Float64,length(pbits[bi]),length(pbits[bi]))
    nMat = zeros(Float64,length(nbits[bi]),length(nbits[bi]))
    lp = length(msps_p);ln = length(msps_n)
    TF=[true];ridx=[-1,-1,-1]
    ret_a = [-1,-1];ret_b = [-1,-1]
    ret_c = [-1,-1];ret_d = [-1,-1]
    idim = tdims[bi]
    l_Nn = length(nbits[bi])
    #### pp
    vec_ani = [false for i = 1:lp]; vec_cre=[false for i = 1:lp]
    for i2 = 1:length(p_sps)
        j2 = p_sps[i2][3]
        for m2 = -j2:2:j2-2
            vec_ani .= false; vec_cre .= false
            for k = 1:lp-1
                if @views msps_p[k][1:4] != p_sps[i2];continue;end
                if msps_p[k][5] == m2
                    vec_ani[k] = true; vec_cre[k+1] = true
                    break
                end
            end
            bitarr_to_int!(vec_ani,ret_d); bitarr_to_int!(vec_cre,ret_c)
            for (SDidx,Phi) in enumerate(pbits[bi])
                TF_connectable_1(Phi,ret_c[1],ret_d[1],TF)
                if TF[1]==false;continue;end
                APhi = ret_c[1] | ((~ret_d[1]) & Phi)
                for i1 = 1:length(p_sps)
                    j1 = p_sps[i1][3]
                    @inbounds for m1 = -j1+2:2:j1
                        vec_ani .= false; vec_cre .= false
                        for k = 2:lp
                            if @views msps_p[k][1:4] != p_sps[i1];continue;end
                            if msps_p[k][5] == m1
                                vec_ani[k] = true; vec_cre[k-1] = true
                                break
                            end
                        end
                        bitarr_to_int!(vec_ani,ret_b)
                        bitarr_to_int!(vec_cre,ret_a)
                        TF_connectable_1(APhi,ret_a[1],ret_b[1],TF)
                        if TF[1]==false;continue;end
                        NPhi = ret_a[1] | ((~ret_b[1]) & APhi)
                        bisearch!(pbits[bi],NPhi,ridx) #ridx <=Npf
                        Npi = SDidx
                        Npf = ridx[1]
                        if Npi > Npf; continue;end
                        fac = ifelse(Npi==Npf,0.5,1.0) .* func_j(j1,j2,m1,m2)
                        pMat[Npi,Npf] += fac
                    end
                end
            end
        end
    end
    #### nn
    vec_ani = [false for i = 1:ln];vec_cre = [false for i = 1:ln]
    vec_ani_2 = [false for i = 1:ln];vec_cre_2 = [false for i = 1:ln]
    for i2 = 1:length(n_sps)
        j2 = n_sps[i2][3]
        for m2 = -j2:2:j2-2
            vec_ani_2 .= false ; vec_cre_2 .= false
            for k = 1:ln-1
                if @views msps_n[k][1:4] != n_sps[i2];continue;end
                if msps_n[k][5] == m2
                    vec_ani_2[k] = true; vec_cre_2[k+1] = true;break
                end
            end
            bitarr_to_int!(vec_ani_2,ret_d); bitarr_to_int!(vec_cre_2,ret_c)
            for (SDidx,Phi) in enumerate(nbits[bi])
                TF_connectable_1(Phi,ret_c[1],ret_d[1],TF)
                if TF[1]==false;continue;end
                APhi = ret_c[1] | ((~ret_d[1]) & Phi)
                for i1 = 1:length(n_sps)
                    j1 = n_sps[i1][3]
                    for m1 = -j1+2:2:j1
                        vec_ani .= false; vec_cre .= false
                        for k = 2:ln
                            if @views msps_n[k][1:4] != n_sps[i1];continue;end
                            if msps_n[k][5] == m1
                                vec_ani[k] = true; vec_cre[k-1]=true;break
                            end
                        end
                        bitarr_to_int!(vec_ani,ret_b)
                        bitarr_to_int!(vec_cre,ret_a)
                        TF_connectable_1(APhi,ret_a[1],ret_b[1],TF)
                        if TF[1]==false;continue;end
                        NPhi = ret_a[1] | ((~ret_b[1]) & APhi)
                        bisearch!(nbits[bi],NPhi,ridx) #ridx <=Npf
                        Nni = SDidx
                        Nnf = ridx[1]
                        if Nni > Nnf; continue;end
                        fac = ifelse(Nni==Nnf,0.5,1.0) .* func_j(j1,j2,m1,m2)
                        nMat[Nni,Nnf] += fac
                    end
                end
            end
        end
    end
    op_p=MiMf[]; op_n=MiMf[]
    m = size(pMat)[1]
    for j=1:m
        for i=1:j
            if pMat[i,j] != 0.0; push!(op_p,MiMf(i,j,pMat[i,j]));end
        end
    end
    m = size(nMat)[1]
    for j=1:m
        for i=1:j
            if nMat[i,j] != 0.0; push!(op_n,MiMf(i,j,nMat[i,j]));end
        end
    end
    oPP[bi] = op_p; oNN[bi] = op_n
    return nothing
end

function prep_J(tdims,p_sps,n_sps,
                msps_p::Array{Array{Int64,1},1},
                msps_n::Array{Array{Int64,1},1},
                pbits::Array{Array{Int64,1},1},
                nbits::Array{Array{Int64,1},1})
    #to = TimerOutput()
    lp=length(msps_p); ln=length(msps_n)
    lblock = length(nbits)
    oPP = [ MiMf[] for i =1:lblock]
    oNN = [ MiMf[] for i =1:lblock]
    #@timeit to "nn//pp" begin
    @inbounds @threads for bi=1:lblock
        JT1(bi,oPP,oNN,p_sps,n_sps,msps_p,msps_n,pbits,nbits,tdims)
    end
    oPNu = [ Jpninfo[]  for i =1:lblock]
    oPNd = [ Jpninfo[]  for i =1:lblock]
    #end
    #@timeit to "pn" begin
    vec_p_ani = [false for i = 1:lp]; vec_p_cre = [false for i = 1:lp]
    vec_n_ani = [false for i = 1:ln]; vec_n_cre = [false for i = 1:ln]
    ret_a = [-1,-1];ret_b = [-1,-1]
    ret_c = [-1,-1];ret_d = [-1,-1]
    birange = [1:lblock-1,2:lblock]
    for pidx = 1:length(p_sps)
        jp = p_sps[pidx][3]
        mp_range = [-jp:2:jp-2,-jp+2:2:jp]
        for nidx = 1:length(n_sps)
            jn = n_sps[nidx][3]
            mn_range = [-jn+2:2:jn,-jn:2:jn-2]
            for ud = 1:2 # up:1 down:2
                for mp in mp_range[ud]
                    for mn in mn_range[ud]
                        fac = 1.0
                        if ud == 1
                            fac = func_j(jn,jp,mn,mp)
                        else
                            fac = func_j(jp,jn,mp,mn)
                        end
                        if fac==0.0;continue;end
                        vec_p_ani .= false; vec_p_cre .= false
                        vec_n_ani .= false; vec_n_cre .= false

                        if ud == 2  ### for pn down
                            for k = 2:lp
                                if msps_p[k][1] != p_sps[pidx][1];continue;end
                                if msps_p[k][2] != p_sps[pidx][2];continue;end
                                if msps_p[k][3] != p_sps[pidx][3];continue;end
                                if msps_p[k][4] != p_sps[pidx][4];continue;end
                                if msps_p[k][5] == mp
                                    vec_p_ani[k] = true; vec_p_cre[k-1] = true
                                    break
                                end
                            end
                            for k = 1:ln-1
                                if msps_n[k][1] != n_sps[nidx][1];continue;end
                                if msps_n[k][2] != n_sps[nidx][2];continue;end
                                if msps_n[k][3] != n_sps[nidx][3];continue;end
                                if msps_n[k][4] != n_sps[nidx][4];continue;end
                                if msps_n[k][5] == mn
                                    vec_n_ani[k] = true; vec_n_cre[k+1] = true
                                    break
                                end
                            end
                        else  ### for p up n down
                            for k = 1:lp-1
                                if msps_p[k][1] != p_sps[pidx][1];continue;end
                                if msps_p[k][2] != p_sps[pidx][2];continue;end
                                if msps_p[k][3] != p_sps[pidx][3];continue;end
                                if msps_p[k][4] != p_sps[pidx][4];continue;end
                                if msps_p[k][5] == mp
                                    vec_p_ani[k] = true;vec_p_cre[k+1] = true;break
                                end
                            end
                            for k = 2:ln
                                if msps_n[k][1] != n_sps[nidx][1];continue;end
                                if msps_n[k][2] != n_sps[nidx][2];continue;end
                                if msps_n[k][3] != n_sps[nidx][3];continue;end
                                if msps_n[k][4] != n_sps[nidx][4];continue;end
                                if msps_n[k][5] == mn
                                    vec_n_ani[k] = true;vec_n_cre[k-1] = true;break
                                end
                            end
                        end
                        bitarr_to_int!(vec_p_ani,ret_c)
                        bitarr_to_int!(vec_p_cre,ret_a)
                        bitarr_to_int!(vec_n_ani,ret_d)
                        bitarr_to_int!(vec_n_cre,ret_b)

                        @threads for bi in birange[ud]
                            TF=[true]; ret_p=[0,0];ret_n=[0,0]
                            ridx_p=[-1,-1,-1];ridx_n=[-1,-1,-1]
                            bf = bi + Int64(2*(1.5-ud))
                            pjump = ifph[]; njump=ifph[]
                            for (Npi,pPhi) in enumerate(pbits[bi])
                                TF_connectable_1(pPhi,ret_a[1],ret_c[1],TF)
                                if TF[1]==false;continue;end
                                calc_phase_1!(pPhi,ret_a[1],ret_c[1],ret_p)
                                bisearch!(pbits[bf],ret_p[2],ridx_p)
                                phase_p = ret_p[1]
                                Npf= ridx_p[1]
                                push!(pjump,ifph(Npi,Npf,phase_p==-1))
                            end
                            for (Nni,nPhi) in enumerate(nbits[bi])
                                TF_connectable_1(nPhi,ret_b[1],ret_d[1],TF)
                                if TF[1]==false;continue;end
                                calc_phase_1!(nPhi,ret_b[1],ret_d[1],ret_n)
                                bisearch!(nbits[bf],ret_n[2],ridx_n)
                                phase_n = ret_n[1]
                                Nnf= ridx_n[1]
                                push!(njump,ifph(Nni,Nnf,phase_n==-1))
                            end
                            if length(pjump)*length(njump) == 0; continue;end
                            if ud == 1
                                push!(oPNu[bi],Jpninfo(fac,pjump,njump))
                            else
                                push!(oPNd[bi],Jpninfo(fac,pjump,njump))
                            end
                        end
                    end
                end
            end
        end
    end
    #end
    return oPP,oNN,oPNu,oPNd
end

function operate_J!(Rvec,Jv,pbits,nbits,tdims,Jidxs,oPP,oNN,oPNu,oPNd,beta_J=1.0)
    #to = TimerOutput()
    lblock=length(pbits)
    #@timeit to "pp/nn " begin
    @inbounds @threads for bi in Jidxs
        if bi==0;continue;end
        idim = tdims[bi]
        lNn = length(nbits[bi])
        lNp = length(pbits[bi])
        opPP = oPP[bi]
        opNN = oNN[bi]
        offset = idim-lNn
        for tmp in opPP
            Npi =tmp.Mi; Npf=tmp.Mf; fac=tmp.fac .* beta_J
            tMi = offset + Npi*lNn
            tMf = offset + Npf*lNn
            for nidx = 1:lNn
                Mi = tMi+nidx; Mf = tMf+nidx
                Jv[Mf] += fac .* Rvec[Mi]
                Jv[Mi] += fac .* Rvec[Mf]
            end
        end
        for tmp in opNN #nn
            Nni =tmp.Mi; Nnf=tmp.Mf; fac=tmp.fac .* beta_J
            tMi = offset + Nni
            tMf = offset + Nnf
            for pidx = 1:lNp
                Mi = tMi+pidx*lNn; Mf = tMf+pidx*lNn
                Jv[Mf] += fac .* Rvec[Mi]
                Jv[Mi] += fac .* Rvec[Mf]
            end
        end
    end
    #end
    #@timeit to "pn down" begin
    @inbounds @threads for bi = 2:lblock
        operator = oPNd[bi]
        bf = bi-1
        l_Nn_i = length(nbits[bi])
        l_Nn_f = length(nbits[bf])
        off_i = tdims[bi] - l_Nn_i
        off_f = tdims[bf] - l_Nn_f
        for top in operator
            pj = top.pjump
            nj = top.njump
            fac =top.fac .* beta_J
            for tmp_p in pj
                phase_p=tmp_p.phase
                tMi = off_i + tmp_p.i * l_Nn_i
                tMf = off_f + tmp_p.f * l_Nn_f
                for tmp_n in nj
                    phase_n=tmp_n.phase
                    Mi = tMi + tmp_n.i
                    Mf = tMf + tmp_n.f
                    Jv[Mf] += ifelse(phase_p!=phase_n,-fac,fac) .* Rvec[Mi]
                end
            end
        end
    end
    @inbounds @threads for bi = 1:lblock-1
        operator = oPNu[bi]
        bf = bi+1
        l_Nn_i = length(nbits[bi])
        l_Nn_f = length(nbits[bf])
        off_i = tdims[bi] - l_Nn_i
        off_f = tdims[bf] - l_Nn_f
        for top in operator
            pj  = top.pjump
            nj  = top.njump
            fac = top.fac .* beta_J
            for tmp_p in pj
                phase_p=tmp_p.phase
                tMi = off_i + tmp_p.i * l_Nn_i
                tMf = off_f + tmp_p.f * l_Nn_f
                for tmp_n in nj
                    Nni = tmp_n.i; Nnf = tmp_n.f; phase_n=tmp_n.phase
                    Mi = tMi + tmp_n.i
                    Mf = tMf + tmp_n.f
                    Jv[Mf] += ifelse(phase_p!=phase_n,-fac,fac) .* Rvec[Mi]
                end
            end
        end
    end
    return nothing
end


mutable struct occ_object
    pbit::Vector{Bool}
    nbit::Vector{Bool}
    p_nljs::Vector{Vector{Int64}}
    n_nljs::Vector{Vector{Int64}}
    width::Int64
end

function make_occ_obj(msps_p,msps_n)
    pbit = [false for i=1:length(msps_p)]
    nbit = [false for i=1:length(msps_n)]
    p_nljs=Vector{Int64}[]
    n_nljs=Vector{Int64}[]
    width = 0
    for (pn,mstates) in enumerate([msps_p,msps_n])
        prev = jmax = 0 
        nljs = Vector{Int64}[ ]    
        for tmp in mstates
            n,l,j = tmp[1:3]
            if prev != tmp[end]
                prev = tmp[end]
                push!(nljs,[n,l,j])
            end
            jmax = max(jmax,j)
        end        
        width = max(jmax + 1 + 2,width)
        if pn == 1
            p_nljs = nljs
        else
            n_nljs = nljs
        end
    end
    return occ_object(pbit,nbit,p_nljs,n_nljs,width)
end

function get_bitarry_from_int(bit_int::Int64,numbit::Int64) 
    return reverse(digits(bit_int, base=2, pad=numbit))
end

function update_occ_obj!(occ_obj,pbit_int,nbit_int)
    pbit_arr = reverse(digits(pbit_int, base=2, pad=length(occ_obj.pbit)))
    nbit_arr = reverse(digits(nbit_int, base=2, pad=length(occ_obj.nbit)))
    occ_obj.pbit .= pbit_arr
    occ_obj.nbit .= nbit_arr
    return nothing
end

function get_occstring(occ_obj,msps_p,msps_n)
    pbit_arr = occ_obj.pbit
    nbit_arr = occ_obj.nbit
    width = occ_obj.width
    for pn = 1:2
        pntext = ""
        jtxt = "-"
        bitarr = ifelse(pn==1,pbit_arr,nbit_arr)
        mstates= ifelse(pn==1,msps_p,msps_n)
        for (bitidx,bit) in enumerate(bitarr)
            n,l,j,tz,mz,oidx = mstates[bitidx]
            mark = ifelse(bit == 1,"x","o")
            jtxt *= mark
            if mz == j
                tmp = jtxt * "-"
                while length(tmp) < width
                    tmp = "-"*tmp*"-"
                end
                tmp *= "  "*ifelse(pn==1,"π","ν")*string(n)*chara_l[l]*string(j)*"/2"
                pntext *= tmp * ifelse(bitidx!=length(pbit_arr),"\n","")
                jtxt = "-"
            end
        end
        println(pntext)
    end
    return nothing
end

"""
visualize all configurations in a given model-space
"""
function visualize_configurations(msps_p,msps_n,pbits,nbits,mdim;mdim_max=100)
    if mdim > mdim_max
        println("You are trying to execute `visualize_configurations`` for Dim.= $(mdim)!")
        println("Are you sure?(y/n, default:n)")
        s = readline( )
        if s != "y"
            println("skiped!")
            return nothing
        end
    end
    println("")
    nth = 0
    nblock = length(pbits)
    occ_obj = make_occ_obj(msps_p,msps_n)
    for block = 1:nblock
        t_pbits = pbits[block]
        t_nbits = nbits[block]
        for pbit in t_pbits
            for nbit in t_nbits
                nth += 1
                nth_txt = nth 
                if mdim > mdim_max
                    println("config: ",nth_txt)
                else 
                    println("config: ",@sprintf("%4i",nth))
                end
                update_occ_obj!(occ_obj,pbit,nbit)
                get_occstring(occ_obj,msps_p,msps_n)
            end
        end
    end
    return nothing
end

function entropy(Rvecs,pbits,nbits,tdims,to)
    pdims = [0];for bi=1:length(pbits);push!(pdims,pdims[end]+length(pbits[bi]));end
    ndims = [0];for bi=1:length(nbits);push!(ndims,ndims[end]+length(nbits[bi]));end
    pdim = pdims[end]; ndim = ndims[end]
    if ndim >= 1.e+6
        @warn  "`entropy` is function is now very naive and slow. mdim >= 1.e+6 is not recommended and it will not be done..."
    end
    rho_p = zeros(Float64,pdim,pdim)
    for (n,Rvec) in enumerate(Rvecs)
        rho_p .= 0.0
        for bi = 1:length(pbits)
            idim = tdims[bi]
            nbit = nbits[bi]
            ln = length(nbit)
            pofst = pdims[bi]
            for (pid_i,pbit_i) in enumerate(pbits[bi])
                tMi = idim + (pid_i-1)*ln
                ri = pofst + pid_i
                for (pid_j,pbit_j) in enumerate(pbits[bi])
                    if pid_j <pid_i;continue;end
                    tMj = idim + (pid_j-1)*ln
                    rj = pofst + pid_j
                    for nidx =1:ln
                        Mi = tMi+nidx; Mj = tMj+nidx
                        tmp = Rvec[Mi] .* Rvec[Mj]
                        rho_p[ri,rj] += tmp
                        if ri != rj;rho_p[rj,ri] = rho_p[ri,rj];end
                    end
                end
            end
        end
        s = -tr( rho_p * log(rho_p) )
        println("n= $n entropy $s")
        if n == 1
            for i = 1:pdim
                print_vec("",rho_p[i,:])
            end
            #svd of rho_p
            SVD = svd(rho_p)
            U = SVD.U; Sig = Diagonal(SVD.S); Vt = SVD.Vt; V = Vt'
            println("U")
            for j = 1:pdim
                print_vec("",U[j,:];long=false)
            end 
            println("Sig")
            for j = 1:pdim
                print_vec("",Sig[j,:];long=false)
            end 
            println("V")
            for j = 1:pdim
                print_vec("",V[j,:];long=false)
            end 
            println("")
        end    
    end
    return nothing
end
