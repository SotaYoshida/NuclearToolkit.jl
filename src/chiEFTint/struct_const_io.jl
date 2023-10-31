const chara_L = ["S","P","D","F","G","H","I","J","K","L","M","N"]
const chara_l = ["s","p","d","f","g","h","i","j","k","l","m","n"]
const jmax = 6
const lmax = jmax +1
const lcmax = jmax + 1
const Mp = 938.272088 
const Mn = 939.565420 
const Mm = (Mp+Mn)/2
const Me = 0.51099895
const Mred = 2*Mp*Mn/(Mp+Mn)
const Ms = [Mp,Mm,Mn]
const Lambchi = 500.0
const hc = 197.327053
const Fpi = 92.4
const hc2 = hc^2
const hc3 = hc^3
const hc4 = hc^4
const gA = -1.29
const gA2 = gA^2
const gA4 = gA^4
const mpis = [139.5702,134.9766,139.5702]
const fsalpha = 7.2973525693* 1.e-3 #fine structure const. 
const l2l = [ wigner3j(Float64,l,2,l,0,0,0) for l=0:8]
const l2lnd =[[ wigner3j(Float64,l1,2,l2,0,0,0) for l2=0:8] for l1=0:8]
const nuclist = [
     "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F", "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]

delta(a::Int64,b::Int64) = ifelse(a==b,1.0,0.0)
delta(a::Float64,b::Float64) = ifelse(a==b,1.0,0.0)
hat(a::Int64) = sqrt(2.0*a+1.0)
hat(a::Float64) = sqrt(2.0*a+1.0)
hahat(a::Int64) = 2.0*a+1.0

function chara_SLJ(S,L,Lp,J) 
    if L==Lp
        return string(Int(2*S+1))*chara_L[L+1]*string(J)
    else
        return string(Int(2*S+1))*chara_L[L+1]*"-"*chara_L[Lp+1]*string(J)
    end
end

"""
    QLs
second kind Legendre functions needed for pion exchange contributions
"""
struct QLs
    QL0s::Vector{Matrix{Float64}}
    QWL2s::Vector{Matrix{Float64}}
    QWL1s::Vector{Matrix{Float64}}
    QWL0s::Vector{Matrix{Float64}}
    QWlss::Vector{Matrix{Float64}}
    ndQW1s::Vector{Vector{Matrix{Float64}}}
    ndQW2s::Vector{Vector{Matrix{Float64}}}
    QLdict::Vector{Dict{Float64,Float64}}
end
"""
    Fis_2n3n
Struct to store F-integrals for 2n3n contributions, Eqs.(A13)-(A16) in [Kohno2013].
"""
struct Fis_2n3n
    F0s::Vector{Float64}
    F1s::Vector{Float64}
    F2s::Vector{Float64}
    F3s::Vector{Float64}
end
"""
    wsyms_j1_1or2
Struct to store Wigner symbols with some special cases for 2n3n.
"""
struct wsyms_j1_1or2
    cg1s::Matrix{Float64}
    cg2s::Matrix{Float64}
    d6_121::Array{Float64,3}
    d6_21::Array{Float64,4}
    d6_222::Array{Float64,3}
    d9_12::Array{Float64,4}
end
"""
    util_2n3n
Struct to utilize `Fis_2n3n`,`QLs`,`wsyms_j1_1or2` structs.
"""
struct util_2n3n
    Fis::Fis_2n3n
    QWs::QLs
    wsyms::wsyms_j1_1or2
end
"""
    LECs
Struct to store low energy constants (LECs).
Since the vector-form is usuful for sampling, dictonaries `idxs`&`dLECs` are defined to get idx from string and vice versa.
"""
struct LECs
    vals::Vector{Float64}
    idxs::Dict{String,Int64}
    dLECs::Dict{String,Float64}
end

"""
    chiEFTparams
mutable struct to specify parameters for chiEFTint
# Fields
- `n_mesh::Int64` # of momentum mesh points
- `pmax_fm::Float64` max momentum in fm``{}^{-1}``
- `emax::Int64` maximum emax (2n+l) for single particle states
- `Nnmax::Int64` Nnmax quanta
- `chi_order::Int64` order of Chiral EFT (0:LO, 1:NLO, 2:NNLO, 3:NNLO, 4:N4LO)
- `calc_NN::Bool` calculate NN potential or not
- `calc_3N::Bool` calculate density-dependent 3NF (called 2n3n in this package) or not
- `coulomb::Bool` calculate coulomb term or not
- `calc_EperA::Bool` calculate E/A with HF(+MBPT) the kF below will be used for density
- `hw::Float64` oscillator parameter in MeV
- `srg::Bool` carrying out SRG or not
- `srg_lambda::Float64` resolution scale for free space SRG in fm``{}^{-1}``
- `tbme_fmt::String` tbme format (snt, snt.bin)
- `fn_tbme::String` filename of output tbme
- `pottype::String` potential type (em500n3lo, emn500n3lo, emn500n4lo)
- `LambdaSFR::Float64` cutoff for spectral function regularization (SFR)
- `target_nlj::Vector{Vector{Int64}}` option to truncate {nlj} in an output snt
- `v_chi_order::Int64` order of valence chiral interaction (usually not used, i.e. 0)
- `n_mesh_P::Int64` # of momentum mesh for valence interaction  
- `Pmax_fm::Float64` momentum mesh for valence interaction
- `kF::Float64` Fermi momentum for 2n3n 
- `BetaCM::Float64` Lawson's beta for Hcm
"""
mutable struct chiEFTparams
    n_mesh::Int64
    pmax_fm::Float64
    emax::Int64
    Nnmax::Int64
    chi_order::Int64 
    calc_NN::Bool
    calc_3N::Bool
    coulomb::Bool
    calc_EperA::Bool
    hw::Float64
    srg::Bool
    srg_lambda::Float64
    tbme_fmt::String
    fn_tbme::String
    pottype::String
    LambdaSFR::Float64   
    target_nlj::Vector{Vector{Int64}}
    v_chi_order::Int64
    n_mesh_P::Int64
    Pmax_fm::Float64
    kF::Float64
    BetaCM::Float64
end

"""
    ChiralEFTobject
# Fields
- `params::chiEFTparams`
- `xr_fm::Vector{Float64}` momentum mesh points in fm``^{-1}`` (0 to `pmax_fm`)
- `xr::Vector{Float64}` momentum mesh points in MeV (`xr_fm` times ``\\hbar c``)
- `wr::Vector{Float64}` weights vector for Gauss-Legendre quadrature
- `Rnl::Array{Float64,3}` arrays for radial functions
- `xrP_fm::Vector{Float64}` "valence" counter part of `xr_fm` (usually not used)
- `xrP::Vector{Float64}` "valence" counter part of `xr` (usually not used)
- `wrP::Vector{Float64}` "valence" counter part of `wr` (usually not used)
- `RNL::Array{Float64,3}` "valence" counter part of `Rnl` (usually not used)
- `lsjs::Vector{Vector{Vector{Int64}}}`
- `tllsj::Vector{Int64}`
- `opfs::Vector{Vector{Float64}}`
- `ts::Vector{Float64}` mesh points for angular integrals in pion-exchange terms
- `ws::Vector{Float64}` weights vectors for angular integrals in pion-exchange terms
- `infos::Vector{Vector{Int64}}` information of two-body channeles: { [``T_z``,parity,J,dim] } 
- `izs_ab::Vector{Vector{Vector{Int64}}}` two-body kets: { [iza,ia,izb,ib] } where `iz*` stuff is isospin (-1 or 1) and `i*` stuff is label of sps. For example, [-1,1,1,1] corresponds to |p0s1/2 n0s1/2>.
- `nTBME::Int64` # of two-body matrix elements (TBMEs)
- `util_2n3n::util_2n3n`
- `LECs::LECs`
- `V12mom::Vector{Matrix{Float64}}` matrix elements of NN int in momentum space for each two-body channnel
- `V12mom_2n3n::Vector{Matrix{Float64}}` matrix elements of 2n3n int in momentum space for each two-body channnel
- `pw_channels::Vector{Vector{Int64}}` list of partial wave two-body channels like {[pnrank,L,L',S,J]}
- `dict_pwch::Vector{Dict{Vector{Int64}, Int64}}`
- `arr_pwch::Vector{Vector{Vector{Vector{Vector{Int64}}}}}`
"""
struct ChiralEFTobject
    params::chiEFTparams
    xr_fm::Vector{Float64}
    xr::Vector{Float64}
    wr::Vector{Float64}
    Rnl::Array{Float64,3}
    xrP_fm::Vector{Float64}
    xrP::Vector{Float64}
    wrP::Vector{Float64}
    RNL::Array{Float64,3}
    lsjs::Vector{Vector{Vector{Int64}}}
    tllsj::Vector{Int64}
    opfs::Vector{Vector{Float64}}
    ts::Vector{Float64}
    ws::Vector{Float64}
    infos::Vector{Vector{Int64}} 
    izs_ab::Vector{Vector{Vector{Int64}}}
    nTBME::Int64
    util_2n3n::util_2n3n
    LECs::LECs
    V12mom::Vector{Matrix{Float64}}
    V12mom_2n3n::Vector{Matrix{Float64}}
    pw_channels::Vector{Vector{Int64}}
    dict_pwch::Vector{Dict{Vector{Int64}, Int64}} 
    arr_pwch::Vector{Vector{Vector{Vector{Vector{Int64}}}}}
end


"""
    init_chiEFTparams(;fn_params="optional_parameters.jl")
constructor of chiEFTparams, see `chiEFTparams` mutable struct for more details.
"""
function init_chiEFTparams(;fn_params="optional_parameters.jl",use_hw_formula = 0,Anum = -1,io=stdout)
    n_mesh = 50
    pmax_fm = 5.0
    emax = 4
    Nnmax= 20
    chi_order = 3
    calc_NN = true
    calc_3N = false 
    hw = 20.0
    if use_hw_formula != 0; hw = hw_formula(Anum,use_hw_formula); end
    srg = true
    srg_lambda = 2.0
    calc_EperA = false
    tx = "bare";if srg; tx ="srg"*string(srg_lambda);end;if calc_3N; tx="2n3n_"*tx;end
    tbme_fmt = "snt.bin"
    fn_tbme = "tbme_em500n3lo_"*tx*"hw"*string(round(Int64,hw))*"emax"*string(emax)*"."*tbme_fmt
    coulomb = true
    pottype = "em500n3lo"
    target_nlj=Vector{Int64}[]
    v_chi_order = 0 
    n_mesh_P = 10
    Pmax_fm = 3.0   
    BetaCM = 0.0 
    kF = 1.35
    LambdaSFR = 0.0
    params = chiEFTparams(n_mesh,pmax_fm,emax,Nnmax,chi_order,calc_NN,calc_3N,coulomb,calc_EperA,
                          hw,srg,srg_lambda,tbme_fmt,fn_tbme,pottype,LambdaSFR,
                          target_nlj,v_chi_order,n_mesh_P,Pmax_fm,kF,BetaCM)
    if !isfile(fn_params)
        println("Since $fn_params is not found, the default parameters will be used.")
        println("You can specify the parameters with optional argument, fn_params like make_chiEFTint(;fn_params=\"PATH_TO_YOUR_FILE\").")
    else
        read_chiEFT_parameter!(fn_params,params;io=io)
    end
    return params
end

"""
    read_chiEFT_parameter!(fn,params)
Function to overwrite params from the parameter file `fn`.
"""
function read_chiEFT_parameter!(fn,params::chiEFTparams;io=stdout)
    include(fn)
    if @isdefined(n_mesh); params.n_mesh = n_mesh ; end
    if @isdefined(pmax_fm); params.pmax_fm = pmax_fm ; end
    if @isdefined(emax); params.emax = emax; end
    if @isdefined(Nnmax); params.Nnmax = Nnmax; end
    if @isdefined(chi_order); params.chi_order = chi_order; end
    if @isdefined(calc_NN); params.calc_NN = calc_NN; end
    if @isdefined(calc_3N); params.calc_3N = calc_3N; end
    if @isdefined(hw); params.hw = hw; end
    if @isdefined(srg); params.srg = srg; end
    if @isdefined(srg_lambda); params.srg_lambda = srg_lambda; end
    if @isdefined(tbme_fmt); params.tbme_fmt = tbme_fmt; end
    if @isdefined(fn_tbme); params.fn_tbme = fn_tbme; end
    if @isdefined(pottype); params.pottype = pottype; end
    if @isdefined(coulomb); params.coulomb = coulomb; end
    if @isdefined(calc_EperA); params.calc_EperA = calc_EperA; end
    if @isdefined(target_nlj); params.target_nlj = target_nlj; end
    if @isdefined(v_chi_order); params.v_chi_order = v_chi_order; end
    if @isdefined(n_mesh_P); params.n_mesh_P = n_mesh_P; end
    if @isdefined(Pmax_fm); params.Pmax_fm = Pmax_fm; end
    if @isdefined(kF); params.kF = kF; end
    if occursin("emn",params.pottype)
       params.LambdaSFR = ifelse(pottype=="emn500n4lo",700.0,650.0)
    end
    if @isdefined(BetaCM); params.BetaCM = BetaCM; end
    if params.pottype =="emn500n4lo" && params.chi_order>4
        params.chi_order=4
        println("chi_order must be <= 4 for pottype=emn500n4lo, chi_order=4 will be used.") 
    end
    if (params.pottype =="emn500n3lo" || params.pottype =="em500n3lo") && params.chi_order > 3
        println("chi_order must be <= 3 for pottype=emn500n3lo or em500n3lo, chi_order=3 will be used")
        params.chi_order=3
    end
    tx = "bare";if params.srg; tx ="srg"*string(params.srg_lambda);end;if params.calc_3N; tx="2n3n_"*tx;end
    params.fn_tbme = "tbme_"*params.pottype*"_"*tx*"hw"*string(round(Int64,params.hw))*"emax"*string(params.emax)*"."*params.tbme_fmt
    if io != nothing
        println(io,"--- chiEFTparameters used ---")
        for fieldname in fieldnames(typeof(params))                 
            println(io,"$fieldname = ",getfield(params,fieldname))
        end
        println(io,"-----------------------------")
    end
    return nothing
end

function show_TimerOutput_results(to;io=stdout,tf=true,alloc=true,compact=false)
    if tf
        show(io,to, allocations = true,compact = false)
        println("")
    end
    return nothing
end

function delta_arr(a,b)
    hit = 0
    for i in eachindex(a)
        hit += ifelse(a[i]==b[i],1,0)
    end
    return ifelse(length(a)==hit,1.0,0.0)
end

function owtkey!(tkey,n,l,j,tz)
    tkey[1]=n; tkey[2]=l; tkey[3]=j; tkey[4]=tz
    return nothing
end

""" 
    readsnt(sntf,Anum;eachA=false,pnfac=1.0) 

to read a sntfile. This is slightly different from readsnt() in ShellModel.jl
"""
function readsnt(sntf,Anum;eachA=false,pnfac=1.0) 
    f = open(sntf,"r");tlines = readlines(f);close(f)
    lines = rm_comment(tlines)
    line = lines[1]
    lp,ln,cp,cn = map(x->parse(Int,x),rm_nan(split(line," ")))
    p_sps = Vector{Int64}[ ]
    n_sps = Vector{Int64}[ ]
    dictsps = Dict{Vector{Int64},Int64}()
    nls = []
    nlhit=0
    for i = 1:lp
        ith,n,l,j,tz = map(x->parse(Int,x),rm_nan(split(lines[1+i]," "))[1:5])
        push!(p_sps,[n,l,j,tz])
        if ([n,l,tz] in nls)==false
            nlhit +=1
            push!(nls,[n,l,tz])
        end
        dictsps[[n,l,j,tz]] = i
    end
    for i = 1:ln
        ith, n,l,j,tz = map(x->parse(Int,x),rm_nan(split(lines[1+i+ln]," "))[1:5])
        if ([n,l,tz] in nls)==false
            nlhit +=1
            push!(nls,[n,l,tz])
        end
        push!(n_sps,[n,l,j,tz])
        dictsps[[n,l,j,tz]] = i + lp
    end
    sps = vcat(p_sps,n_sps)    
    nsp,zero = map(x->parse(Int,x),rm_nan(split(lines[1+ln+lp+1]," "))[1:2])
    SPEs = [ [0.0 for i=1:lp],[0.0 for i=1:ln]]
    for i = 1:nsp
        idx=0; j=i
        if i<=lp;idx=1;else;idx=2;j-=lp;end
        SPEs[idx][j] =parse(Float64,rm_nan(split(lines[1+ln+lp+1+i]," "))[3])
    end
    ntbme = 0; massop = 0; Aref = 0; p=0
    tmp = rm_nan(split(lines[1+ln+lp+1+nsp+1]," "))
    if length(tmp) == 3
        ntbme,massop,hw = tmp
        ntbme = parse(Int,ntbme)
        massop=parse(Int,massop)
        hw = parse(Float64,hw)
    else
        ntbme,massop,Aref,p = tmp
        ntbme = parse(Int,ntbme);massop=parse(Int,massop)
        Aref=parse(Int,Aref); p=parse(Float64,p)
    end
    dictTBMEs=[ Dict{Vector{Int64},Float64}() for pnrank=1:3]
    for ith = 1:ntbme
        i,j,k,l,totJ,TBME= rm_nan(split(lines[1+ln+lp+1+nsp+1+ith], " "))
        i = parse(Int,i);j = parse(Int,j);k = parse(Int,k);l = parse(Int,l);
        totJ = parse(Int,totJ)
        nth = 0
        if i<=lp && j<=lp
            nth = 1
        elseif i>lp && j > lp
            nth = 3
        elseif i<=lp && j>lp
            nth = 2
        else
            println("i $i j $j k $k l $l totJ $totJ TBME $TBME")
            println("err");exit()
        end
        TBME = parse(Float64,TBME)
        if eachA && massop == 1
            TBME *= (Anum/Aref)^p 
        end        
        ## snt file must be "ordered"; a<=b & c=d & a<=c
        dictTBMEs[nth][[i,j,k,l,totJ]] = TBME *ifelse(nth==2,pnfac,1.0)
    end
    return sps,dictsps,dictTBMEs
end

"""
    write_tbme(chiEFTobj,io,ndim,izs,Jtot,vv,nljsnt,nljdict,tkeys,dWS;ofst=0)

write tbme in myg/snt(snt.bin) format
"""
function write_tbme(params::chiEFTparams,io,ndim,izs,Jtot,vv,vv_2n3n,nljsnt,nljdict,dWS,tkeys;ofst=0)
 
    tbme_fmt = params.tbme_fmt
    target_nlj = params.target_nlj
    
    @inbounds for i = 1:ndim
        iza,ia,izb,ib = izs[i]
        na,la,ja = nljsnt[ia]
        nb,lb,jb = nljsnt[ib]
        for j = 1:i
            izc,ic,izd,id= izs[j]
            nc,lc,jc = nljsnt[ic]
            nd,ld,jd = nljsnt[id]
            if tbme_fmt == "myg"
                owtkey!(tkeys[1],na,la,ja,iza)
                owtkey!(tkeys[2],nb,lb,jb,izb)
                owtkey!(tkeys[3],nc,lc,jc,izc)
                owtkey!(tkeys[4],nd,ld,jd,izd)
                vpp = kinetic_tb(tkeys[1],tkeys[2],tkeys[3],tkeys[4],Jtot,dWS)
                print(io,@sprintf("%4i", iza), @sprintf("%4i", ia), @sprintf("%4i", izb), @sprintf("%4i", ib))
                print(io,@sprintf("%4i", izc), @sprintf("%4i", ic), @sprintf("%4i", izd), @sprintf("%4i", id))
                print(io,@sprintf("%4i", Jtot),@sprintf("%20.10e", vv[i,j]),@sprintf("%20.10e", vv_2n3n[i,j]))
                println(io,@sprintf("%20.10e", vpp))
            elseif tbme_fmt == "snt" ||  tbme_fmt == "snt.bin"
                tv = vv[i,j]
                tv2n3n = vv_2n3n[i,j]
                a = ifelse(iza==-1,ia,ia+ofst)
                b = ifelse(izb==-1,ib,ib+ofst)
                c = ifelse(izc==-1,ic,ic+ofst)
                d = ifelse(izd==-1,id,id+ofst)          
                if length(target_nlj) != 0
                    if get(nljdict,a,0) == 0; continue;end
                    if get(nljdict,b,0) == 0; continue;end
                    if get(nljdict,c,0) == 0; continue;end
                    if get(nljdict,d,0) == 0; continue;end
                end         
                ra = a; rb = b; rc=c; rd=d
                fac = 1.0
                if a > b
                    ra = b; rb = a
                    fac *= (-1.0)^(div(nljsnt[a][3]+nljsnt[b][3],2)+Jtot+1)
                end
                if c > d 
                    ra = d; rd = c
                    fac *= (-1.0)^(div(nljsnt[c][3]+nljsnt[d][3],2)+Jtot+1)
                end
                fa = ra; fb =rb; fc=rc; fd=rd
                if ra > rc 
                    fa=rc;fb=rd;fc=ra;fd=rb
                else
                    if rb > rd
                        fa=rc;fb=rd;fc=ra;fd=rb
                    end
                end
                if length(target_nlj) != 0
                    fa = nljdict[fa]; fb = nljdict[fb]
                    fc = nljdict[fc]; fd = nljdict[fd]
                end               
                owtkey!(tkeys[1],na,la,ja,iza)
                owtkey!(tkeys[2],nb,lb,jb,izb)
                owtkey!(tkeys[3],nc,lc,jc,izc)
                owtkey!(tkeys[4],nd,ld,jd,izd)
                vpp = kinetic_tb(tkeys[1],tkeys[2],tkeys[3],tkeys[4],Jtot,dWS)
                if tbme_fmt == "snt"
                    print(io,@sprintf("%5i", fa),@sprintf("%5i", fb),@sprintf("%5i", fc),@sprintf("%5i", fd))
                    print(io,@sprintf("%6i", Jtot),@sprintf("%18.10f", tv),@sprintf("%18.10f", tv2n3n))
                    println(io,@sprintf("%20.10e", vpp))
                elseif tbme_fmt == "snt.bin"
                    write(io,Int16(fa));write(io,Int16(fb))
                    write(io,Int16(fc));write(io,Int16(fd));write(io, Int16(Jtot))
                    write(io,Float32(tv));write(io,Float32(tv2n3n)); write(io,Float32(vpp))
                elseif tbme_fmt == "snt.bin64"
                    write(io,Int16(fa));write(io,Int16(fb))
                    write(io,Int16(fc));write(io,Int16(fd));write(io, Int16(Jtot))
                    write(io,tv); write(io,tv2n3n);write(io,vpp)        
                end
            end
        end
    end
    return nothing
end

function select_io(MPIcomm::Bool,optimizer::String,nucs;use_stdout=false,fn="")
    io = stdout
    if MPIcomm
        @assert optimizer == "MCMC" "when using MPI for make_chiEFTint function, optimizer should be \"MCMC\""
        @assert nucs!=[] "nucs must not be empty if you set MPIcomm=true"
        MPI.Init()
        myrank = MPI.Comm_rank(MPI.COMM_WORLD)
        io = open("./mpilog_rank$(myrank).dat","w")
    elseif fn !=""
        io = open(fn,"w")
    else
        io = open("logfile.dat","w")
    end
    if use_stdout; io = stdout; end
    return io
end

function write_spes(chiEFTobj::chiEFTparams,io,nljsnt,lp,nTBME,nljdict;bin=false)
    hw = chiEFTobj.hw
    target_nlj = chiEFTobj.target_nlj
    ## header part
    if length(target_nlj)!=0
        ln = length(target_nlj)
        if bin 
            write(io, ln); write(io, ln); write(io, 0);write(io, 0)
        else
            println(io,@sprintf("%4i", ln),@sprintf("%4i",ln),
                    @sprintf("%3i", 0),@sprintf("%3i", 0))
        end
    else
        if bin 
            write(io, lp); write(io,lp); write(io,0); write(io,0)
        else
            println(io,@sprintf("%4i", lp),@sprintf("%4i",lp),
                    @sprintf("%3i", 0),@sprintf("%3i", 0))
        end 
    end
    for tz = -1:2:1
        for i in eachindex(nljsnt)
            n,l,j = nljsnt[i]
            if target_nlj != []
                ii = ifelse(tz==-1,i,i+lp)
                t = get(nljdict,ii,0)
                if t==0; continue;end
                if bin 
                    write(io,nljdict[ii])
                    write(io,n); write(io,l); write(io,j); write(io,tz)
                else
                    println(io,@sprintf("%4i",nljdict[ii]),
                            @sprintf("%4i", n),@sprintf("%4i", l),
                            @sprintf("%4i", j), @sprintf("%4i", tz))
                end
            else
                if bin 
                    write(io,ifelse(tz==-1,i,i+lp))
                    write(io,n);write(io,l);write(io,j);write(io,tz)
                else
                    println(io,@sprintf("%4i", ifelse(tz==-1,i,i+lp)),
                            @sprintf("%4i", n),@sprintf("%4i", l),
                            @sprintf("%4i", j), @sprintf("%4i", tz))
                end
            end
        end
    end
    tx = ""
    hit = 0 
    for tz = -1:2:1
        for (i,tmp) in enumerate(nljsnt)
            n,l,j=tmp
            SPE = (2*n+l+3/2) * 0.5
            if target_nlj != []
                ii = ifelse(tz==-1,i,i+lp)
                t = get(nljdict,ii,0)
                if t==0; continue;end
            end
            hit += 1
            for (k,tmp2) in enumerate(nljsnt)
                if k <= i;continue;end
                nn,ll,jj=tmp2
                if nn!=n && ll==l && jj==j 
                    hit +=1
                end
            end
        end
    end
    if bin; write(io,hit); write(io,10);write(io,hw); end
    for tz = -1:2:1
        for (i,tmp) in enumerate(nljsnt)
            n,l,j=tmp
            SPE = (2*n+l+3/2) * 0.5
            fi = i
            if target_nlj != []
                ii = ifelse(tz==-1,i,i+lp)
                t = get(nljdict,ii,0)
                if t==0; continue;end
                fi = nljdict[ii]
                if bin
                    write(io,fi);write(io,fi);write(io,SPE)
                else
                    tx *= @sprintf("%4i", fi)* @sprintf("%4i", fi)* @sprintf("%16.8f",SPE)*"\n"
                end
            else
                fi = ifelse(tz==-1,i,i+lp)
                if bin 
                    write(io,fi);write(io,fi);write(io,SPE)
                else  
                    tx *= @sprintf("%4i", fi)* @sprintf("%4i", fi)* @sprintf("%16.8f",SPE)*"\n"
                end
            end
            for (k,tmp2) in enumerate(nljsnt)
                if k <= i;continue;end
                nn,ll,jj=tmp2
                if nn!=n && ll==l && jj==j 
                    kin = kinetic_ob(tmp,tmp2)
                    if bin 
                        write(io,ifelse(tz==-1,i,i+lp));write(io,ifelse(tz==-1,k,k+lp));write(io,kin)
                    else
                        tx *= @sprintf("%4i", ifelse(tz==-1,i,i+lp))* @sprintf("%4i", ifelse(tz==-1,k,k+lp))* @sprintf("%16.8f",kin)*"\n"
                    end
                end
            end
        end
    end
    if bin         
        write(io,nTBME); write(io,10);write(io,hw)
    else
        println(io,@sprintf("%5i",hit)*@sprintf("%5i",10)*@sprintf("%8.3f", hw))
        print(io,tx)
        println(io,@sprintf("%5i",nTBME)*@sprintf("%5i",10)*@sprintf("%8.3f", hw))
    end       
    return nothing   
end

function svd_vmom(chiEFTobj::ChiralEFTobject,target_LSJ;pnrank=2,verbose=false,truncation_rank=20)   
    dict_pwch = chiEFTobj.dict_pwch[pnrank]
    V12mom = chiEFTobj.V12mom
    itt = 2 *(pnrank -2)
    for tmp in target_LSJ
        L,Lp,S,J = tmp
        V12idx = get(dict_pwch,[itt,L,Lp,S,J],-1)
        vmat = V12mom[V12idx]
        fullrank = rank(vmat)
        tx = chara_SLJ(S,L,Lp,J) 
        SVD = LinearAlgebra.svd(vmat)
        U = SVD.U; Sig = SVD.S; Vt = SVD.Vt
        Sig[truncation_rank+1:end] .= 0.0
        SV = Diagonal(Sig)* Vt
        Vtilde = BLAS.gemm('N','N',1.0,U,SV)
        if verbose 
            println("Tz $itt $tx rank ",fullrank)
            for trank = 1:fullrank
                println("rank=",@sprintf("%4i",trank), " sval ",Sig[trank])
            end
            println("norm V-V' ",norm(vmat-Vtilde,2))
        end
        vmat .= Vtilde
    end
end


""" 
    write_onshell_vmom(chiEFTobj::ChiralEFTobject,pnrank::Int;label="")

"""
function write_onshell_vmom(chiEFTobj::ChiralEFTobject,pnrank::Int,target_LSJ;label="")
    n_mesh = chiEFTobj.params.n_mesh
    dict_pwch = chiEFTobj.dict_pwch[pnrank]
    xr = chiEFTobj.xr_fm
    V12mom = chiEFTobj.V12mom
    tx1d=""
    itt = 2 *(pnrank -2)
    for i= 1:n_mesh
        x = xr[i]
        tx = @sprintf("%18.8e", x)
        for tmp in target_LSJ
            L,Lp,S,J = tmp
            V12idx = get(dict_pwch,[itt,L,Lp,S,J],-1)
            if V12idx==-1
                tx *= @sprintf("%18.8e",0.0)
            else
                v = V12mom[V12idx][i,i]
                tx *= @sprintf("%18.8e",v)
            end
        end
        tx *= "\n"
        tx1d *= tx    
    end        
    io = open("vmom_1d_"*label*".dat","w")
    print(io,"                  ")
    for tmp in target_LSJ
        L,Lp,S,J = tmp
        ctext = string(2*S+1) * chara_L[L+1] * string(J) 
        print(io,"   "*ctext,tmp)
    end
    println(io,"")
    println(io,rstrip(tx1d))
    close(io)
    return nothing
end

"""
    print_vec(s,v;ine=false)

function to print float vectors more readable. This is usuful for debug.
"""
function print_vec(s,v,io=stdout;ine=false,long=false)
    s *= " "
    for i = 1:length(v)
        if ine
            s *= @sprintf "%9.1e" v[i]
        else
            if long
                s *= @sprintf "%15.8f" v[i]
            else  
                s *= @sprintf "%10.4f" v[i]
            end
    	end
    end
    println(io,s)
end

"""
    momplot(xr,V12mom,tdict,pnrank,llpSJ_s;ctext="",fpath="")

plot nn potential in partial wave over relative momentum space
"""
function momplot(xr,V12mom,tdict,pnrank,llpSJ_s;ctext="",fpath="")
    tfdat = []
    if fpath != ""; xf,yfs = compfdat(fpath); end    
    itt = 2 *(pnrank -2)
    tv = zeros(Float64,n_mesh)
    for vidx = 1:7
        l,lp,S,J= llpSJ_s[vidx]
        V12idx = get(tdict,[itt,l,lp,S,J],-1)
        if V12idx==-1;println("V12idx==$V12idx");continue;end
        V = V12mom[V12idx]
        tx = ""
        cS = ifelse(S==0,"1","3")
        if l==lp
            tx=ctext*cS*chara_L[l+1]*string(J)
        else
            tx =ctext*cS*chara_L[l+1]*string(J)
            tx *= "_"*cS*chara_L[lp+1]*string(J)
        end
        for i =1:n_mesh;  tv[i] = V[i,i]; end
        if fpath != ""; tfdat = [xf,yfs[vidx]];end
        pw_plt(tx,xr,V,tv,pnrank;fdat=tfdat)
    end
    return nothing
end

"""
    pw_plt(tx,xr,z,zds,pnrank;fdat=[])

    fdat: normally no need to specify. optional array to specify text file for cross check
"""
function pw_plt(tx,xr,z,zds,pnrank;fdat=[])
    tls = ["pp","pn","nn"]
    cpnrank = tls[pnrank]    
    xr *= 1.0/hc; yr = copy(xr)
    fig = plt.figure(figsize=(10,4))
    axs = [fig.add_subplot(121),fig.add_subplot(122)]
    axs[2].set_xlabel(latexstring("p ")*" [fm"*latexstring("^{-1}")*"]")
    axs[2].set_ylabel(latexstring("p'")*" [fm"*latexstring("^{-1}")*"]")
    CS = axs[2].contourf(xr, yr, z)
    fig.colorbar(CS)
    axs[1].set_xlabel(latexstring("p=p' ")*" [fm"*latexstring("^{-1}")*"]")
    axs[1].set_ylabel(latexstring("V(p,p)")*" [MeV fm"*latexstring("^3")*"]")
    axs[1].plot(xr,zds,marker="o",markersize=2)
    if fdat != []
        axs[1].plot(fdat[1],fdat[2],marker="x",markersize=2,alpha=0.4,label="Fortran")
        axs[1].legend()
    end
    axs[1].grid(color="gray",linestyle="dotted")
    plt.savefig("pic/chiEFT_"*tx*"_"*cpnrank*".pdf",pad_inches=0)
    plt.close()
end

function compfdat(inpf)
    f = open(inpf,"r");lines = readlines(f);close(f)
    xs = Float64[]
    ys = [ Float64[] for i=1:7]
    for line in lines
        tl = split(line)
        tmp = map(x->parse(Float64,x),tl)
        push!(xs,tmp[1])
        for i=2:length(tmp); push!(ys[i-1],tmp[i]);end
    end
    return xs,ys
end
