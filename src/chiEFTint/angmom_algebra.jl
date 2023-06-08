"""
    wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)

calculate Wigner's 9j symbols, all j should be given as integer (0,1,...) or halfinteger (3/2, 3//2,...)
"""
function wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)    
    s= 0.0
    xmin = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))
    xmax = min(    j1+j9 ,    j2+j6 ,    j4+j8)
    sig = j1+j2+j3+j4+j5+j6+j7+j8+j9
    if j1 == j4 && j2 == j5 && j3 == j6 && sig % 2 == 1; return 0.0;end
    if j7 == j4 && j8 == j5 && j9 == j6 && sig % 2 == 1; return 0.0;end
    if j7 == j1 && j8 == j2 && j9 == j3 && sig % 2 == 1; return 0.0;end
    if j1 == j2 && j4 == j5 && j7 == j8 && sig % 2 == 1; return 0.0;end
    if j1 == j3 && j4 == j6 && j7 == j9 && sig % 2 == 1; return 0.0;end
    if j2 == j3 && j5 == j6 && j8 == j9 && sig % 2 == 1; return 0.0;end   
    for x =xmin:xmax
        t  = (-1)^(2*x) * (2*x+1)
        t *= wigner6j(Float64,j1,j4,j7,j8,j9,x)
        t *= wigner6j(Float64,j2,j5,j8,j4,x,j6)
        t *= wigner6j(Float64,j3,j6,j9,x,j1,j2)
        s += t 
    end
    return s
end

function wigner9j_d6jint(j1::Int64,j2::Int64,j3::Int64,
                         j4::Int64,j5::Int64,j6::Int64,
                         j7::Int64,j8::Int64,j9::Int64,d6j_int)
    sig = j1+j2+j3+j4+j5+j6+j7+j8+j9
    if j1 == j4 && j2 == j5 && j3 == j6 && sig % 2 == 1; return 0.0;end
    if j7 == j4 && j8 == j5 && j9 == j6 && sig % 2 == 1; return 0.0;end
    if j7 == j1 && j8 == j2 && j9 == j3 && sig % 2 == 1; return 0.0;end
    if j1 == j2 && j4 == j5 && j7 == j8 && sig % 2 == 1; return 0.0;end
    if j1 == j3 && j4 == j6 && j7 == j9 && sig % 2 == 1; return 0.0;end
    if j2 == j3 && j5 == j6 && j8 == j9 && sig % 2 == 1; return 0.0;end
    s= 0.0
    xmin = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))
    xmax = min(abs(j1+j9),abs(j2+j6),abs(j4+j8))
    for x =xmin:2:xmax
        t  = (-1)^(x) * (x+1)
        t *= call_d6j(j1,j4,j7,j8,j9,x,d6j_int)
        t *= call_d6j(j2,j5,j8,j4,x,j6,d6j_int)
        t *= call_d6j(j3,j6,j9,x,j1,j2,d6j_int)
        s += t 
    end
    return s
end

"""
    TMtrans(chiEFTobj::ChiralEFTobject,dWS,to;writesnt=true)

Function to carry out Talmi-Mochinsky transformation for NN interaction in HO space and to write out an sinput file.
"""
function TMtrans(chiEFTobj::ChiralEFTobject,dWS,to;writesnt=true)
    V12ab = Vrel(chiEFTobj,to)
    V12ab_2n3n = Vrel(chiEFTobj,to;calc_2n3n=true) 
    params = chiEFTobj.params
    nTBME = chiEFTobj.nTBME; infos = chiEFTobj.infos; izs_ab = chiEFTobj.izs_ab
    emax = chiEFTobj.params.emax 
    nofst = 0
    nljsnt = Vector{Int64}[] 
    for temax = 0:emax 
        for l = temax%2:2:temax
            n = div(temax-l,2)
            jmin = 2*l-1
            if jmin < 1; jmin=1;end
            for j=jmin:2:2*l+1
                push!(nljsnt,[n,l,j])
                nofst += 1
            end
        end
    end 
    nljdict = Dict{Int64,Int64}()
    if writesnt 
        target_nlj = params.target_nlj
        tbme_fmt = params.tbme_fmt
        io=open(params.fn_tbme,"w")
        if tbme_fmt  == "myg"
            println(io, nTBME)
            println(io,"tza,  a, tzb,  b, tzc,  c, tzd,  d, J12,   <ab|v|cd>,   <ab|p_ip_j|cd>/mhw")
        elseif tbme_fmt=="snt"
            if length(target_nlj)!=0
                nljdict = def_sps_snt(params.emax,target_nlj)[2]
            end
            write_spes(params,io,nljsnt,nofst,nTBME,nljdict)
        elseif tbme_fmt=="snt.bin"
            if length(target_nlj)!=0
                nljdict = def_sps_snt(emax,target_nlj)[2]
            end
            write_spes(params,io,nljsnt,nofst,nTBME,nljdict;bin=true)
        end
    end
    tbmes = [ Dict{Vector{Int64},Vector{Vector{Float64}}}() for pnrank=1:3]
    tkeys = [zeros(Int64,4) for i=1:4]
    @inbounds for ich in eachindex(infos) 
        izz,ip,Jtot,ndim=infos[ich]
        pnrank = Int(div(izz,2))+2
        izs = izs_ab[ich]
        vv = zeros(Float64,ndim,ndim)
        vv_2n3n = zeros(Float64,ndim,ndim)
        @timeit to "vtrans" @inbounds @floop for i = 1:ndim
            iza,ia,izb,ib = izs[i]
            @inbounds for j = 1:i
                izc,ic,izd,id= izs[j]
                v12,v12_2n3n = vtrans(chiEFTobj,dWS,pnrank,izz,ip,Jtot,iza,ia,izb,ib,izc,ic,izd,id,nljsnt,V12ab,V12ab_2n3n,to)
                vv[i, j] = vv[j, i] = v12    
                vv_2n3n[i, j] = vv_2n3n[j, i] = v12_2n3n    
            end
        end       
        if writesnt 
            write_tbme(params,io,ndim,izs,Jtot,vv,vv_2n3n,nljsnt,nljdict,dWS,tkeys;ofst=nofst)
        else
            set_tbme(chiEFTobj,tbmes[pnrank],ndim,izs,Jtot,vv,vv_2n3n,nljsnt,nljdict,dWS,tkeys,to;ofst=nofst)
        end
    end
    if writesnt 
        close(io)
        return nothing
    else
        return tbmes
    end
end

function set_tbme(chiEFTobj::ChiralEFTobject,tbmes,ndim,izs,Jtot,vv,vv_2n3n,nljsnt,nljdict,dWS,tkeys,to;ofst=0)
    target_nlj = chiEFTobj.params.target_nlj
    @inbounds for i = 1:ndim
        iza,ia,izb,ib = izs[i]
        na,la,ja = nljsnt[ia]
        nb,lb,jb = nljsnt[ib]
        for j = 1:i
            tv = vv[i,j]
            tv_2n3n = vv_2n3n[i,j]
            izc,ic,izd,id= izs[j]
            nc,lc,jc = nljsnt[ic]
            nd,ld,jd = nljsnt[id]
            a = ifelse(iza==-1,ia,ia+ofst)
            b = ifelse(izb==-1,ib,ib+ofst)
            c = ifelse(izc==-1,ic,ic+ofst)
            d = ifelse(izd==-1,id,id+ofst)          
            if length(target_nlj) !=0
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
            if length(target_nlj) !=0
                fa = nljdict[fa]; fb = nljdict[fb]
                fc = nljdict[fc]; fd = nljdict[fd]
            end                       
            owtkey!(tkeys[1],na,la,ja,iza)
            owtkey!(tkeys[2],nb,lb,jb,izb)
            owtkey!(tkeys[3],nc,lc,jc,izc)
            owtkey!(tkeys[4],nd,ld,jd,izd)
            vpp = kinetic_tb(tkeys[1],tkeys[2],tkeys[3],tkeys[4],Jtot,dWS)            
            owtkey!(tkeys[1],fa,fb,fc,fd)
            t = get(tbmes,tkeys[1],false)
            if t == false                           
                tbmes[copy(tkeys[1])] = [ [1.0*Jtot, 0.0, tv, tv_2n3n, vpp] ]
            else
                push!(tbmes[tkeys[1]], [1.0*Jtot, 0.0, tv, tv_2n3n,vpp])
            end 
        end
    end
    return nothing
end

"""
    red_nabla_l(n1,l1,n2,l2)
returns ``b \\langle n_1,l_1|| \\nabla ||n_2,l_2 \\rangle`` in ``l``-reduced matrix element.
"""
function red_nabla_l(n1,l1,n2,l2)
    if n1 == n2 && l1 == l2+1
        return -sqrt((l2 + 1.0)*(n2 + l2 + 1.5))
    elseif n1 == n2-1 && l1 == l2+1
        return -sqrt((l2 + 1.0)*n2)
    elseif n1 == n2 && l1 == l2-1
        return -sqrt(l2*(n2 + l2 + 0.5))
    elseif n1 == n2+1 && l1==l2-1
        return -sqrt(l2*(n2 + 1.0))
    else
        return 0.0
    end
end

"""
    red_nabla_j(nlj1,nlj2,d6j,key6j) 
returns ``b \\langle j || \\nabla || j_2\\rangle``
Note that ``l_1,l_2`` in `nlj1`&`nlj2` are not doubled.
"""
function red_nabla_j(nlj1,nlj2,d6j_lj) 
    n1, l1, j1 = nlj1; n2, l2, j2 = nlj2
    if tri_check(j1//2,j2//2,1)==false;return 0.0;end                    
    if tri_check(l1,l2,1)==false;return 0.0;end                    
    if tri_check(j1//2,1//2,l1)==false;return 0.0;end   
    ret = (-1)^((3+2*l1+j2)//2) * sqrt(1.0*(j1+1)*(j2+1)) *
          call_d6j(j2,j1,2,l1*2,l2*2,1,d6j_lj) * red_nabla_l(n1,l1,n2,l2)
    return ret
end

"""
    kinetic_ob(nlj1, nlj2)

calc. kinetic one-body contribution ```\\langle j_1 |T/\\habr\\omega| j_2 \\rangle```
"""
function kinetic_ob(nlj1, nlj2)
    n1, l1, j1 = nlj1
    n2, l2, j2 = nlj2
    ret = 0.0
    if l1 != l2 || j1 != j2 || abs(n1-n2)>2 || abs(l1-l2)>2
        return float(ret)
    end
    for n =maximum([min(n1-1,n2-1), 0]): maximum([n1+1,n2+1])+1
        for l = maximum([minimum([l1-1,l2-1]), 0]) : maximum([l1+1,l2+1])+1
            ret += -0.5 / sqrt((2*l1+1)*(2*l2+1)) * (-1)^(l-l2) * red_nabla_l(n1,l1,n,l) * red_nabla_l(n,l,n2,l2)
        end
    end
    return float(ret)
end

"""
    kinetic_tb(nljtz1, nljtz2, nljtz3, nljtz4, J, dWS)

calc. kinetic two-body contribution ``\\langle j_1j_2|| -p_1 p_2/\\hbar\\omega ||j_3j_4\\rangle_J`` using preallocated 6j-symbols.
"""
function kinetic_tb(nljtz1,nljtz2,nljtz3,nljtz4,J,dWS)
    n1, l1, j1, tz1 = nljtz1;  n2, l2, j2, tz2 = nljtz2
    n3, l3, j3, tz3 = nljtz3;  n4, l4, j4, tz4 = nljtz4
    nlj1 = @view nljtz1[1:3] 
    nlj2 = @view nljtz2[1:3] 
    nlj3 = @view nljtz3[1:3] 
    nlj4 = @view nljtz4[1:3] 
    norm = 1.0; ret = 1.0
    if nljtz1 == nljtz2; norm *= 1.0/sqrt(2.0);end
    if nljtz3 == nljtz4; norm *= 1.0/sqrt(2.0);end
    t6j = t6j_2 = 0.0
    d6j_lj = dWS.d6j_lj
    if tri_check(j1/2,j3/2,1) && tri_check(j4/2,j3/2,J) && tri_check(j4/2,j2/2,1)
        t6j = call_d6j(j1,j2,J*2,j4,j3,2,d6j_lj)
    end
    if tri_check(j1/2,j4/2,1) && tri_check(j4/2,j3/2,J) && tri_check(j3/2,j2/2,1)
        t6j_2 = call_d6j(j1,j2,J*2,j3,j4,2,d6j_lj)
    end
    nabla1324 = red_nabla_j(nlj1,nlj3,d6j_lj) * red_nabla_j(nlj2,nlj4,d6j_lj)
    
    if tz1==tz2==tz3==tz4
        nabla1423 = red_nabla_j(nlj1, nlj4,d6j_lj) * red_nabla_j(nlj2,nlj3,d6j_lj)
        ret = (-1)^((j2+j3)//2+J) * t6j * nabla1324 - (-1)^((j3+j4)//2-J) * (-1)^((j2+j4)//2+J) * t6j_2 * nabla1423
    elseif ( tz1==-1 && tz2 ==1 && tz3 == -1 && tz4 ==1) ||
           ( tz1==1 && tz2 ==-1 && tz3 == 1 && tz4 ==-1)
        ret = (-1)^((j2+j3)//2+J) * t6j * nabla1324
    else
        @error "error in kinetic_tb"
    end
    return norm*ret
end

function Nfac_jj(na,la,ja,nb,lb,jb)
    s =1.0/ sqrt(1.0+delta(na,nb)*delta(la,lb)*delta(ja,jb))
    return s
end

function mydoublefactorial(n::Integer)
    if n ==0 || n==1 || n==-1; return 1.0;end
    if n < 0
        throw(DomainError(n, "n must be nonnegative"))
    end
    z = Ref{BigInt}(0)
    ccall((:__gmpz_2fac_ui, :libgmp), Cvoid, (Ref{BigInt}, UInt), z, UInt(n))
    return Float64(z[])
end

function Ghob(e1, l1, ea, la, eb, lb,dtri,dcgm0;to=nothing) 
    cg = call_dcgm0(la,lb,l1,dcgm0)
    r  = cg * sqrt( (2.0*la + 1.0) * (2.0*lb + 1.0) ) 
    r *= sqrt(dtri[get_nkey3(e1-l1,ea-la,eb-lb)])
    r *= sqrt(dtri[get_nkey3(e1+l1+1,ea+la+1,eb+lb+1)])
    return r
end 

function trinomial(n1,na,nb)
    return mydoublefactorial(n1) /( mydoublefactorial(na) * mydoublefactorial(nb))
end

function tri_check(ja,jb,jc)
    if ja+jb < jc ; return false;end
    if jb+jc < ja ; return false;end
    if jc+ja < jb ; return false;end
    if abs(ja-jb) > jc ; return false;end
    if abs(jb-jc) > ja ; return false;end
    if abs(jc-ja) > jb ; return false;end
    return true
end

"""
    HObracket(nl, ll, nr, lr, n1, l1, n2, l2, Lam, d::Float64,dWs,tkey9j,dict9j_HOB,to)
To calc. generalized harmonic oscillator brackets (HOBs), ``<<n_l\\ell_l,n_r\\ell_r:\\Lambda|n_1\\ell_1,n_2\\ell_2:\\Lambda>>_d`` from the preallocated 9j dictionary.

## Reference:
- [1] B.Buck& A.C.merchant, Nucl.Phys. A600 (1996) 387-402
- [2] G.P.Kamuntavicius et al., Nucl.Phys. A695 (2001) 191-201
"""
function HObracket_d6j(nl, ll, nr, lr, n1, l1, n2, l2, Lam, d::Float64, dtri, dcgm0, d6j, to)
    r = 0.0
    ee = 2*nl + ll
    er = 2*nr + lr
    e1 = 2*n1 + l1
    e2 = 2*n2 + l2
    if ee + er != e1 + e2; return r;end
    if !tri_check(ll, lr, Lam);return r;end
    if !tri_check(l1, l2, Lam);return r;end
    phase = (-1.0)^(n1 + n2 + nr + nl) 
    t = sqrt( ( d^(e1 - er)) / ((1.0 + d)^(e1 + e2)))
    m = min(er, e2)    
    for ed = 0:m
        eb = er - ed
        ec = e2 - ed
        ea = e1 - er + ed
        for ld = ed:-2:0
            for lb = eb:-2:0
                if tri_check(ld,lb,lr)==false;continue;end
                for lc = ec:-2:0
                    if !tri_check(ld,lc,l2) ;continue;end
                    for la = ea:-2:0
                        if !tri_check(la,lb,l1);continue;end
                        if !tri_check(la,ll,lc);continue;end
                        t9j = wigner9j_d6jint(la*2,lb*2,l1*2,lc*2,ld*2,l2*2,ll*2,lr*2,Lam*2,d6j)
                        tmp = ((-d)^ed)  * t 
                        tmp *= t9j 
                        tmp *= Ghob(e1, l1, ea, la, eb, lb, dtri, dcgm0)
                        tmp *= Ghob(e2, l2, ec, lc, ed, ld, dtri, dcgm0)
                        tmp *= Ghob(ee, ll, ea, la, ec, lc, dtri, dcgm0)
                        tmp *= Ghob(er, lr, eb, lb, ed, ld, dtri, dcgm0)
                        r += tmp
                    end
                end
            end
        end
    end
    return r * phase 
end

""" 
    vtrans(chiEFTobj,pnrank,izz,ip,Jtot,iza,ia,izb,ib,izc,ic,izd,id,nljsnt,V12ab,t5v,to)

Function to calculate V in pn-formalism:
```math
\\langle ab;JTz|V| cd;JTz \\rangle = N_{ab} N_{cd}  \\sum_{\\Lambda S \\Lambda' S'} \\sum_{n\\ell N L}\\sum_{n'\\ell' N' L'}  
\\sum_{J_\\mathrm{rel}J'_\\mathrm{rel}} [\\Lambda][\\Lambda'] \\hat{S}\\hat{S'} \\hat{J}_\\mathrm{rel}\\hat{J}'_\\mathrm{rel} 
\\hat{j}_a \\hat{j}_b \\hat{j}_c \\hat{j}_d
(-1)^{\\ell + S + J_\\mathrm{rel} + L} (-1)^{\\ell' + S' + J'_\\mathrm{rel} + L'}
\\langle n N [ (\\ell L)\\Lambda S] J| n_a n_b [ (\\ell_a \\ell_b)\\Lambda (\\tfrac{1}{2}\\tfrac{1}{2})S ]J \\rangle_{d=1}
\\langle n' N' [ (\\ell' L')\\Lambda' S'] J| n_c n_d [ (\\ell_c \\ell_d)\\Lambda' (\\tfrac{1}{2}\\tfrac{1}{2})S' ]J \\rangle_{d=1}
\\left\\{ \\begin{matrix} \\ell_a & \\ell_b & \\Lambda \\\\ 1/2 & 1/2 & S \\\\ j_a & j_b & J \\end{matrix} \\right\\} 
\\left\\{ \\begin{matrix} \\ell_c & \\ell_d & \\Lambda' \\\\ 1/2 & 1/2 & S' \\\\ j_c & j_d & J \\end{matrix} \\right\\} 
\\left\\{ \\begin{matrix} L & \\ell & \\Lambda \\\\ S & J & J_\\mathrm{rel} \\end{matrix} \\right\\}
\\left\\{ \\begin{matrix} L' & \\ell' & \\Lambda' \\\\ S' & J & J'_\\mathrm{rel} \\end{matrix} \\right\\} 
\\langle n\\ell S J_\\mathrm{rel} T|V_\\mathrm{NN}|n'\\ell' S' J'_\\mathrm{rel} T\\rangle
```
"""
function vtrans(chiEFTobj::ChiralEFTobject,dWS,pnrank,izz,ip,Jtot,iza,ia,izb,ib,izc,ic,izd,id,nljsnt,V12ab,V12ab_2n3n,to)
    d6j_int = dWS.d6j_int
    d9j_lsj = dWS.d9j_lsj
    arr_pwch = chiEFTobj.arr_pwch
    ret = ret_2n3n = 0.0
    na,la,jda = nljsnt[ia]; nb,lb,jdb = nljsnt[ib]
    nc,lc,jdc = nljsnt[ic]; nd,ld,jdd = nljsnt[id]
    lrmax = jmax + 1
    Eab=2*na+la+2*nb+lb; Ecd=2*nc+lc+2*nd+ld
    TF = false
    if izz != iza+izb || izz != izc+izd; TF=true;end
    if Jtot > div(jda+jdb,2) || Jtot < abs(div(jda-jdb,2)); TF=true;end
    if Jtot > div(jdc+jdd,2) || Jtot < abs(div(jdc-jdd,2)); TF=true;end
    if (-1)^(la+lb) != ip || (-1)^(lc+ld) != ip; TF=true;end
    if (izz==2 || izz==-2) && ia==ib && (-1)^(Jtot)==-1; TF=true;end
    if (izz==2 || izz==-2) && ic==id && (-1)^(Jtot)==-1; TF=true;end
    if TF; return ret,ret_2n3n;end   
    for S=0:1
        tarr_pwch = arr_pwch[pnrank][S+1]
        lmax1=min(Jtot+S,la+lb); lmin1=max(abs(Jtot-S),abs(la-lb)); if lmin1 > lmax1;continue;end
        lmax2=min(Jtot+S,lc+ld); lmin2=max(abs(Jtot-S),abs(lc-ld)); if lmin2 > lmax2;continue;end
        @inbounds for Lam=lmin1:lmax1
            x1 = sqrt((jda+1.0)*(jdb+1.0)) * hat(Lam) *hat(S) * (-1)^Lam
            x1 *= call_d9j_lsj(la*2,lb*2,Lam*2,1,1,S*2,jda,jdb,Jtot*2,d9j_lsj)     
            if x1 == 0.0; continue;end
            for Lamp=lmin2:lmax2
                x2 = sqrt((jdc+1.0)*(jdd+1.0)) * hat(Lamp) *hat(S) * (-1)^Lamp
                x2 *= call_d9j_lsj(lc*2,ld*2,Lamp*2,1,1,S*2,jdc,jdd,Jtot*2,d9j_lsj)  
                if x2 == 0.0; continue;end
                for Ncm=0:div(Eab,2)
                    Lcm_max =min(Eab-2*Ncm, Ecd-2*Ncm)
                    if Lcm_max < 0;continue;end
                    for Lcm=0:Lcm_max
                        for lr1=abs(Lcm-Lam):min(lrmax,Lcm+Lam)
                            nx1=Eab-2*Ncm-(lr1+Lcm); nr1=div(nx1,2); if nr1 < 0 || nx1!=2*nr1;continue;end
                            y1 = get_dictHOB(nr1,lr1,Ncm,Lcm,na,la,nb,lb,Lam,dWS.dictHOB)
                            if abs(y1) < 1.e-10;continue;end
                            for lr2=abs(Lcm-Lamp):Lcm+Lamp
                                if lr1%2 != lr2%2;continue;end
                                if lr2 > lrmax;continue;end
                                nx2=Ecd-2*Ncm-(lr2+Lcm); nr2=div(nx2,2)
                                if nr2 < 0 || (nx2!=2*nr2);continue;end
                                y2 = get_dictHOB(nr2,lr2,Ncm,Lcm,nc,lc,nd,ld,Lamp,dWS.dictHOB)
                                if abs(y2) < 1.e-10;continue;end

                                mj1=abs(lr1-S); mj2=abs(lr2-S); mj3=abs(Jtot-Lcm)
                                kjmin=max(mj1,mj2,mj3)
                                kjmax=min(lr1+S,lr2+S,Jtot+Lcm,6)
                                if kjmin > kjmax;continue;end
                                sumv=sumv_2n3n=0.0
                                @inbounds for Jrel=kjmin:kjmax
                                    zu1  = hat(Lam)*hat(Jrel)*(-1.0)^(lr1+S+Jrel+Lcm)
                                    zu1 *= call_d6j_nond(Lcm,lr1,Lam,S,Jtot,Jrel,d6j_int)
                                    if abs(zu1) < 1.e-10;continue;end

                                    zu2  = hat(Lamp)*hat(Jrel)*(-1.0)^(lr2+S+Jrel+Lcm)
                                    zu2 *= call_d6j_nond(Lcm,lr2,Lamp,S,Jtot,Jrel,d6j_int)
                                    if abs(zu2) < 1.e-10;continue;end

                                    izfac=0
                                    if izz==-2 || izz==2
                                        izfac=1+(-1)^(lr1+S)
                                    end
                                    if izz==0;izfac=1;end
                                    rv12=rv12_2n3n=0.0
                                    if izfac!=0
                                        num= tarr_pwch[Jrel+1][lr1-abs(Jrel-S)+1][lr2-abs(Jrel-S)+1]
                                        rv12=V12ab[num][nr1+1,nr2+1]
                                        rv12_2n3n=V12ab_2n3n[num][nr1+1,nr2+1]
                                    end
                                    sumv += zu1*zu2*izfac*rv12
                                    sumv_2n3n += zu1*zu2*izfac*rv12_2n3n
                                end
                                zxy=x1*x2*y1*y2
                                ret += sumv*zxy
                                ret_2n3n += sumv_2n3n*zxy
                            end
                        end
                    end
                end
            end
        end
    end    
    if pnrank!= 2
        Nab = Nfac_jj(na,la,jda,nb,lb,jdb)
        Ncd = Nfac_jj(nc,lc,jdc,nd,ld,jdd)
        ret *= Nab*Ncd
        ret_2n3n *= Nab*Ncd
    end
    return ret,ret_2n3n
end

"""
    prep_wsyms()

preparing Clebsch-Gordan coefficients for some special cases: cg1s = (1,0,l,0|l',0), cg2s = (2,0,l,0|l',0)
"""
function prep_wsyms(;lmax=7)
    dim = lmax+3
    cg1s = zeros(Float64,dim,dim)
    cg2s = zeros(Float64,dim,dim)
    d6_121 = zeros(Float64,dim,dim,dim)
    d6_21 = zeros(Float64,dim,dim,dim,dim)
    d6_222 = zeros(Float64,dim,dim,dim)
    d9_12 = zeros(Float64,dim,dim,dim,dim)
    for li = 0:dim-1
        for lj =0:dim-1
            cg1s[li+1,lj+1] = clebschgordan(Float64,1,0,li,0,lj,0)
            cg2s[li+1,lj+1] = clebschgordan(Float64,2,0,li,0,lj,0)
            for j = 0:dim-1
                d6_121[li+1,lj+1,j+1] = wigner6j(Float64,1,2,1,li,lj,j)
                d6_222[li+1,lj+1,j+1] = wigner6j(Float64,2,li,lj,j,2,2)
                for jp = 0:dim-1
                    d6_21[li+1,lj+1,j+1,jp+1] = wigner6j(Float64,li,lj,2,j,jp,1)
                    d9_12[li+1,lj+1,j+1,jp+1] = wigner9j(1,li,lj,1,j,jp,2,2,2)
                end
            end
        end
    end
    return wsyms_j1_1or2(cg1s,cg2s,d6_121,d6_21,d6_222,d9_12)
end

# function jj_std(sps,dictsps,dictTBMEs;fname="")
#     if fname=="";fname = "monopole"*fname*".dat"
#     else;fname = "monopole_"*fname*".dat";end
#     ln = length(sps)
#     nmax = maximum([sps[i][1] for i=1:ln])
#     nmin = minimum([sps[i][1] for i=1:ln])
#     lmax = maximum([sps[i][2] for i=1:ln])
#     lmin = minimum([sps[i][2] for i=1:ln])
#     jmax = maximum([sps[i][3] for i=1:ln])
#     jmin = minimum([sps[i][3] for i=1:ln])
#     emax = maximum( [ 2*sps[i][1]+sps[i][2] for i=1:length(sps)])
#     emin = minimum( [ 2*sps[i][1]+sps[i][2] for i=1:length(sps)])
#     Jmax = div(2*jmax+1,2)
#     Lmin = 0; Lmax = 2*lmax
#     monodict = Dict(["",""]=> [[0.0,[0.0]]]);delete!(monodict,["",""])
#     cpnrank=["pp","pn","nn"]

#     keya = [0,0,0,0]; keyb = [0,0,0,0];keyc=[0,0,0,0]; keyd=[0,0,0,0]   # nljtz
#     key_abcdJ = [0,0,0,0,0]
#     tvjj = zeros(Float64,5)
#     for pnrank =1:3
#         tdict = dictTBMEs[pnrank]
#         for tkey in keys(tdict)
#             a,b,c,d,J = tkey
#             na,la,ja2,iza = sps[a]; nb,lb,jb2,izb = sps[b]
#             nc,lc,jc2,izc = sps[c]; nd,ld,jd2,izd = sps[d]
#             keya[1]=na; keya[2]=la; keya[3]=ja2; keya[4]=iza
#             keyb[1]=nb; keyb[2]=lb; keyb[3]=jb2; keyb[4]=izb
#             keyc[1]=nc; keyc[2]=lc; keyc[3]=jc2; keyc[4]=izc
#             keyd[1]=nd; keyd[2]=ld; keyd[3]=jd2; keyd[4]=izd
#             v_target = tdict[tkey]            
#             tvjj .= 0.0
#             hats = (-1)^J * sqrt((ja2+1.0)*(jb2+1.0)*(jc2+1.0)*(jd2+1.0))
#             for L = abs(la-lb):la+lb
#                 for S=0:1
#                     for Lp = abs(lc-ld):lc+ld
#                         for Sp = 0:1
#                             Jfac  = wigner9j(la,1//2,ja2//2,lb,1//2,jb2//2,L,S,J)
#                             Jfac *= wigner9j(lc,1//2,jc2//2,ld,1//2,jd2//2,Lp,Sp,J)
#                             for jpa2 = 2*la-1:2:2*la+1
#                                 if jpa2 <=0;continue;end
#                                 for jpb2 = 2*lb-1:2:2*lb+1
#                                     if jpb2 <=0;continue;end
#                                     for jpc2 = 2*lc-1:2:2*lc+1
#                                         if jpc2 <=0;continue;end
#                                         for jpd2 = 2*ld-1:2:2*ld+1
#                                             if jpd2 <=0;continue;end
#                                             hatps = sqrt((jpa2+1.0)*(jpb2+1.0)*(jpc2+1.0)*(jpd2+1.0))
#                                             LLSS = (2*L+1)*(2*Lp+1)*(2*S+1)*(2*Sp+1)
#                                             for Jp = max(abs(Lp-Sp),abs(L-S)):min(Lp+Sp,L+S)
#                                                 d9j_ab = wigner9j(la,1//2,jpa2//2,lb,1//2,jpb2//2,L,S,Jp)
#                                                 d9j_cd = wigner9j(lc,1//2,jpc2//2,ld,1//2,jpd2//2,Lp,Sp,Jp)
#                                                 keya[3] = jpa2; keyb[3]=jpb2;keyc[3] = jpc2; keyd[3]=jpd2                                                
#                                                 ap = dictsps[keya];bp = dictsps[keyb];cp = dictsps[keyc];dp = dictsps[keyd]
#                                                 key_abcdJ[1] = ap;key_abcdJ[2] = bp
#                                                 key_abcdJ[3] = cp;key_abcdJ[4] = dp;key_abcdJ[5] = Jp
#                                                 ta=ap; tb=bp;tc=cp;td=dp
#                                                 phase = 1.0
#                                                 if ap > bp
#                                                     ta = bp; tb = ap
#                                                     phase *= (-1)^(div(sps[ap][3]+sps[bp][3],2)+Jp+1)
#                                                 end
#                                                 if cp > dp
#                                                     tc = dp; td = cp
#                                                     phase *= (-1)^(div(sps[cp][3]+sps[dp][3],2)+Jp+1)
#                                                 end
#                                                 fa=ta; fb=tb;fc=tc;fd=td
#                                                 if fa > fc || (fa==fc && fb > fd)
#                                                     fa=tc;fb=td;fc=ta;fd=tb
#                                                 end
#                                                 key_abcdJ[1] = fa;key_abcdJ[2] = fb
#                                                 key_abcdJ[3] = fc;key_abcdJ[4] = fd;key_abcdJ[5] = Jp
#                                                 vjj = get(tdict,key_abcdJ,false)                                                
#                                                 if vjj==false;continue; end
#                                                 #println("abcd' $ap $bp $cp $dp J' $Jp $nvjj $vjj")
#                                                 Jpfac = (-1)^Jp * (2*Jp+1) * d9j_ab * d9j_cd
#                                                 for k = 0:2 
#                                                     tidx = k+1; if S==Sp && k==1;tidx=4;end                                                
#                                                     kfac = (2*k+1) * wigner6j(Float64,L,S,J,Sp,Lp,k) * wigner6j(Float64,L,S,Jp,Sp,Lp,k)
#                                                     tvjj[tidx] += hats*hatps*kfac*Jfac*Jpfac * LLSS * vjj *phase
#                                                 end
#                                             end
#                                         end
#                                     end
#                                 end                                
#                             end
#                         end
#                     end
#                 end
#             end
#             tvjj[5] = sum( @views tvjj[1:4])
#             if abs(v_target-tvjj[5]) > 1.e-9
#                 println("Error TBME(in jj) mismatch!: $tkey  v ",
#                         @sprintf("%12.4e",v_target),"; sum ", @sprintf("%12.4e",tvjj[5]))
#             end
#             if a==c && b==d && ( a <= b )
#                 c = cpnrank[pnrank]
#                 ca = string(na)*chara_l[la+1]*string(ja2)
#                 cb = string(nb)*chara_l[lb+1]*string(jb2)
#                 if get(monodict,[c,ca,cb],false) == false
#                     monodict[[c,ca,cb]] = [ [J*1.0,copy(tvjj)] ]
#                 else
#                     push!(monodict[[c,ca,cb]],[J*1.0,copy(tvjj)])
#                 end
#             end
#         end
#         monopole(monodict,fname)
#     end
#     return nothing
# end

function monopole(monodict,fname)
    io = open(fname,"w")
    nume = zeros(Float64,5)
    deno = 0.0
    for key in keys(monodict)
        tarr = monodict[key]
        nume .= 0.0
        deno = 0.0
        for tmp in tarr
            J, V =tmp
            nume .+= (2*J+1.0) .* V
            deno += 2*J +1.0
        end
        vm = nume ./ deno
        tx = key[1]*":"*key[2]*"-"*key[3]*@sprintf("%20.12f", vm[1])*@sprintf("%20.12f", vm[2])*@sprintf("%20.12f", vm[3])*@sprintf("%20.12f", vm[4])*@sprintf("%20.12f", vm[5])
        println(io,tx)
    end
    close(io)
    return nothing
end

function nlj_from_str(str)
    lstr = match(r"[a-z]",str).match
    n,j = split(str,lstr)
    l = -1
    for k=1:length(chara_L)
        if lstr==chara_l[k];l=k;break;end
    end    
    l = l-1
    n = parse(Int64,n)
    j = parse(Int64,j)
    return n,l,j   
end

function get_nkey2_arr(tmp;ofst=10^3)
    return tmp[1] + tmp[2] * ofst
end
function get_nkey2(i,j;ofst=10^3)
    return i + j * ofst
end
function get_nkey2_u(i,j;ofst=10^3)::UInt64
    return UInt64(i + j * ofst)
end

function get_nkey3(i,j,k;ofst=10^3)
    return i + ofst * j + ofst^2 * k
end

function get_nkey3_JPT(arr::Vector{Int64};ofst=10^3)
    aTz,P,J = arr
    return UInt64(aTz + ofst *(P+2) + ofst^2 * J)
end
function get_nkey3_JPT(aTz,P,J;ofst=10^3)
    return UInt64(aTz + ofst *(P+2) + ofst^2 * J)
end

function get_nkey4(i,j,k,l;ofst=10^3)
    return i + ofst * j + ofst^2 * k + ofst^3 * l
end

function get_nkey4_ketJT(i,j,J,Tz;ofst=10^3)::UInt64
    return UInt64(i + ofst * j + ofst^2 * J + ofst^3 * (Tz+2))
end

function get_nkey6(i,j,k,l,m,n;ofst=10^3)::UInt64
    return UInt64(i + ofst * j + ofst^2 * k + ofst^3 * l+ ofst^4 * m + ofst^5 * n)
end

"""
    get_nkey_from_abcdarr(tkey;ofst=1000)

To get integer key from an Int array (with length greater than equal 4)
"""
function get_nkey_from_abcdarr(tkey;ofst=1000)
    return tkey[1] + tkey[2] * ofst + tkey[3] * ofst^2 + tkey[4] * ofst^3 
end

function call_dcgm0(l1,l2,l,dcgm0)
    l1p = l1 
    l2p = l2
    l3p = l
    fac = 1.0
    lmin = min(l1p,min(l2p,l3p))
    if l1p != lmin
        if l2p == lmin
            l1p,l2p = l2p,l1p
            fac *= (-1)^(l1+l2-l)
        elseif l3p == lmin
            fac *= hat(l3p) / hat(l1p) * (-1)^l2p
            l1p,l3p = l3p,l1p    
        end
    end
    if l2p > l3p
        fac *= hat(l3p) / hat(l2p) * (-1)^l1p
        l2p,l3p = l3p,l2p
    end
    return fac * dcgm0[get_nkey3(l1p,l2p,l3p)]    
end

function call_d6j_nond(j1,j2,j3,j4,j5,j6,d6j)
    return d6j[get_key6j_sym(j1*2,j2*2,j3*2,j4*2,j5*2,j6*2)]
end

function get_phase_j3(j1::Int64,j2::Int64,j3::Int64)::Float64
    p = div(j1+j2+j3,2) % 2
    return ifelse(p==1,-1.0,1.0)
end

function call_d6j(j1::Int64,j2::Int64,j3::Int64,j4::Int64,j5::Int64,j6::Int64,d6j::Dict{UInt64,Float64})::Float64
    if j1 == 0; return delta(j5,j6) * delta(j2,j3) * get_phase_j3(j4,j5,j3) / sqrt((j3+1)*(j6+1)); end
    if j2 == 0; return delta(j1,j3) * delta(j4,j6) * get_phase_j3(j5,j1,j6) / sqrt((j1+1)*(j4+1)); end
    if j3 == 0; return delta(j1,j2) * delta(j4,j5) * get_phase_j3(j2,j4,j6) / sqrt((j2+1)*(j5+1)); end
    if j4 == 0; return delta(j3,j5) * delta(j2,j6) * get_phase_j3(j1,j2,j3) / sqrt((j3+1)*(j6+1)); end
    if j5 == 0; return delta(j3,j4) * delta(j1,j6) * get_phase_j3(j1,j2,j3) / sqrt((j1+1)*(j4+1)); end
    if j6 == 0; return delta(j1,j5) * delta(j2,j4) * get_phase_j3(j1,j2,j3) / sqrt((j2+1)*(j5+1)); end
    nkey = get_key6j_sym(j1,j2,j3,j4,j5,j6)
    return d6j[nkey]
end

function j_col_score(j1::Int64,j2::Int64)::Int64
    100 * (j1 + j2) + min(j1,j2)
end

function get_canonical_order_6j(j1::Int64,j2::Int64,j3::Int64,j4::Int64,j5::Int64,j6::Int64)::Int64
    o = 0
    cscore_12 = j_col_score(j1,j2)
    cscore_34 = j_col_score(j3,j4)
    cscore_56 = j_col_score(j5,j6)
    if cscore_12 <= cscore_34 #12 is fixed => 312, 132, 123
        if cscore_34 <= cscore_56 
            o = 123            
        elseif cscore_56 < cscore_12 # 3 < 12
            o = 312
        else
            o = 132
        end
    else # 2 < 1 is fixed 
        if cscore_56 <= cscore_34
            o = 321
        elseif cscore_12 <= cscore_56
            o = 213
        else
            o = 231
        end
    end
    return o
end

function get_key6j_sym(j1::Int64,j3::Int64,j5::Int64,j2::Int64,j4::Int64,j6::Int64)::UInt64
    tj1 = j1; tj3 = j3; tj5 = j5
    tj2 = j2; tj4 = j4; tj6 = j6
    column_order = get_canonical_order_6j(j1,j2,j3,j4,j5,j6)
    if column_order == 231
        tj1,tj2,tj3,tj4,tj5,tj6 = tj3,tj4,tj5,tj6,tj1,tj2
    elseif column_order == 132
        tj3,tj4,tj5,tj6 = tj5,tj6,tj3,tj4
    elseif column_order == 213
        tj1,tj2,tj3,tj4 = tj3,tj4,tj1,tj2
    elseif column_order == 312
        tj1,tj2,tj3,tj4,tj5,tj6 = tj5,tj6,tj1,tj2,tj3,tj4
    elseif column_order == 321
        tj1,tj2,tj5,tj6 = tj5,tj6,tj1,tj2
    end
 
    # if ju = jl in any column, we can change orders to satisfy ju<=jl for all the columns with no additional factor
    if (tj1 == tj2 || tj3 == tj4 || tj5 == tj6)
        if tj1 > tj2 
            tj1,tj2=tj2,tj1
        end
        if tj3 > tj4
            tj3,tj4=tj4,tj3
        end
        if tj5 > tj6
            tj5,tj6=tj6,tj5
        end
        return get_nkey6(tj1,tj3,tj5,tj2,tj4,tj6)
    else 
        tint = ifelse(tj1<tj2,1,0) + ifelse(tj3<tj4,1,0) + ifelse(tj5<tj6,1,0)
        if tint == 0 || tint == 3 
            return get_nkey6(tj1,tj3,tj5,tj2,tj4,tj6)
        end
        if (tj5 < tj6 && tint == 1) || (tj5 > tj6 && tint == 2)            
            #tj1,tj2=tj2,tj1; tj3,tj4=tj4,tj3
            return get_nkey6(tj2,tj4,tj5,tj1,tj3,tj6)
        elseif (tj3 < tj4 && tint == 1) || ( tj3 > tj4 && tint == 2)
            #tj1,tj2=tj2,tj1;tj5,tj6=tj6,tj5
            return get_nkey6(tj2,tj3,tj6,tj1,tj4,tj5)
        elseif (tj1 < tj2 && tint == 1) || (tj1 > tj2 && tint == 2)
            #tj3,tj4=tj4,tj3;tj5,tj6=tj6,tj5
            return get_nkey6(tj1,tj4,tj6,tj2,tj3,tj5)
        end
    end
end

function call_d9j_int(j1,j2,j3,j4,j5,j6,j7,j8,j9,d9j)
    key1 = get_nkey3(j1,j2,j3)
    key2 = get_nkey3(j4,j5,j6)
    key3 = get_nkey3(j7,j8,j9)
    return d9j[key1][key2][key3]
end

function call_d9j_int_notdoubled(j1,j2,j3,j4,j5,j6,j7,j8,j9,d9j)
    key1 = get_nkey3(j1*2,j2*2,j3*2)
    key2 = get_nkey3(j4*2,j5*2,j6*2)
    key3 = get_nkey3(j7*2,j8*2,j9*2)
    return d9j[key1][key2][key3]
end

function get_key9j_lsj(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    key1 = key2 = key3 = -1
    # check columns order
    if j1%2==j4%2==j7%2==0
        key1 = get_nkey3(j1,j2,j3)
        key2 = get_nkey3(j4,j5,j6)
        key3 = get_nkey3(j7,j8,j9)        
    elseif j2%2==j5%2==j8%2==0
        key1 = get_nkey3(j2,j3,j1)
        key2 = get_nkey3(j5,j6,j4)
        key3 = get_nkey3(j8,j9,j7)
    elseif j3%2==j6%2==j9%2==0
        key1 = get_nkey3(j3,j1,j2)
        key2 = get_nkey3(j6,j4,j5)
        key3 = get_nkey3(j9,j7,j8)
    else
        @error "This case must not happen in get_keyt9j_lsj\n {$j1 $j2 $j3\n $j4 $j5 $j6\n $j7 $j8 $j9}"
    end
    # check rows order
    tkey1 = key1; tkey2 = key2; tkey3 = key3
    if (j1%2==j2%2==j3%2==0); key1 = tkey1; key2 = tkey2; key3 = tkey3;end 
    if (j4%2==j5%2==j6%2==0); key1 = tkey2; key2 = tkey3; key3 = tkey1;end
    if (j7%2==j8%2==j9%2==0); key1 = tkey3; key2 = tkey1; key3 = tkey2;end
    return key1,key2,key3
end

function call_d9j_lsj(j1,j2,j3,j4,j5,j6,j7,j8,j9,d9j)
    key1,key2,key3 = get_key9j_lsj(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    return d9j[key1][key2][key3]
end

"""
`dWS2n` struct 
"""
function prep_dWS2n(params,to;emax_calc=0)
    emax = ifelse(emax_calc!=0,emax_calc,params.emax)
    e2max = 2 * emax
    Nnmax = params.Nnmax
    jmax = 2*emax + 1
    lmax = max(e2max,Nnmax,jmax) * 2
    dtri = prep_dtri(lmax+1)
    dcgm0 = prep_dcgm0(lmax)
    @timeit to "d6j_int" d6j_int = prep_d6j_int(emax,jmax,to)
    @timeit to "d6j_lj" d6j_lj  = prep_d6j_lj(jmax)
    @timeit to "d9j" d9j_lsj = prep_d9j_lsj(jmax,lmax)
    @timeit to "HOB" dictHOB = prep_dictHOB(e2max,dtri,dcgm0,d6j_int,to)

    println("size of dWS (jmax $jmax lmax $lmax e2max $e2max Nnmax $Nnmax):\n",
            show_size_inMB("   dtri",dtri), show_size_inMB("  dcgm0",dcgm0),
            show_size_inMB("d6j_int",d6j_int),show_size_inMB(" d6j_lj",d6j_lj),"\n",
            show_size_inMB("d9j_lsj",d9j_lsj),show_size_inMB("dictHOB",dictHOB))            
    return dWS2n(dtri,dcgm0,d6j_int,d6j_lj,d9j_lsj,dictHOB)
end

#{la 1/2 ja;lb 1/2 jb;L S J}
function prep_d9j_lsj(jmax2,Jmax)
    d9j_lsj = Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}()
    for J = 0:2:Jmax
        for S = 0:2:2
            for L = abs(J-S):2:J+S
               for ja = 1:2:jmax2
                    for la = abs(ja-1):2:ja+1
                        for jb = abs(ja-J):2:ja+J
                            for lb = abs(jb-1):2:jb+1
                                t9j = wigner9j(la/2,1/2,ja/2,lb/2,1/2,jb/2,L/2,S/2,J/2)
                                key1,key2,key3= get_key9j_lsj(la,lb,L,1,1,S,ja,jb,J)
                                if !haskey(d9j_lsj,key1); d9j_lsj[key1] =Dict{Int64,Dict{Int64,Float64}}();end
                                if !haskey(d9j_lsj[key1],key2); d9j_lsj[key1][key2] = Dict{Int64,Float64}();end
                                d9j_lsj[key1][key2][key3]  = t9j
                            end
                        end
                    end                
                end
            end
        end
    end
    return d9j_lsj
end

function prep_dictHOB(e2max,dtri,dcgm0,d6j_int,to)
    dictHOB = Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}()
    for Lam = 0:e2max
        dictHOB[Lam] = Dict{Int64,Dict{Int64,Float64}}()
    end
    for E = 0:e2max
        for e_nl = 0:E
            e_NL = E - e_nl
            for n = 0:div(e_nl,2)
                l = e_nl - 2*n
                for N = 0:div(e_NL,2)
                    L = e_NL - 2*N                           
                    for e_nlp = 0:E
                        e_NLp = E - e_nlp
                        for np = 0:div(e_nlp,2)
                            lp = e_nlp - 2*np
                            for Np = 0:div(e_NLp,2)
                                Lp = e_NLp - 2*Np
                                for Lam = max(abs(l-L),abs(lp-Lp)):min(l+L,lp+Lp)
                                    key1, key2 = get_HOB_nlkey(n, l, N, L, np, lp, Np, Lp)
                                    if !haskey(dictHOB[Lam],key1); dictHOB[Lam][key1] = Dict{Int64,Float64}();end
                                    tHOB = HObracket_d6j(n, l, N, L, np, lp, Np, Lp, Lam, 1.0, dtri, dcgm0, d6j_int, to)
                                    dictHOB[Lam][key1][key2] = tHOB * (-1)^(l+Lp)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    return dictHOB
end
function check0_wigner6j(j1,j2,j3,j4,j5,j6)
    if j1*j2*j3*j4*j5*j6 ==0
        if j1 == 0; return delta(j5,j6) * delta(j2,j3) * (-1)^((j4+j5+j3)/2) / sqrt((j3+1)*(j6+1)); end
        if j2 == 0; return delta(j1,j3) * delta(j4,j6) * (-1)^((j5+j1+j6)/2) / sqrt((j1+1)*(j4+1)); end
        if j3 == 0; return delta(j1,j2) * delta(j4,j5) * (-1)^((j2+j4+j6)/2) / sqrt((j2+1)*(j5+1)); end
        if j4 == 0; return delta(j3,j5) * delta(j2,j6) * (-1)^((j1+j2+j3)/2) / sqrt((j3+1)*(j6+1)); end
        if j5 == 0; return delta(j3,j4) * delta(j1,j6) * (-1)^((j1+j2+j3)/2) / sqrt((j1+1)*(j4+1)); end
        if j6 == 0; return delta(j1,j5) * delta(j2,j4) * (-1)^((j1+j2+j3)/2) / sqrt((j2+1)*(j5+1)); end
    else
         return wigner6j(Float64,j1/2,j2/2,j3/2,j4/2,j5/2,j6/2)
    end
end

function prep_d6j_int(emax,jmax_in,to)
    jmax = jmax_in * 2
    d6j_int = Dict{UInt64,Float64}()
    for J12 = 0:2:jmax
        for j1 = 0:2:jmax_in
            for j2 = abs(J12-j1):2:j1+J12
                for j3 = 0:2:jmax
                    for J23 = abs(j2-j3):2:j2+j3
                        for J = abs(J23-j1):2:J23+j1
                            if !tri_check(J/2,J12/2,j3/2); continue;end
                            if !(j1+j3<=j2+J<=J12+J23); continue;end
                            #tj1, tj2, tj3, tj4, tj5, tj6 = get_key6j_sym(j1,j2,J12,j3,J,J23)
                            #d6j_int[get_nkey6(tj1,tj2,tj3,tj4,tj5,tj6)] = wigner6j(Float64,j1/2,j2/2,J12/2,j3/2,J/2,J23/2)
                            nkey = get_key6j_sym(j1,j2,J12,j3,J,J23)
                            d6j_int[nkey] = wigner6j(Float64,j1/2,j2/2,J12/2,j3/2,J/2,J23/2)
                        end
                    end
                end
            end
        end
    end
    return d6j_int
end

function prep_d6j_lj(jmax2)
    d6j_lj = Dict{UInt64,Float64}()
    for j1 = 1:2:jmax2
        for j2 = 1:2:jmax2
            for J12 = abs(j1-j2):2:j1+j2
                for j3 = 1:2:jmax2
                    for J23 = abs(j2-j3):2:j2+j3
                        for J = max(abs(j1-J23),abs(j3-J12)):2:min(j1+J23,j3+J12)
                            nkey = get_key6j_sym(j1,j2,J12,j3,J,J23)
                            d6j_lj[nkey] =  wigner6j(Float64,j1/2,j2/2,J12/2,j3/2,J/2,J23/2)
                        end
                    end
                end
            end
        end
    end
    # some special case for kinetic_tb
    for j2 = 1:2:jmax2
        J12 = 1*2
        J23 = 1
        for j1 = abs(j2-J12):2:j2+J12
            for l1 = abs(j1-1):2:j1+1
                for l2 = abs(j2-1):2:(j2+1)
                    nkey = get_key6j_sym(j2,j1,J12,l1,l2,J23)
                    d6j_lj[nkey] = wigner6j(Float64,j2/2,j1/2,J12/2,l1/2,l2/2,J23/2)
                end
            end
        end
    end
    return d6j_lj
end

function prep_dtri(lmax)
    dtri = Dict{Int64,Float64}()
    for l1 = 0:lmax
        for l2 = 0:lmax
            for l = 0:lmax
                dtri[get_nkey3(l1,l2,l)] = trinomial(l1,l2,l)
            end
        end
    end
    return dtri
end

function prep_dcgm0(lmax)
    dcgm0 = Dict{Int64,Float64}()
    for l1 = 0:lmax
        for l2 = l1:lmax
            for l = abs(l1-l2):l1+l2                    
                if l1 <= l2 <= l
                    tkey = get_nkey3(l1,l2,l)
                    tcg = clebschgordan(Float64,l1,0,l2,0,l,0)
                    dcgm0[tkey] = tcg
                end
            end
        end
    end
    return dcgm0
end

"""
it doesn't work for now
"""
function prep_dWS3N(N3max,J12max,j3max,to;debug=false)
    jmax_9j = lmax_6j = 2 * N3max
    lmax_9j = jmax_9j 
    lmax = jmax_9j * 2
    println("N3max $N3max J12max $J12max jmax_9j $jmax_9j ")  
    @timeit to "trinomial&cgm0" begin
        dtri = prep_dtri(lmax)
        dcgm0 = prep_dcgm0(lmax)
    end    
    #  {j1, j2, j12; j3 J j23}
    @timeit to "6j" begin
        hit6j_int = hit6j_lj = 0      
        @timeit to "const:d6j_int" begin
            # {l' l X; j j' s} type
            d6j_int = Dict{UInt64,Float64}()
            for lp = 0:2:lmax_6j
                for l = 0:2:lmax_6j
                    for s = 0:2:2
                        for j = abs(l-s):l+s
                            for jp = abs(lp-s):lp+s
                                for X = abs(l-lp):2:l+lp
                                    nkey = get_key6j_sym(lp,l,X,j,jp,s)
                                    d6j_int[nkey] = wigner6j(Float64,lp//2,l//2,X//2,j//2,jp//2,s//2)
                                    hit6j_int += 1
                                end
                            end
                        end
                    end
                end
            end
            # {L' K1 X; K2 L K3}-type (L,L',X can be large and K's are at most 2)
            for Lp = 0:2:lmax
                for L = 0:2:lmax
                    for K3 = 0:2:4
                        if !tri_check(K3/2,L/2,Lp/2);continue;end
                        for K1 = 0:2:4
                            for K2=0:2:4
                                if !tri_check(K1/2,K2/2,K3/2); continue;end
                                for X = max(abs(Lp-K1),abs(L-K2)):min(Lp+K1,L+K2)
                                    nkey = get_key6j_sym(Lp,K1,X,K2,L,K3)
                                    d6j_int[nkey] = wigner6j(Float64,Lp//2,K1//2,X//2,K2//2,L//2,K3//2)
                                    hit6j_int += 1
                                end
                            end
                        end
                    end
                end
            end
        end
        @timeit to "const:6j_lj" begin      
            d6j_lj  = Dict{UInt64,Float64}()      
            # {1/2 1/2 s;1/2 S123 s}-type
            j1 = j2 = j3 = 1
            for s12 = 0:2:2
                for s45 = 0:2:2
                    for S = 1:2:3
                        nkey = get_key6j_sym(j1,j2,s12,j3,S,s45)
                        d6j_lj[nkey] = wigner6j(Float64,j1//2,j2//2,s12//2,j3//2,S//2,s45//2)
                        hit6j_lj += 1
                    end
                end
            end
            # {j j' JJ; J3' J3 J123}-type includinig {t t' T;1/2 1/2 1/2}-type
            for j = 0:2:2*N3max
                for jp = 0:2:2*N3max
                    for JJ = abs(j-jp):2:j+jp
                        for dJ123 = 1:2:j3max
                            for j3p = abs(jp-dJ123):2:jp+dJ123
                                for j3 = abs(j-dJ123):2:j+dJ123
                                    if !tri_check(j3p/2,j3/2,JJ/2);continue;end
                                    nkey = get_key6j_sym(j,jp,JJ,j3p,j3,dJ123)
                                    d6j_lj[nkey] = wigner6j(Float64,j//2,jp//2,JJ//2,j3p//2,j3//2,dJ123//2)
                                    hit6j_lj += 1
                                end
                            end
                        end
                    end
                end
            end
        end
        println("hit $hit6j_int $hit6j_lj size ",show_size_inMB("d6j_int",d6j_int),show_size_inMB("d6j_lj",d6j_lj))
    end

    ## prep wigner 9j
    d9j_int= Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}()
    d9j_lsj= Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}()
    @timeit to "9j" begin
        for j4=0:2:lmax_9j
            for j7 = 0:2:lmax_9j
                for j1 = abs(j7-j4):j7+j4
                    if (j1+j4+j7)%2 == 1;continue;end                    
                    @timeit to "const:9j_int" begin
                        for j2 = 0:2:4
                            for j3 = abs(j1-j2):j1+j2
                                if (j1+j2+j3)%2==1; continue; end
                                key1 = get_nkey3(j1,j2,j3)
                                if !haskey(d9j_int,key1); d9j_int[key1] = Dict{Int64,Dict{Int64,Float64}}();end
                                for j5=0:2:2
                                    for j6=abs(j4-j5):min(j4+j5,J12max*2)
                                        if (j4+j5+j6)%2 == 1; continue; end
                                        key2 = get_nkey3(j4,j5,j6)
                                        if !haskey(d9j_int[key1],key2); d9j_int[key1][key2] = Dict{Int64,Float64}();end 
                                        for j8 = abs(j2-j5):j2+j5
                                            if (j2+j5+j8)%2 == 1;continue;end
                                            for j9 = max(abs(j3-j6),abs(j7-j8)):min(j3+j6,j7+j8,J12max*2)
                                                if (j3+j6+j9) %2 == 1; continue;end
                                                key3 = get_nkey3(j7,j8,j9)                        
                                                t9j = wigner9j(j1//2,j2//2,j3//2,j4//2,j5//2,j6//2,j7//2,j8//2,j9//2)
                                                d9j_int[key1][key2][key3] = t9j
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end   
                    @timeit to "const:9j_lsj" begin                 
                        j5 = 1
                        for j2 = 0:2
                            for j3 = abs(j1-j2):j1+j2
                                if (j1+j2+j3)%2==1; continue; end
                                key1 = get_nkey3(j1,j2,j3)                            
                                if !haskey(d9j_lsj,key1); d9j_lsj[key1] = Dict{Int64,Dict{Int64,Float64}}();end
                                for j6=abs(j4-j5):j4+j5
                                    if (j4+j5+j6)%2 == 1; continue; end
                                    key2 = get_nkey3(j4,j5,j6)
                                    if !haskey(d9j_lsj[key1],key2); d9j_lsj[key1][key2] = Dict{Int64,Float64}();end 
                                    for j8 = abs(j2-j5):j2+j5
                                        if (j2+j5+j8)%2 == 1;continue;end
                                        for j9 = max(abs(j3-j6),abs(j7-j8)):min(j3+j6,j7+j8)
                                            if (j3+j6+j9) %2 == 1; continue;end
                                            key3 = get_nkey3(j7,j8,j9)
                                            t9j = wigner9j(j1//2,j2//2,j3//2,j4//2,j5//2,j6//2,j7//2,j8//2,j9//2)
                                            d9j_lsj[key1][key2][key3] = t9j
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end    
        end
        println(show_size_inMB("d9j_int",d9j_int),show_size_inMB("d9j_lsj",d9j_lsj))

        @timeit to "HOB" dictHOB = const_d9j_HOB_3NF(N3max,J12max,j3max,dtri,dcgm0,to)
    end
    return dWS3N(dtri,dcgm0,d6j_int,d6j_lj,d9j_int,d9j_lsj,dictHOB)
end 

function const_d9j_HOB_3NF(N3max,J12max,dJ3max,dtri,dcgm0,to;mode="3NF")
    dictHOB= Dict{Int64,Dict{Int64,Dict{Int64,Float64}}}()    
    if mode == "3NF"
        l12max = J12max+2 
        l3max = dJ3max+1
        Lammax = N3max #dJ3max + 3
        lmin = max(N3max,min(l12max,l3max))
        lmax = max(l12max,l3max)
        @timeit to "indirect 9j" begin
            #{j1 j2 j3; j4 j5 j6}
            d6j_int = Dict{UInt64,Float64}()
            jmax = max(l12max,l3max,Lammax)
            num6j = 0
            maxj36 = 0
            @timeit to "6j" for j1 = 0:2:jmax
                for j3 = 0:2:jmax*2
                    for j2 = abs(j1-j3):2:j1+j3
                        for j4 = j1:2:jmax*2
                            for j6 = abs(j2-j4):2:jmax*2
                                if !tri_check(j6/2,j2/2,j4/2); continue; end
                                for j5 = abs(j1-j6):2:j1+j6    
                                    if !tri_check(j4/2,j5/2,j3/2); continue; end
                                    if !tri_check(j1/2,j5/2,j6/2); continue; end
                                    if !(j1+j4<=j2+j5<=j3+j6);continue;end
                                    maxj36 = max(maxj36,j3+j6)
                                    nkey = get_key6j_sym(j1,j2,j3,j4,j5,j6)
                                    if !haskey(d6j_int,nkey)
                                        d6j_int[nkey] = wigner6j(Float64,j1/2,j2/2,j3/2,j4/2,j5/2,j6/2)
                                        num6j += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end
            println("num6j $num6j maxj36 $maxj36 lmin $lmin lmax $lmax jmax $jmax Lammax $Lammax")
            hitHOB = 0
            @timeit to "Store HOB" begin                
                for Lam = 0:Lammax
                    dictHOB[Lam] = Dict{Int64,Dict{Int64,Float64}}()
                end
                for E = 0:N3max
                    for e_nl = 0:E
                        e_NL = E - e_nl
                        for n = 0:div(e_nl,2)
                            l = e_nl - 2*n
                            for N = 0:div(e_NL,2)
                                L = e_NL - 2*N                           
                                for e_nlp = 0:E
                                    e_NLp = E - e_nlp
                                    for np = 0:div(e_nlp,2)
                                        lp = e_nlp - 2*np
                                        for Np = 0:div(e_NLp,2)
                                            Lp = e_NLp - 2*Np
                                            for Lam = max(abs(l-L),abs(lp-Lp)):min(l+L,lp+Lp)
                                                key1, key2 = get_HOB_nlkey(n, l, N, L, np, lp, Np, Lp)
                                                if !haskey(dictHOB[Lam],key1); dictHOB[Lam][key1] = Dict{Int64,Float64}();end    
                                                tHOB = HObracket_d6j(n, l, N, L, np, lp, Np, Lp, Lam, 1.0/3.0, dtri, dcgm0, d6j_int, to)
                                                dictHOB[Lam][key1][key2] = tHOB
                                                hitHOB += 1
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif mode == "LabTrans"

    else
        @error "mode=$mode is not supported!"
    end

    println("HOB # ",hitHOB, " size ",show_size_inMB("HOB",dictHOB))
    return  dictHOB
end

function show_size_inMB(label,mydict)
    tval = Base.summarysize(mydict)/(10^6)
    if tval < 0.01 
        return "  $label "*@sprintf("%7.2f",tval*10^3)*" KB"
    elseif tval < 10000
        return "  $label "*@sprintf("%7.2f",tval)*" MB"
    else 
        return "  $label "*@sprintf("%7.2f",tval/10^3)*" GB"
    end
end

function zero_9j_check(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    j1 = Int(j1); j2 = Int(j2); j3 = Int(j3)
    j4 = Int(j4); j5 = Int(j5); j6 = Int(j6)
    j7 = Int(j7); j8 = Int(j8); j9 = Int(j9)
    tf = false
    sig = j1 + j2 + j3 + j4 + j5 + j6 + j7 + j8 + j9
    if j1 == j4 && j2 == j5 && j3 == j6 && sig % 2 == 1; tf = true;end
    if j7 == j4 && j8 == j5 && j9 == j6 && sig % 2 == 1; tf = true;end
    if j7 == j1 && j8 == j2 && j9 == j3 && sig % 2 == 1; tf = true;end
    if j1 == j2 && j4 == j5 && j7 == j8 && sig % 2 == 1; tf = true;end
    if j1 == j3 && j4 == j6 && j7 == j9 && sig % 2 == 1; tf = true;end
    if j2 == j3 && j5 == j6 && j8 == j9 && sig % 2 == 1; tf = true;end
    if tf ; return nothing;end

    xmin = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))
    xmax = min(    j1+j9 ,    j2+j6 ,    j4+j8)
    println("{$j1 $j2 $j3\n $j4 $j5 $j6\n $j7 $j8 $j9} sig $sig Set ",
            Set([j1,j2,j3,j4,j5,j6,j7,j8,j9]), "  xmin $xmin xmax $xmax")
    s = 0.0
    for x =xmin:xmax
        t  = (-1)^(2*x) * (2*x+1)
        w1 = wigner6j(Float64,j1,j4,j7,j8,j9,x)
        w2 = wigner6j(Float64,j2,j5,j8,j4,x,j6)
        w3 = wigner6j(Float64,j3,j6,j9,x,j1,j2)
        s += t * w1 * w2 * w3
        println("x $x w1 $w1 w2 $w2 w3 $w3 s $s")
    end
end

function get_HOB_nlkey(n, l, N, L, np, lp, Np, Lp)
    key1 = get_nkey4(n,l,N,L)
    key2 = get_nkey4(np,lp,Np,Lp)
    return key1, key2
end

function get_dictHOB(n12, l12, n3, l3, n45, l45, n6, l6, lambda, dictHOB)
    key1,key2 = get_HOB_nlkey(n12,l12,n3,l3,n45,l45,n6,l6)
    # try
    #     return dictHOB[lambda][key1][key2]
    # catch
    #     println("HOB: Lam $lambda nls $n12 $l12 $n3 $l3, $n45 $l45 $n6 $l6")
    #     println("key1 $key1 key2 $key2")
    # end
    return dictHOB[lambda][key1][key2]
end

function HObracket_naiv(nl, ll, nr, lr, n1, l1, n2, l2, Lam, d::Float64, dtri, dcgm0, to)
    r = 0.0
    ee = 2*nl + ll
    er = 2*nr + lr
    e1 = 2*n1 + l1
    e2 = 2*n2 + l2
    if ee + er != e1 + e2; return r;end
    if !tri_check(ll, lr, Lam);return r;end
    if !tri_check(l1, l2, Lam);return r;end
    phase = (-1.0)^(n1 + n2 + nr + nl) 
    t = sqrt( ( d^(e1 - er)) / ((1.0 + d)^(e1 + e2)))
    m = min(er, e2)    
    for ed = 0:m
        eb = er - ed
        ec = e2 - ed
        ea = e1 - er + ed
        for ld = ed:-2:0
            for lb = eb:-2:0
                if tri_check(ld,lb,lr)==false;continue;end
                for lc = ec:-2:0
                    if !tri_check(ld,lc,l2) ;continue;end
                    for la = ea:-2:0
                        if !tri_check(la,lb,l1);continue;end
                        if !tri_check(la,ll,lc);continue;end
                        t9j = wigner9j(la,lb,l1,lc,ld,l2,ll,lr,Lam)
                        tmp = ((-d)^ed)  * t 
                        tmp *= t9j 
                        tmp *= Ghob(e1, l1, ea, la, eb, lb, dtri, dcgm0)
                        tmp *= Ghob(e2, l2, ec, lc, ed, ld, dtri, dcgm0)
                        tmp *= Ghob(ee, ll, ea, la, ec, lc, dtri, dcgm0)
                        tmp *= Ghob(er, lr, eb, lb, ed, ld, dtri, dcgm0)
                        r += tmp
                    end
                end
            end
        end
    end
    return r * phase
end

"""
associated Legendre Polynomials are calculated with AssociatedLegendrePolynomials.jl
"""
function spherical_harmonics(l, m, θ, φ)
    Ylm = (-1)^m * sqrt( (2*l+1) * factorial(l-m) / ( 4*pi * factorial(l+m)) ) * Plm(l, m, cos(θ)) * (cos(m*φ)+im*sin(m*φ))
    return Ylm
end
