"""
    wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)

calculate Wigner 9j symbols, all j should be given as integer (0,1,...) or halfinteger (3/2, 3//2,...)
"""
function wigner9j(j1,j2,j3,j4,j5,j6,j7,j8,j9)
    s= 0.0
    xmin = max(abs(j1-j9),abs(j2-j6),abs(j4-j8))
    xmax = min(abs(j1+j9),abs(j2+j6),abs(j4+j8))
    @inbounds for x =xmin:xmax
        t  = (-1)^(2*x) * (2*x+1)
        t *= wigner6j(Float64,j1,j4,j7,j8,j9,x)
        t *= wigner6j(Float64,j2,j5,j8,j4,x,j6)
        t *= wigner6j(Float64,j3,j6,j9,x,j1,j2)
        s += t 
    end
    return s
end

"""
    s_wigner9j(j1,j3,j4,j6,j7,j8,j9) 

to calc. wigner9j for specific cases with j2=j5=1/2
"""
function s_wigner9j(j1,j3,j4,j6,j7,j8,j9) 
    s= 0.0
    xmin = max(abs(j1-j9),abs(Int(1//2-j6)),abs(j4-j8))
    xmax = min(abs(j1+j9),abs(Int(1//2+j6)),abs(j4+j8))
    @inbounds for x =xmin:xmax
        t  = (-1.0)^(2*x) * (2.0*x+1.0)
        t *= wigner6j(Float64,j1,j4,j7,j8,j9,x)
        t *= wigner6j(Float64,1//2,1//2,j8,j4,x,j6)
        t *= wigner6j(Float64,j3,j6,j9,x,j1,1//2)
        s += t
    end
    return s
end

"""
    TMtrans(dLECs,xr,wr,xrP,wrP,Rnl,RNL,nTBME,infos,izs_ab,Numpn,V12ab,arr_numst,dict6j,d6j_nabla,X9,U6,to;calc_relcm=false,writesnt=true)

Function to carry out Talmi-Mochinsky transformation
"""
function TMtrans(chiEFTobj,dLECs,xr,wr,xrP,wrP,Rnl,RNL,nTBME,infos,izs_ab,
                 Numpn,V12ab,arr_numst,dict6j,d6j_nabla, X9,U6,to;calc_relcm=false,writesnt=true)
    emax = chiEFTobj.emax
    Nrmax = 2*emax
    nume = zeros(Int64,Nrmax+1,Nrmax+1)
    numknn = [ [ [ zeros(Int64,div(K1,2)+1,div(KK-K1,2)+1) for K1=0:KK] for Lam=0:KK] for KK=0:Nrmax]
    Ndim = [ zeros(Int64,i) for i=1:Nrmax+1]
    Transbk= Float64[]
    
    sp_P5_9j = [ wigner9j(1,1,0,1//2,1//2,0,1//2,1//2,0) wigner9j(1,1,0,1//2,1//2,1,1//2,1//2,1);
                 wigner9j(1,1,2,1//2,1//2,0,1//2,1//2,0) wigner9j(1,1,2,1//2,1//2,1,1//2,1//2,1)]
    cg1s = [clebschgordan(Float64,1,0,1,0,0,0),clebschgordan(Float64,1,0,1,0,2,0)]
    nofst = 0
    nljsnt = [ [0,0] ]; deleteat!(nljsnt,1)
    for temax = 0:emax  #2n +l = temax
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
    
    f_mb,g_mb,w_mb = def_fgw()    
    num=0
    for KK=0:Nrmax
        for Lam=0:KK
            nume[KK+1,Lam+1]=num #nume(KK,Lam)=num
            for K1=0:KK
                K2=KK-K1
                if K2 < K1; continue;end
                for N1=0:div(K1,2) 
                    L1=K1-2*N1
                    for N2=0:div(K2,2)
                        L2=K2-2*N2
                        for K3=0:KK
                            K4=KK-K3
                            if K4 < K3;continue;end
                            for N3=0:div(K3,2)#N3=0,K3/2
                                L3=K3-2*N3
                                for N4=0:div(K4,2)
                                    L4=K4-2*N4
                                    if abs(L1-L2) > Lam || (L1+L2) < Lam;continue;end
                                    if 2*N1+L1+2*N2+L2 != KK;continue;end
                                    if abs(L3-L4) > Lam || (L3+L4) < Lam;continue;end
                                    if 2*N3+L3+2*N4+L4 != KK;continue;end
                                    num=num+1
                                    if KK <= Nrmax                                        
                                        fmosh = gmosh(N1,L1,N2,L2,N3,L3,N4,L4,Lam,1.0,f_mb,g_mb,w_mb)
                                        push!(Transbk,fmosh)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    #println("dim. of HO brackets = $num")
    for KK=0:Nrmax
        tarr = numknn[KK+1]        
        for Lam=0:KK
            num = 0
            tar2 = tarr[Lam+1]
            for K1=0:KK
                K2=KK-K1
                tar3 = tar2[K1+1]
                if K2 < K1; continue;end
                for N1=0:div(K1,2) # do N2=0,K2/2
                    L1=K1-2*N1
                    for N2=0:div(K2,2)
                        L2=K2-2*N2
                        if abs(L1-L2) >Lam || L1+L2 < Lam;continue;end
                        if 2*N1+L1+2*N2+L2 != KK; continue;end
                        num=num+1
                        tar3[N1+1,N2+1]=num
                    end
                end
            end
            Ndim[KK+1][Lam+1]=num
        end
    end
    nljdict = Dict(0=>0); delete!(nljdict,0)
    if writesnt 
        target_nlj = chiEFTobj.target_nlj
        tbme_fmt = chiEFTobj.tbme_fmt
        io=open(chiEFTobj.fn_tbme,"w")
        if tbme_fmt  == "myg"
            println(io, nTBME)
            println(io,"tza,  a, tzb,  b, tzc,  c, tzd,  d, J12,   <ab|v|cd>,   <ab|p_ip_j|cd>/mhw")
        elseif tbme_fmt=="snt"
            if length(target_nlj)!=0
                nljdict = def_sps_snt(emax,chiEFTobj.target_nlj)[2]
            end
            write_spes(chiEFTobj,io,nljsnt,nofst,nTBME,nljdict)
        elseif tbme_fmt=="snt.bin"
            if length(target_nlj)!=0
                nljdict = def_sps_snt(emax,target_nlj)[2]
            end
            write_spes(chiEFTobj,io,nljsnt,nofst,nTBME,nljdict;bin=true)
        end
    end
    tbmes = [ Dict{Vector{Int64},Vector{Vector{Float64}}}() for pnrank=1:3]
    tkeys = [zeros(Int64,4) for i=1:4]
    key6j = zeros(Int64,5)
    @inbounds for ich=1:length(infos)
        @timeit to "misc" begin
            izz,ip,Jtot,ndim=infos[ich]
            pnrank = Int(div(izz,2))+2
            izs = izs_ab[ich]
            vv = zeros(Float64,ndim,ndim)
        end
        @timeit to "vtrans" @inbounds @qthreads for i = 1:ndim
            t2v=zeros(Int64,2)
            t5v=zeros(Int64,5)
            iza,ia,izb,ib = izs[i]
            @inbounds for j = 1:i
                izc,ic,izd,id= izs[j]                
                v12 = vtrans(chiEFTobj,pnrank,izz,ip,Jtot,iza,ia,izb,ib,izc,ic,izd,id,
                             nljsnt,V12ab,X9,U6,nume,numknn,Ndim,Transbk,arr_numst,t2v,t5v,to)
                if chiEFTobj.v_chi_order >= 1
                    vvs = zeros(Float64,5)
                    NLOvs(chiEFTobj,dLECs,vvs,xr,wr,xrP,wrP,Rnl,RNL,cg1s,sp_P5_9j,nljsnt,pnrank,ip,
                          X9,U6,t5v,Jtot,iza,ia,izb,ib,izc,ic,izd,id,to)
                    for i=1:5
                        v12 += vvs[i]
                    end
                end
                vv[i, j] = v12; vv[j, i] = v12                
            end
        end       
        if writesnt 
            write_tbme(chiEFTobj,io,ndim,izs,Jtot,vv,nljsnt,nljdict,tkeys,dict6j,d6j_nabla,key6j;ofst=nofst)
        else
            set_tbme(chiEFTobj,tbmes[pnrank],ndim,izs,Jtot,vv,nljsnt,nljdict,tkeys,dict6j,d6j_nabla,key6j,to;ofst=nofst)
        end
    end
    if writesnt 
        close(io)
        return nothing
    else
        return tbmes
    end
end

function set_tbme(chiEFTobj,tbmes,ndim,izs,Jtot,vv,nljsnt,nljdict,tkeys,dict6j,d6j_nabla,key6j,to;ofst=0)
    target_nlj = chiEFTobj.target_nlj
    @inbounds for i = 1:ndim
        iza,ia,izb,ib = izs[i]
        na,la,ja = nljsnt[ia]
        nb,lb,jb = nljsnt[ib]
        for j = 1:i
            tv = vv[i,j]
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
            vpp = kinetic_tb(tkeys[1],tkeys[2],tkeys[3],tkeys[4],Jtot,dict6j,d6j_nabla,key6j)            
            owtkey!(tkeys[1],fa,fb,fc,fd)
            t = get(tbmes,tkeys[1],false)
            if t == false                           
                tbmes[copy(tkeys[1])] = [ [1.0*Jtot, 0.0, tv, vpp] ]
            else
                push!(tbmes[tkeys[1]], [1.0*Jtot, 0.0, tv, vpp])
            end 
        end
    end
    return nothing
end

function red_nabla_l(n1,l1,n2,l2)
    #""" b*<n1,l1|| nabla ||n2,l2> in l-reduced m.e.
    #N.B.  l1, l2 are not żd """
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

function red_nabla_j(nlj1, nlj2) 
    # b*<j|| nabla ||j2>  N.B.  l1, l2 are not doubled """
    n1, l1, j1 = nlj1
    n2, l2, j2 = nlj2
    ret = (-1)^((3+2*l1+j2)//2) * sqrt(1.0*(j1+1)*(j2+1)) *
        wigner6j(Float64,j1//2, 1, j2//2, l2, 1//2, l1) * red_nabla_l(n1,l1,n2,l2)
    return ret
end
function red_nabla_j(nlj1,nlj2,d6j,key6j) 
    # b*<j|| nabla ||j2>  N.B.  l1, l2 are not doubled """
    # j1//2  j2//2      1;  l2     l1   1//2
    n1, l1, j1 = nlj1
    n2, l2, j2 = nlj2
    if tri_check(j1//2,j2//2,1)==false;return 0.0;end                    
    if tri_check(l1,l2,1)==false;return 0.0;end                    
    if tri_check(j1//2,1//2,l1)==false;return 0.0;end   
    key6j[1] = j1; key6j[2] = j2; key6j[3] = l2; key6j[4] = l1; key6j[5]=0
    t6j = d6j[key6j]
    ret = (-1)^((3+2*l1+j2)//2) * sqrt(1.0*(j1+1)*(j2+1)) *
           t6j * red_nabla_l(n1,l1,n2,l2)
    return ret
end

"""
    kinetic_ob(nlj1, nlj2)

calc. kinetic one-body contribution
"""
function kinetic_ob(nlj1, nlj2)
    #""" <j1 | T/hw | j2> """
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
    kinetic_tb(nljtz1, nljtz2, nljtz3, nljtz4,J, dict6j,d6j_nabla,key6j)

calc. kinetic two-body contribution using preallocated 6j Dict
"""
function kinetic_tb(nljtz1, nljtz2, nljtz3, nljtz4,J, dict6j,d6j_nabla,key6j)
    #""" <j1j2|| -p1*p2/hw ||j3j4>_J """
    n1, l1, j1, tz1 = nljtz1;  n2, l2, j2, tz2 = nljtz2
    n3, l3, j3, tz3 = nljtz3;  n4, l4, j4, tz4 = nljtz4
    nlj1 = @view nljtz1[1:3] 
    nlj2 = @view nljtz2[1:3] 
    nlj3 = @view nljtz3[1:3] 
    nlj4 = @view nljtz4[1:3] 
    norm = 1.0; ret = 1.0
    if nljtz1 == nljtz2; norm *= 1.0/sqrt(2.0);end
    if nljtz3 == nljtz4; norm *= 1.0/sqrt(2.0);end
    key6j[1] = j1; key6j[2] = j2; key6j[3] = j4; key6j[4] = j3; key6j[5] =1
    t = get(dict6j[J+1],key6j,false)
    t6j = 0.0; 
    if t!=false;t6j = dict6j[J+1][key6j]; end     
    key6j[1] = j1; key6j[2] = j2; key6j[3] = j3; key6j[4] = j4; key6j[5] =1
    t = get(dict6j[J+1],key6j,false)
    t6j_2 = 0.0; if t!=false; t6j_2 = dict6j[J+1][key6j]; end 
    
    nabla1324 = red_nabla_j(nlj1,nlj3,d6j_nabla,key6j) * red_nabla_j(nlj2,nlj4,d6j_nabla,key6j)
    if tz1==tz2==tz3==tz4
        nabla1423 = red_nabla_j(nlj1, nlj4,d6j_nabla,key6j) * red_nabla_j(nlj2,nlj3,d6j_nabla,key6j)
        ret = (-1)^((j2+j3)//2+J) * t6j * nabla1324 - (-1)^((j3+j4)//2-J) * (-1)^((j2+j4)//2+J) * t6j_2 * nabla1423
    elseif ( tz1==-1 && tz2 ==1 && tz3 == -1 && tz4 ==1) ||
           ( tz1==1 && tz2 ==-1 && tz3 == 1 && tz4 ==-1)
        ret = (-1)^((j2+j3)//2+J) * t6j * nabla1324
    else
        println("error in kinetic_tb")
        exit()
    end
    return norm*ret
end

function Nfac_jj(na,la,ja,nb,lb,jb)
    s =1.0/ sqrt(1.0+delta(na,nb)*delta(la,lb)*delta(ja,jb))
    return s
end

"""
    HObracket(n1,l1,n2,l2,n,l,N,L,Lam;d=1)

To calc. harmonic oscillator brackets, ``<E_L,e_l:Lam|e_1l_1,e_2l_2:Lam>_d``.
The phase ``(-1)^(N+n+n1+n2)`` is needed to reproduce Tab.1 of [1].

## Reference:
[1] B.Buck& A.C.merchant, Nucl.Phys. A600 (1996) 387-402

[2] G.P.Kamuntavicius et al., Nucl.Phys. A695 (2001) 191-201
"""
function HObracket(n1,l1,n2,l2,n,l,N,L,Lam;d=1)
    e1 = 2*n1+l1
    e2 = 2*n2+l2
    e  = 2*n+l
    E  = 2*N+L
    t1 = d^((e1-e)/2) * (1+d)^(-(e1+e2)/2)
    phase = (-1)^(N+n+n1+n2)
    sum = 0.0
    for ed = 0:min(e2,e)
        ec = e2-ed
        ea = E-ec
        eb = e1-ea
        for la = ea:-2:0
            na = div(ea-la,2)
            for lb = eb:-2:0
                nb = div(eb-lb,2)
                for lc = ec:-2:0
                    nc = div(ec-l,2)
                    for ld = ed:-2:0
                        nd = div(ed-ld,2)
                        t2 = (-d)^ed * wigner9j(la,lb,l1,lc,ld,l2,L,l,Lam)
                        if t2 == 0.0;continue;end
                        t2 *= Ghob(e1,l1,ea,la,eb,lb) * Ghob(e2,l2,ec,lc,ed,ld)
                        if t2 == 0.0;continue;end
                        t2 *= Ghob(E,L,ea,la,ec,lc) *  Ghob(e,l,eb,lb,ed,ld)
                        if t2 == 0.0;continue;end
                        sum += t2
                    end
                end
            end
        end
    end
    return Float64(t1*sum)*phase
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

"""
    CGm0(l1,l2,l)

CG coefficients for special case (l1,m1=0,l2,m2=0|l,m=0)
"""
function CGm0(l1,l2,l)
    if (l1+l2+l) % 2 == 1;return 0.0;end
    g = div(l1+l2+l,2)
    r = (-1)^(g-l)
    r *= sqrt( (2*l+1) * factorial(big(g))/factorial(big(g-l))
    * mydoublefactorial(2*g-2*l1-1) / factorial(big(g-l1))
    * mydoublefactorial(2*g-2*l2-1) / factorial(big(g-l2))
    * mydoublefactorial(2*g-2*l-1) /  doublefactorial(2*g+1) )
    # r *= sqrt( (2*l+1) * factorial(g)/factorial(g-l)
    #            * mydoublefactorial(2*g-2*l1-1) / factorial(g-l1)
    #            * mydoublefactorial(2*g-2*l2-1) / factorial(g-l2)
    #            * mydoublefactorial(2*g-2*l-1) /  doublefactorial(2*g+1) )
    return r
end
function Ghob(e1, l1, ea, la, eb, lb,dtri,dcgm0,keycg;to=nothing) 
    keycg[3] = l1
    cg = 0.0
    if la <= lb
        keycg[1] = la; keycg[2] = lb; cg = dcgm0[keycg]
    else
        keycg[1] = lb; keycg[2] = la; cg = dcgm0[keycg] *( (-1)^(la+lb-l1))
    end
    r  = cg * sqrt( (2.0*la + 1.0) * (2.0*lb + 1.0) ) 
    keycg[1]=e1-l1; keycg[2]=ea-la; keycg[3]=eb-lb; r *= sqrt( dtri[keycg] )
    keycg[1]=e1+l1+1; keycg[2]=ea+la+1; keycg[3]=eb+lb+1; r *= sqrt( dtri[keycg] )
    return r
end 
function Ghob(e1, l1, ea, la, eb, lb;to=nothing)
    cg = CGm0(la,lb,l1)
    r  = cg * sqrt( (2*la + 1) * (2*lb + 1) 
                   * trinomial(e1-l1, ea-la, eb-lb) 
                   * trinomial(e1+l1+1, ea+la+1, eb+lb+1) )
    return r
end 

function trinomial(n1,na,nb)
    return mydoublefactorial(n1) /( mydoublefactorial(na) * mydoublefactorial(nb))
end

function tri_check(ja,jb,jc)
    TF= true
    if ja+jb < jc ; return false;end
    if jb+jc < ja ; return false;end
    if jc+ja < jb ; return false;end
    if abs(ja-jb) > jc ; return false;end
    if abs(jb-jc) > ja ; return false;end
    if abs(jc-ja) > jb ; return false;end
    return true
end

function gmosh2(nl, ll, nr, lr, n1, l1, n2, l2, lm, d::Float64,dWs,to)
    # dWS: dWS2n or dWS3n struct
    # NN => mainly d=1,3N => d=1/3 
    r = 0.0
    ee = 2*nl + ll
    er = 2*nr + lr
    e1 = 2*n1 + l1
    e2 = 2*n2 + l2
    if ee + er != e1 + e2; return r;end
    if !tri_check(ll, lr, lm);return r;end
    if !tri_check(l1, l2, lm);return r;end
    keycg = dWs.keycg[threadid()]
    dcgm0 = dWs.dcgm0
    dtri = dWs.dtri
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
                        t9j = wigner9j(la,lb,l1,lc,ld,l2,ll,lr,lm)
                        #t9j = 1.e-6
                        tmp = ((-d)^ed)  * t                        
                        tmp *= t9j 
                        tmp *= Ghob(e1, l1, ea, la, eb, lb, dtri, dcgm0, keycg)
                        tmp *= Ghob(e2, l2, ec, lc, ed, ld, dtri, dcgm0, keycg) 
                        tmp *= Ghob(ee, ll, ea, la, ec, lc, dtri, dcgm0, keycg)
                        tmp *= Ghob(er, lr, eb, lb, ed, ld, dtri, dcgm0, keycg)
                        r += tmp
                    end
                end
            end
        end
    end
    r = r * phase
    return r
end

# #to prepare Mochinsky brakets
function gmosh(n,l,nc,lc,n1,l1,n2,l2,lr,d,f_mb,g_mb,w_mb)
    ret = 0.0
    if n+n+nc+nc+l+lc-n1-n1-n2-n2-l1-l2 != 0 ;return ret;end
    if l+lc-lr < 0 || l1+l2-lr < 0 ;return ret;end
    if abs(l-lc)-lr > 0 || abs(l1-l2)-lr > 0  ;return ret;end
    dl=log(d)
    d1l=log(d+1.0)

    bb =  f_mb[n1+1]+f_mb[n2+1]+f_mb[n+1]-f_mb[nc+1]
    bb += g_mb[n1+l1+1]+g_mb[n2+l2+1]-g_mb[n+l+1]-g_mb[nc+lc+1]
    ba =  w_mb[l1+1]+w_mb[l2+1]+w_mb[lc+1]+w_mb[l+1]
    ba += f_mb[l1+l2-lr+1]+f_mb[l+lc+lr+2]+f_mb[l+lc-lr+1]+f_mb[lc+lr-l+1]
    ba += f_mb[lr+l-lc+1]-f_mb[l1+l2+lr+2]-f_mb[l1+lr-l2+1]-f_mb[l2+lr-l1+1]-l*d1l

    ip=lr+n+n1+n2
    p=1+2*(div(ip,2)*2-ip)

    anorm=p*exp(0.5*(bb+ba))
    y=0.0
    j1f=l+1

    for j1=1:j1f
        j2=l+2-j1
        k1i=abs(l1-j1+1)+1
        k1f=l1+j1
        for k1=k1i:2:k1f
            m1f=n1-div(j1+k1-l1,2)+2
            if m1f-1 < 0 ;continue;end
            k2i=max(abs(l2-j2+1),abs(lc-k1+1))+1
            k2f=min(l2+j2,lc+k1)
            if k2i-k2f > 0 ;continue;end
            for k2=k2i:2:k2f
                m2f=n2-div(j2+k2-l2,2)+2
                if m2f-1 < 0 ;continue;end
                ip=j2-1+div(l1+k1+j1+l2+k2+j2,2)
                p=1+2*(div(ip,2)*2-ip)
                bc = 0.5*((k1+j2-2)*dl-(k1+k2-2)*d1l)
                bc +=  f_mb[k1+l1-j1+1]+f_mb[k1+k2-lc-1]+f_mb[k2+l2-j2+1]-f_mb[k1+l1+j1]-f_mb[k1+k2+lc]
                bc += -f_mb[k2+l2+j2]+w_mb[k1]+w_mb[k2]+f_mb[div(k1+l1+j1,2)]
                bc +=  f_mb[div(k1+k2+lc,2)]+f_mb[div(k2+l2+j2,2)]              
                bc += -f_mb[div(k1+l1-j1,2)+1]-f_mb[div(l1+j1-k1,2)+1]
                bc += -f_mb[div(j1+k1-l1,2)]-f_mb[div(k1+k2-lc,2)]
                bc += -f_mb[div(k2+lc-k1,2)+1]-f_mb[div(lc+k1-k2,2)+1] 
                bc += -f_mb[div(k2+l2-j2,2)+1]-f_mb[div(l2+j2-k2,2)+1]
                bc += -f_mb[div(j2+k2-l2,2)]
                #println("bc $bc")

                cfac=p*exp(bc)
                sxy=0.0
                ixf=min(k1+k1,k1+k2-lc)-1
                for ix=1:ixf
                    iyi=max(1,ix+j1+l2-k1-lr)
                    iyf=min(l2+l2+1,l1+l2-lr+1,l2+lc+ix-k1-j2+2)
                    if iyi-iyf > 0;continue;end
                    for iy=iyi:iyf
                        ip=ix+iy
                        p=1+2*(div(ip,2)*2-ip)
                        bxy =  f_mb[k1+k1-ix]+f_mb[l2+l2-iy+2]
                        bxy += f_mb[k2+lc-k1+ix]+f_mb[l1+lr-l2+iy] 
                        bxy += -f_mb[ix]-f_mb[iy]-f_mb[k1+k2-lc-ix]
                        bxy += -f_mb[l1+l2-lr-iy+2]-f_mb[k1-l2+lr-j1+iy-ix+1]
                        bxy += -f_mb[l2-k1+lc-j2+ix-iy+3]
                        sxy += p*exp(bxy)
                    end
                end
                s=cfac*sxy
                sm=0.0
                for m1=1:m1f
                    m2i=max(1,nc-m1-div(k1+k2-lc,2)+3)
                    if m2i-m2f > 0;continue;end
                    for m2=m2i:m2f
                        ip=m1+m2
                        p=1+2*(div(ip,2)*2-ip)
                        bm =  (m1-1)*dl-(m1+m2-2)*d1l+g_mb[1]
                        bm += g_mb[m1+m2+div(k1+k2+lc,2)-2]-g_mb[k1+m1-1]
                        bm += -g_mb[k2+m2-1]+f_mb[m1+m2+div(k1+k2-lc,2)-2]
                        bm += -f_mb[m1]-f_mb[m2]-f_mb[n1-m1-div(j1+k1-l1,2)+3]
                        bm += -f_mb[n2-m2-div(j2+k2-l2,2)+3]-f_mb[m1+m2-nc+div(k1+k2-lc,2)-2]
                        sm=sm+p*exp(bm)
                    end
                end
                y += s*sm
            end
        end
    end
    ret = anorm*y
    return ret
end

function def_fgw(;maxjj=200)
    f_mb = zeros(Float64,maxjj)
    g_mb = zeros(Float64,maxjj)
    w_mb = zeros(Float64,maxjj)
    f_mb[1]=0.0
    g_mb[1]=log(0.5)
    w_mb[1]=0.0
    for i=2:maxjj
       a=i-1
       f_mb[i]=f_mb[i-1]+log(a)
       g_mb[i]=g_mb[i-1]+log(a+0.5)
       w_mb[i]=log(a+a+1.0)
    end
    return f_mb,g_mb,w_mb
end

function vtrans(chiEFTobj,pnrank,izz,ip,Jtot,iza,ia,izb,ib,izc,ic,izd,id,
                nljsnt,V12ab,X9,U6,nume,numknn,Ndim,Transbk,arr_numst,t2v,t5v,to)
    emax = chiEFTobj.emax
    ret = 0.0
    na,la,jda = nljsnt[ia]; nb,lb,jdb = nljsnt[ib]
    nc,lc,jdc = nljsnt[ic]; nd,ld,jdd = nljsnt[id]
    lrmax = jmax + 1
    Nrmax = 2*emax   
    mab=2*na+la+2*nb+lb; mcd=2*nc+lc+2*nd+ld
    TF = false
    if izz != iza+izb || izz != izc+izd; TF;end
    if Jtot > div(jda+jdb,2) || Jtot < abs(div(jda-jdb,2)); TF=true;end
    if Jtot > div(jdc+jdd,2) || Jtot < abs(div(jdc-jdd,2)); TF=true;end
    if (-1)^(la+lb) != ip || (-1)^(lc+ld) != ip; TF=true;end
    if (izz==2 || izz==-2) && ia==ib && (-1)^(Jtot)==-1; TF=true;end
    if (izz==2 || izz==-2) && ic==id && (-1)^(Jtot)==-1; TF=true;end
    if TF; return ret;end   
    U6_j = U6[Jtot+1]    
    for S=0:1
        tX9 = X9[S+1]
        U6_s = U6_j[S+1]
        tarr_numst = arr_numst[pnrank][S+1]
        lmax1=min(Jtot+S,la+lb)
        lmin1=max(abs(Jtot-S),abs(la-lb))
        if lmin1 > lmax1;continue;end
        lmax2=min(Jtot+S,lc+ld)
        lmin2=max(abs(Jtot-S),abs(lc-ld))
        if lmin2 > lmax2;continue;end
        @inbounds for Lam=lmin1:lmax1
            Ja = jda-2*la; Jb = jdb-2*lb
            t5v[1] = la;t5v[2] = Ja;t5v[3] = lb;
            t5v[4] = Jb;t5v[5] = Lam
            ttX9 = tX9[Jtot+1]            
            x1= get(ttX9,t5v,0.0) * (-1)^Lam
            U6_lam1 = U6_s[Lam+1]
            @inbounds for Lamp=lmin2:lmax2
                Jc = jdc-2*lc; Jd = jdd-2*ld
                t5v[1] = lc;t5v[2] = Jc;t5v[3] = ld
                t5v[4] = Jd;t5v[5] = Lamp
                x2=get(ttX9,t5v,0.0) * (-1)^Lamp
                kncu=div(mab,2)
                U6_lam2 = U6_s[Lamp+1]                
                @inbounds for Ncm=0:kncu
                    klcmax=min((mab-2*Ncm),(mcd-2*Ncm))
                    if klcmax < 0;continue;end
                    @inbounds for Lcm=0:klcmax
                        klrmi1=abs(Lcm-Lam)
                        klrma1=min(lrmax,Lcm+Lam)
                        U6_Lcm1 = U6_lam1[Lcm+1]
                        U6_Lcm2 = U6_lam2[Lcm+1]
                        @inbounds for lr1=klrmi1:klrma1
                            nx1=mab-2*Ncm-(lr1+Lcm)
                            nr1=div(nx1,2)
                            if nr1 < 0 || nx1!=2*nr1;continue;end
                            y1 = 0.0
                            if 2*nr1+lr1+2*Ncm+Lcm <= Nrmax && mab <= Nrmax
                                KK = 2*nr1+lr1 + 2*Ncm+Lcm 
                                s_numknn = numknn[KK+1][Lam+1]
                                ndim = Ndim[KK+1][Lam+1]
                                y1=trbknum(nr1,lr1,Ncm,Lcm,
                                           na,la,nb,lb,Lam,
                                           nume,s_numknn,ndim,Transbk,t2v)
                            end                            
                            if abs(y1) < 1.e-10;continue;end
                            U6_lr1 = U6_Lcm1[lr1+1]
                            klrmi2=abs(Lcm-Lamp); klrma2=Lcm+Lamp
                            @inbounds for lr2=klrmi2:klrma2
                                if lr1%2 != lr2%2;continue;end
                                if lr2 > lrmax;continue;end
                                nx2=mcd-2*Ncm-(lr2+Lcm); nr2=div(nx2,2)
                                if nr2 < 0 || (nx2!=2*nr2);continue;end
                                y2 = 0.0
                                if 2*nr2+lr2+2*Ncm+Lcm <= Nrmax && mcd <= Nrmax
                                    KK = 2*nr2+lr2 + 2*Ncm+Lcm 
                                    s_numknn = numknn[KK+1][Lamp+1]
                                    ndim = Ndim[KK+1][Lamp+1]
                                    y2=trbknum(nr2,lr2,Ncm,Lcm,
                                               nc,lc,nd,ld,Lamp,
                                               nume,s_numknn,ndim,Transbk,t2v)
                                end
                                if abs(y2) < 1.e-10;continue;end
                                mj1=abs(lr1-S); mj2=abs(lr2-S); mj3=abs(Jtot-Lcm)
                                kjmin=max(mj1,mj2,mj3)
                                kjmax=min(lr1+S,lr2+S,Jtot+Lcm,6)
                                if kjmin > kjmax;continue;end
                                U6_lr2 = U6_Lcm2[lr2+1]
                                sumv=0.0
                                @inbounds for Jrel=kjmin:kjmax
                                    zu1=U6_lr1[Jrel-abs(lr1-S)+1]
                                    if abs(zu1) < 1.e-10;continue;end
                                    zu2=U6_lr2[Jrel-abs(lr2-S)+1]
                                    if abs(zu2) < 1.e-10;continue;end
                                    izfac=0
                                    if izz==-2 || izz==2
                                        izfac=1+(-1)^(lr1+S)
                                    end
                                    if izz==0;izfac=1;end
                                    rv12=0.0
                                    if izfac!=0
                                        num= tarr_numst[Jrel+1][lr1-abs(Jrel-S)+1][lr2-abs(Jrel-S)+1]
                                        rv12=V12ab[num][nr1+1,nr2+1]
                                    end
                                    sumv += zu1*zu2*izfac*rv12
                                end
                                zxy=x1*x2*y1*y2
                                ret += sumv*zxy
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
    end
    return ret
end

"""
    prepareX9U6(Nrmax;to=nothing)

return 6j/9j dict:
for 6j => `jj`->`S`->`lam`->`lc`->`lr`, for 9j => `S`->`J`->`key`= [`la`,`nja`,`lb`,`njb`,`lam`]
"""
function prepareX9U6(Nrmax;to=nothing)
    jrange = max(Nrmax+1,2*jmax+2)
    X9 = [ [ Dict( [0,0] => 0.0) for J=0:jrange ] for S=0:1]
    for S=0:1;for J=0:jrange; delete!(X9[S+1][J+1], [0,0]);end;end
    jrmax = jmax
    lrmax = jrmax+1
    hit6 = 0; hit9=0
    U6 = [[[[[ zeros(Float64,lr+iss-abs(lr-iss)+1) for lr=0:lrmax ] for lc=0:Nrmax] for lam=0:Nrmax] for iss=0:1] for jj=0:Nrmax+1]
    for lc =0:Nrmax
        for lr =0:lrmax
            for lam=0:Nrmax
                if lam < abs(lc-lr) || lam > lc+lr; continue;end
                for iss=0:1
                    for jj=abs(lam-iss):lam+iss
                        for jr=abs(lr-iss):lr+iss
                            if jr > jrmax;continue;end
                            if jr < abs(lc-jj) || jr > lc+jj;continue;end
                            sfac = sqrt((2.0*lam+1.0)*(2.0*jr+1.0))*(-1.0)^(lc+lr+iss+jj)
                            tmp =  sfac * wigner6j(Float64,lc,lr,lam,iss,jj,jr)
                            U6[jj+1][iss+1][lam+1][lc+1][lr+1][jr-abs(lr-iss)+1] = tmp
                            hit6 += 1
                        end
                    end
                end
            end
        end
    end
    for la=0:Nrmax
        for lb=0:Nrmax
            for nja=-1:2:1
                jda = 2*la + nja 
                if jda < 0;continue;end
                for njb=-1:2:1
                    jdb = 2*lb + njb 
                    if jdb < 0;continue;end
                    for lam=abs(la-lb):la+lb
                        if lam > Nrmax;continue;end
                        for iss=0:1
                            for jj=abs(lam-iss):lam+iss
                                sfac = sqrt((jda+1.0)*(jdb+1.0)*(2.0*lam+1.0)*(2.0*iss+1.0))
                                X9[iss+1][jj+1][
                                     [la,nja,lb,njb,lam]
                                 ] = sfac .* s_wigner9j(la,jda//2,
                                                        lb,jdb//2,
                                                        lam,iss,jj)
                                hit9 += 1
                            end
                        end
                    end
                end
            end
        end
    end
    return X9,U6
end   


function overwritekeyHOB!(key,N,Lam,n,lam,n1,l1,n2,l2,L)
    key[1] =N; key[2] =Lam; key[3] =n;
    key[4] =lam; key[5] =n1; key[6] =l1;
    key[7] =n2; key[8] =l2; key[9] =L
    phase = 1.0
    if  (2*N+Lam > 2*n+lam) && (2*n1+l1 > 2*n2+l2)
        key[1]=n; key[2]=lam; key[3]=N; key[4]=Lam
        key[5]=n2; key[6]=l2; key[7]=n1; key[8]=l1
        phase = (-1)^(Lam+l2)
    end
    return phase
end 

"""
    PreCalcHOB(chiEFTobj,to)

see also "struct HarmonicOscillatorBrackets" in hartreefock.jl/def_struct.jl

In the early version, dict9j is defined as
dict9j = [ [ Dict{Vector{Int64},Float64}() for S=0:1 ] for J=0:Jmax] with key~[la,ja,lb,jb,L]
Dict using array as key is slow, so this was fixed to...
dict9j => nested array J=>S=>L=>ja=>la=>jb=>lb
Note that div(j,2)+1 will be used as idx for ja&jb.

The same can be said for HOBs
HOBs => nested array N=>n=>Lam=>lam=>L=>na=>nb=>la (lb is automatically determined)
"""
function PreCalcHOB(emax,chiEFTobj,to;io=stdout)
    Nnmax = chiEFTobj.Nnmax
    Nmax = max(2*emax,Nnmax)
    if emax >= 10; Nmax = Nnmax + 10;end
    Jmax = jmax2 = 2*emax + 1
    Lmax = Jmax + 1
    lmax = emax + 1
    e2max = emax * 2

    ### trinomial
    dtri = Dict{Vector{Int64},Float64}()
    for l = 0:2*Nmax
        for l1 = 0:2*Nmax
            for l2 = 0:2*Nmax
                key = [l1,l2,l]
                dtri[key] = trinomial(l1,l2,l)
            end
        end
    end
    ### CG coeff for special case
    Nmax2 = Nmax*2
    dcgm0 = Dict{Vector{Int64},Float64}()
    hitCG = 0
    for l = 0:Nmax2
        for l1 = 0:Nmax2
            for l2 = abs(l-l1):l+l1
                if !tri_check(l1,l2,l);continue;end
                if l1 > l2;continue;end
                key = zeros(Int64,3)
                key[1] = l1; key[2]=l2; key[3]=l
                dcgm0[key] = CGm0(l1,l2,l)
                hitCG +=1
            end
        end
    end
    keycg = [ zeros(Float64,3) for i=1:nthreads()]        
    dWS = dWS2n(dtri,dcgm0,keycg)

    #### 9j with {la 1/2 ja; lb 1/2 jb; L S J} structure
    num9j = 0
    dict9j = [ [ [ [ [ [ [ 0.0 for L=0:Lmax] for lb=0:lmax] for jb=1:2:jmax2] for la=0:lmax] for ja=1:2:jmax2] for S=0:1] for J=0:Jmax]
    for J = 0:Jmax
        tJ = dict9j[J+1]
        for S = 0:1
            tS = tJ[S+1]
            for ja = 1:2:jmax2
                tja = tS[div(ja,2)+1]
                for la = div(ja-1,2):div(ja+1,2)
                    if !tri_check(2*la,1,ja);continue;end
                    tla = tja[la+1]
                    for jb = 1:2:jmax2
                        if tri_check(ja,jb,J*2)==false;continue;end
                        tjb = tla[div(jb,2)+1]
                        for lb = div(jb-1,2):div(jb+1,2)
                            if !tri_check(2*lb,1,jb);continue;end
                            tlb = tjb[lb+1]
                            for L = abs(la-lb):la+lb
                                if !tri_check(L,S,J);continue;end
                                t9j = wigner9j(la,1//2,ja//2,lb,1//2,jb//2,L,S,J)
                                tlb[L+1] = t9j
                                num9j +=1
                            end
                        end
                    end                
                end
            end
        end
    end

    ### Calc. HObrackets
    # To reduce total number of HOBs stored, 2n_i+l_i > 2n_j+l_j case is not considered
    #HOBs = [[[[[[[[ 0.0 for l1 = 0:2*N+Lam+2*n+lam-2*n1-2*n2] for n2 =0:div(2*N+Lam+2*n+lam,2)-n1] for n1=0:div(2*N+Lam+2*n+lam,2)] for L = 0:e2max-2*N-2*n] for lam = 0:e2max-2*N-2*n-Lam] for Lam =0:e2max-2*N-2*n] for n=0:e2max-N] for N=0:e2max]
    # used order n_ab,N_ab,lam_ab,Lam_ab,Lab+1,nb+1,na+1,lb+1 or N_ab,n_ab,Lam_ab,lam_ab,Lab+1,na+1,nb+1,la+1
    #new (faster, but more memory greedy)
    HOBs = Dict{Int64, Dict{Int64,Float64}}()
    HOBkeys = Vector{Int64}[ ]
    hit = 0 
    for N=0:e2max
        for n = 0:e2max-N
            Lam_max = e2max-2*N-2*n
            for Lam = 0:Lam_max
                lam_max = Lam_max-Lam
                for lam = 0:lam_max
                    e2 = 2*N+Lam + 2*n+lam 
                    nkey1 = get_nkey_from_key6j(N,n,Lam,lam,0) 
                    defined = get(HOBs,nkey1,false)
                    if defined == false
                        HOBs[nkey1] = Dict{Int64,Float64}()
                    end
                    for L = abs(Lam-lam):Lam+lam
                        for n1=0:div(e2,2)
                            for n2 = 0:div(e2,2)-n1
                                l1max = e2-2*n1-2*n2
                                for l1 = 0:l1max
                                    l2 = e2-2*n1-2*n2-l1
                                    e_1 = 2*n1 + l1; e_2 = 2*n2 + l2
                                    if (e_1 > e_2) && (2*N+Lam > 2*n+lam); continue;end
                                    if (l1+l2+lam+Lam)%2 > 0;continue;end
                                    if !tri_check(l1,l2,L);continue;end
                                    nkey2 = get_nkey_from_key6j(L,n1,n2,l1,0)
                                    push!(HOBkeys,[nkey1,nkey2])
                                    hit += 1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    tkeys = [ zeros(Int64,4) for i=1:nthreads()]        
    @inbounds @threads for i = 1:hit
        nkey1,nkey2 = HOBkeys[i]
        tkey = tkeys[threadid()]
        get_abcdarr_from_intkey!(nkey1,tkey)
        N = tkey[1]; n = tkey[2]; Lam = tkey[3]; lam = tkey[4]
        get_abcdarr_from_intkey!(nkey2,tkey)
        L = tkey[1]; n1 = tkey[2]; n2 = tkey[3]; l1= tkey[4]
        e2 = 2*N+Lam+2*n+lam; l2 = e2-2*n1-2*n2-l1
        tHOB = gmosh2(N,Lam,n,lam,n1,l1,n2,l2,L,1.0,dWS,to)
        HOBs[nkey1][nkey2] = tHOB
    end  
    println(io,"@emax $emax ","hitCG $hitCG dWS  =>", @sprintf("%7.2f",Base.summarysize(dWS)/1024/1024)," MB ",
            "  9j($num9j) =>", @sprintf("%7.2f",Base.summarysize(dict9j)/1024/1024)," MB ",
            "  HOB ($hit)=>",@sprintf("%7.2f",Base.summarysize(HOBs)/1024/1024), " MB")
    dWS = nothing; dcgm0 =nothing; dtri=nothing 
    return dict9j,HOBs
end

function trbknum(Nr,Lr,Nc,Lc,Na,La,Nb,Lb,Lam,nume,s_numknn,ndim,Transbk,t2v)
    ret = 0.0    
    L1=0;L2=0; L3=0;L4=0; N1=0; N2=0;N3=0; N4=0
    phase=0.0
    Kr=2*Nr+Lr
    Kc=2*Nc+Lc
    Ka=2*Na+La
    Kb=2*Nb+Lb
    if (Kr+Kc != Ka+Kb) || abs(Lr-Lc) > Lam || Lr+Lc < Lam || abs(La-Lb) > Lam || La+Lb < Lam
        return ret
    end
    if Kr <= Kc && Ka <= Kb
        N1=Nr; L1=Lr; N2=Nc;L2=Lc; N3=Na; L3=La; N4=Nb; L4=Lb; phase=1.0
    elseif Kr > Kc && Ka <= Kb
        N1=Nc; L1=Lc; N2=Nr;L2=Lr; N3=Na; L3=La; N4=Nb; L4=Lb; phase=(-1.0)^(La+Lam)
    elseif Kr <= Kc && Ka > Kb
        N1=Nr; L1=Lr; N2=Nc;L2=Lc; N3=Nb; L3=Lb; N4=Na; L4=La; phase=(-1.0)^(Lc+Lam)
    elseif Kr > Kc && Ka > Kb
        N1=Nc; L1=Lc; N2=Nr;L2=Lr; N3=Nb; L3=Lb; N4=Na; L4=La; phase=(-1.0)^(Lr+La)
    end
    K1=2*N1+L1; K3=2*N3+L3; KK=Kr+Kc
    numknn_1 = s_numknn[K1+1][N1+1,N2+1] 
    numknn_2 = s_numknn[K3+1][N3+1,N4+1]
    num=nume[KK+1,Lam+1]+ndim*(numknn_1-1) +numknn_2
    ret = Transbk[num]*phase
    return ret
end 

const l2l = [ wigner3j(Float64,l,2,l,0,0,0) for l=0:8]
const l2lnd =[[ wigner3j(Float64,l1,2,l2,0,0,0) for l2=0:8] for l1=0:8]

"""
    prep_wsyms()

cg1s = (1,0,l,0|l',0), cg2s = (2,0,l,0|l',0)
"""
function prep_wsyms()
    lmax = 7
    dim = lmax+3
    cg1s = zeros(Float64,dim,dim)
    cg2s = zeros(Float64,dim,dim)
    d6_121 = zeros(Float64,dim,dim,dim)
    d6_222 = zeros(Float64,dim,dim,dim)
    d6_21 = zeros(Float64,dim,dim,dim,dim)
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
    return cg1s,cg2s,d6_121,d6_21,d6_222,d9_12
end

""" 
    readsnt(sntf,Anum;eachA=false,pnfac=1.0) 

to read sntfile. This is slightly different from readsnt() in ShellModel.jl
"""
function readsnt(sntf,Anum;eachA=false,pnfac=1.0) 
    f = open(sntf,"r");tlines = readlines(f);close(f)
    lines = rm_comment(tlines)
    line = lines[1]
    lp,ln,cp,cn = map(x->parse(Int,x),rm_nan(split(line," ")))
    p_sps = [[0,0]];deleteat!(p_sps,1)
    n_sps = [[0,0]];deleteat!(n_sps,1)
    dictsps = Dict([0,0,0,0]=>0);delete!(dictsps,[0,0,0,0])
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
    dictTBMEs=[ Dict([0,0,0,0,0]=>0.0) for pnrank=1:3]
    for i=1:3
        delete!(dictTBMEs[i],[0,0,0,0,0])
    end
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

function jj_std(sps,dictsps,dictTBMEs;fname="")
    if fname=="";fname = "monopole"*fname*".dat"
    else;fname = "monopole_"*fname*".dat";end
    ln = length(sps)
    nmax = maximum([sps[i][1] for i=1:ln])
    nmin = minimum([sps[i][1] for i=1:ln])
    lmax = maximum([sps[i][2] for i=1:ln])
    lmin = minimum([sps[i][2] for i=1:ln])
    jmax = maximum([sps[i][3] for i=1:ln])
    jmin = minimum([sps[i][3] for i=1:ln])
    emax = maximum( [ 2*sps[i][1]+sps[i][2] for i=1:length(sps)])
    emin = minimum( [ 2*sps[i][1]+sps[i][2] for i=1:length(sps)])
    Jmax = div(2*jmax+1,2)
    Lmin = 0; Lmax = 2*lmax
    monodict = Dict(["",""]=> [[0.0,[0.0]]]);delete!(monodict,["",""])
    cpnrank=["pp","pn","nn"]

    keya = [0,0,0,0]; keyb = [0,0,0,0];keyc=[0,0,0,0]; keyd=[0,0,0,0]   # nljtz
    key_abcdJ = [0,0,0,0,0]
    tvjj = zeros(Float64,5)
    for pnrank =1:3
        tdict = dictTBMEs[pnrank]
        for tkey in keys(tdict)
            a,b,c,d,J = tkey
            na,la,ja2,iza = sps[a]; nb,lb,jb2,izb = sps[b]
            nc,lc,jc2,izc = sps[c]; nd,ld,jd2,izd = sps[d]
            keya[1]=na; keya[2]=la; keya[3]=ja2; keya[4]=iza
            keyb[1]=nb; keyb[2]=lb; keyb[3]=jb2; keyb[4]=izb
            keyc[1]=nc; keyc[2]=lc; keyc[3]=jc2; keyc[4]=izc
            keyd[1]=nd; keyd[2]=ld; keyd[3]=jd2; keyd[4]=izd
            v_target = tdict[tkey]            
            tvjj .= 0.0
            hats = (-1)^J * sqrt((ja2+1.0)*(jb2+1.0)*(jc2+1.0)*(jd2+1.0))
            for L = abs(la-lb):la+lb
                for S=0:1
                    for Lp = abs(lc-ld):lc+ld
                        for Sp = 0:1
                            Jfac  = wigner9j(la,1//2,ja2//2,lb,1//2,jb2//2,L,S,J)
                            Jfac *= wigner9j(lc,1//2,jc2//2,ld,1//2,jd2//2,Lp,Sp,J)
                            for jpa2 = 2*la-1:2:2*la+1
                                if jpa2 <=0;continue;end
                                for jpb2 = 2*lb-1:2:2*lb+1
                                    if jpb2 <=0;continue;end
                                    for jpc2 = 2*lc-1:2:2*lc+1
                                        if jpc2 <=0;continue;end
                                        for jpd2 = 2*ld-1:2:2*ld+1
                                            if jpd2 <=0;continue;end
                                            hatps = sqrt((jpa2+1.0)*(jpb2+1.0)*(jpc2+1.0)*(jpd2+1.0))
                                            LLSS = (2*L+1)*(2*Lp+1)*(2*S+1)*(2*Sp+1)
                                            for Jp = max(abs(Lp-Sp),abs(L-S)):min(Lp+Sp,L+S)
                                                d9j_ab = wigner9j(la,1//2,jpa2//2,lb,1//2,jpb2//2,L,S,Jp)
                                                d9j_cd = wigner9j(lc,1//2,jpc2//2,ld,1//2,jpd2//2,Lp,Sp,Jp)
                                                keya[3] = jpa2; keyb[3]=jpb2;keyc[3] = jpc2; keyd[3]=jpd2                                                
                                                ap = dictsps[keya];bp = dictsps[keyb];cp = dictsps[keyc];dp = dictsps[keyd]
                                                key_abcdJ[1] = ap;key_abcdJ[2] = bp
                                                key_abcdJ[3] = cp;key_abcdJ[4] = dp;key_abcdJ[5] = Jp
                                                ta=ap; tb=bp;tc=cp;td=dp
                                                phase = 1.0
                                                if ap > bp
                                                    ta = bp; tb = ap
                                                    phase *= (-1)^(div(sps[ap][3]+sps[bp][3],2)+Jp+1)
                                                end
                                                if cp > dp
                                                    tc = dp; td = cp
                                                    phase *= (-1)^(div(sps[cp][3]+sps[dp][3],2)+Jp+1)
                                                end
                                                fa=ta; fb=tb;fc=tc;fd=td
                                                if fa > fc || (fa==fc && fb > fd)
                                                    fa=tc;fb=td;fc=ta;fd=tb
                                                end
                                                key_abcdJ[1] = fa;key_abcdJ[2] = fb
                                                key_abcdJ[3] = fc;key_abcdJ[4] = fd;key_abcdJ[5] = Jp
                                                vjj = get(tdict,key_abcdJ,false)                                                
                                                if vjj==false;continue; end
                                                #println("abcd' $ap $bp $cp $dp J' $Jp $nvjj $vjj")
                                                Jpfac = (-1)^Jp * (2*Jp+1) * d9j_ab * d9j_cd
                                                for k = 0:2 
                                                    tidx = k+1; if S==Sp && k==1;tidx=4;end                                                
                                                    kfac = (2*k+1) * wigner6j(Float64,L,S,J,Sp,Lp,k) * wigner6j(Float64,L,S,Jp,Sp,Lp,k)
                                                    tvjj[tidx] += hats*hatps*kfac*Jfac*Jpfac * LLSS * vjj *phase
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
            tvjj[5] = sum( @views tvjj[1:4])
            if abs(v_target-tvjj[5]) > 1.e-9
                println("Error TBME(in jj) mismatch!: $tkey  v ",
                        @sprintf("%12.4e",v_target),"; sum ", @sprintf("%12.4e",tvjj[5]))
            end
            if a==c && b==d && ( a <= b )
                c = cpnrank[pnrank]
                ca = string(na)*chara_l[la+1]*string(ja2)
                cb = string(nb)*chara_l[lb+1]*string(jb2)
                if get(monodict,[c,ca,cb],false) == false
                    monodict[[c,ca,cb]] = [ [J*1.0,copy(tvjj)] ]
                else
                    push!(monodict[[c,ca,cb]],[J*1.0,copy(tvjj)])
                end
            end
        end
        monopole(monodict,fname)
    end
    return nothing
end


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

function sample_std()
    sps,dictsps,dictTBMEs = readsnt("usdb.snt",18)
    spin_tensor_decomposition(sps,dictsps,dictTBMEs)
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

"""
    PreCalc6j(emax)

preallocate 6j used for HF
```math
\\begin{Bmatrix} 
j_a & j_b & J \\\\
j_c & j_d & J'
\\end{Bmatrix}
```
``ja,jb,jc,jd`` are half-integer extended for kinetic_tb
```math
\\begin{Bmatrix} 
j_1/2&  j_2/2&     1 \\\\
  l_2&    l_1&   1/2
\\end{Bmatrix}
```
are needed to get `dict6j[J][key]` with `key = [ja,jb,jd,jc,Jp]`.
Using array as key is in general slow, so the key is transformed into integer in the IMSRG part (But, this sometimes reduces readability...)
"""
function PreCalc6j(emax)
    Jmax = 2*emax+1 # = jmax *2    
    d6j = [ Dict{Vector{Int64},Float64}() for i=0:Jmax]
    for totJ = 0:Jmax
        tdict = d6j[totJ+1]        
        for ja = 1:2:Jmax
            for jb = 1:2:Jmax
                if tri_check(ja,jb,totJ*2)==false;continue;end
                for jd =1:2:Jmax
                    for Jp = 0:Jmax
                        if tri_check(jb,jd,Jp*2)==false;continue;end
                        for jc = 1:2:Jmax
                            if tri_check(jd,totJ*2,jc)==false;continue;end
                            if tri_check(ja,Jp*2,jc)==false;continue;end
                            tdict[[ja,jb,jd,jc,Jp]] = wigner6j(Float64,ja//2,jb//2,totJ,jd//2,jc//2,Jp)
                        end
                    end
                end                
            end
        end
    end
    d6j_nabla = Dict{Vector{Int64},Float64}() 
    Jp2 = 1; totJ=1
    for ja = 1:2:Jmax
        for jb = 1:2:Jmax
            if tri_check(ja,jb,2)==false;continue;end
            for l2 =0:emax
                if tri_check(jb//2,l2,1//2)==false;continue;end
                for l1 = 0:emax
                    if tri_check(l2,totJ,l1)==false;continue;end                    
                    if tri_check(ja//2,1//2,l1)==false;continue;end
                    d6j_nabla[[ja,jb,l2,l1,0]] = wigner6j(Float64,ja//2,jb//2,1,l2,l1,1//2)
                end
            end    
        end
    end
    return d6j,d6j_nabla
end
