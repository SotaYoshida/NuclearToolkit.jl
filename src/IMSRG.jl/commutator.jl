"""
    OpCommutator!(X::Op,Y::Op,ret::Op,HFobj,Chan1b,Chan2bD,dictMono,dict6j,PandyaObj,to) where{Op <: Operator}
    
overwrite `ret` operator by the commutator ``[X,Y]`` 
"""
function OpCommutator!(X::Op,Y::Op,ret::Op,HFobj::HamiltonianNormalOrdered,Chan1b::chan1b,Chan2bD::Chan2bD,dictMono,dict6j,PandyaObj::PandyaObject,to) where{Op <: Operator}
    Chan2b = Chan2bD.Chan2b  
    if X.hermite && Y.hermite; ret.hermite=false; ret.antihermite=true
    elseif X.hermite && Y.antihermite; ret.hermite=true;ret.antihermite=false
    elseif X.antihermite && Y.hermite; ret.hermite=true;ret.antihermite=false;
    else ret.hermite = false;end
    ## zero-body pieces
    if !ret.antihermite
        comm110ss!(X,Y,ret,HFobj)
        comm220ss!(X,Y,ret,HFobj,Chan2b)
    end    
    ## one-body pieces
    @timeit to "comm111" comm111ss!(X,Y,ret)
    @timeit to "comm121" comm121ss!(X,Y,ret,HFobj,Chan1b,dictMono,PandyaObj)    
    @timeit to "comm221" comm221ss!(X,Y,ret,HFobj,Chan1b,Chan2bD,PandyaObj)
 
    ## two-body pieces
    @timeit to "comm122" comm122ss!(X,Y,ret,Chan2b,PandyaObj,to)
    @timeit to "comm222pphh" comm222_pphh_ss!(X,Y,ret,HFobj,Chan2bD,PandyaObj,to)
    @timeit to "comm222ph" comm222ph_ss!(X,Y,ret,HFobj,Chan2bD,dict6j,PandyaObj,to)
    return nothing
end

"""
    BCH_Product(X::Op,Y::Op,Z::Op,tmpOp::Op,Nested::Op,ncomm::Vector{Int64},norms::Vector{Float64},Chan1b::chan1b,Chan2bD::Chan2bD,HFobj::HamiltonianNormalOrdered,dictMono,dict6j,PandyaObj::PandyaObject,to;tol=1.e-4) where Op <:Operator

returns ``Z``  to satisfy: ``e^Z = e^X e^Y``.

``Z`` is calculated with Baker–Campbell–Hausdorff (BCH) formula:
``Z = X + Y + 1/2[X, Y]  + 1/12 [X,[X,Y]] + 1/12 [Y,[Y,X]] -1/24 [Y,[X,[X,Y]]] -1/720 [Y,[Y,[Y,[Y,X]]]] -1/720 [X,[X,[X,[X,Y]]]] +...``

For IMSRG flow of ``H(s)``, ``X=\\eta(s) ds``, ``Y=\\Omega(s)``, and ``Z=\\Omega(s+ds)`` 
"""
function BCH_Product(X::Op,Y::Op,Z::Op,tmpOp::Op,Nested::Op,ncomm::Vector{Int64},norms::Vector{Float64},Chan1b::chan1b,Chan2bD::Chan2bD,HFobj::HamiltonianNormalOrdered,
                     dictMono,dict6j,PandyaObj::PandyaObject,to;tol=1.e-4) where Op <:Operator
    Nested.hermite = true; Nested.antihermite=false
    berfac = [-0.5, 1.0/12.0, 0.0, -1.0/720.0, 0.0, 1.0/30240.0, 0.0, 1.0/1209600.0 ]
    Chan2b = Chan2bD.Chan2b
    p_sps = HFobj.modelspace.p_sps; n_sps = HFobj.modelspace.n_sps   
    ## comm0 term: X+Y -> Z
    Ncomm = 0
    aOp1_p_bOp2_Op3!(X,Y,Z,1.0,1.0,0.0)
    ## comm1 term: 1/2[X,Y]
    ## norm check to determine whther to go further or not 
    Ncomm +=1  
    aOp!(Nested,0.0) ## clear old history
    @timeit to "comm." OpCommutator!(Y,X,Nested,HFobj,Chan1b,Chan2bD,dictMono,dict6j,PandyaObj,to)
    nx  = getNorm(X,p_sps,n_sps,Chan2b)
    sqrt(getNorm1b(X.onebody,p_sps,n_sps)^2 + getNorm2b(X.twobody,Chan2b)^2)
    norm1b = getNorm1b(Nested.onebody,p_sps,n_sps)
    norm2b = getNorm2b(Nested.twobody,Chan2b)
    nxy = sqrt(norm1b^2+norm2b^2)

    if nxy *nx > tol
        Ncomm += 1
        aOp!(tmpOp,0.0)
        @timeit to "comm." OpCommutator!(Nested,X,tmpOp,HFobj,Chan1b,Chan2bD,dictMono,dict6j,PandyaObj,to)
        aOp1_p_bOp2!(tmpOp,Z,1.0/12.0,1.0)
    end
    k=1
    while getNorm(Nested,p_sps,n_sps,Chan2b) > tol && k < 3
        if k<=2 || k%2==1
            aOp1_p_bOp2!(Nested,Z,berfac[k],1.0)
        end
        aOp!(tmpOp,0.0)
        @timeit to "comm." OpCommutator!(Y,Nested,tmpOp,HFobj,Chan1b,Chan2bD,dictMono,dict6j,PandyaObj,to)
        aOp1_p_bOp2!(tmpOp,Nested,1.0,0.0)
        Ncomm += 1
        k +=1
    end

    ### update norms for Omega(s+1)
    ncomm[1] += Ncomm
    norms[1] = getNorm1b(Z.onebody,p_sps,n_sps)
    norms[2] = getNorm2b(Z.twobody,Chan2b)
    return nothing
end

"""
    BCH_Transform(Omega,O0,Hs,tOp,Nested,ncomm,norms,Chan1b,Chan2bD,HFobj,dictMono,dict6j,PandyaObj,to;tol=1.e-9,maxit=50,verbose=false)   

Update ``B`` (assumed to be ``H`` or ``O``) via ``e^ABe^{-A} =B+[A,B]+1/2![A,[A,B]]+1/3![A,[A,[A,B]]]+...``

Note that the `ret` and `tOp` are also overwritten.
`BCH_Transform(nOmega,H0,IMSRGObj,tmpOp,norms,Chan2b)`
e.g.,
`Omega`: ``\\Omega(s+ds)``, `H0`: ``H(s=0)`` or ``O(s=0)``, and `ret`: ``H(s+ds)``
"""
function BCH_Transform(Omega::Op,O0::Op,Hs::Op,tOp::Op,Nested::Op,ncomm,norms,Chan1b::chan1b,Chan2bD::Chan2bD,HFobj::HamiltonianNormalOrdered,
                       dictMono,dict6j,PandyaObj::PandyaObject,to;tol=1.e-9,maxit=100,verbose=false) where Op<:Operator
    if Omega.hermite && O0.hermite; tOp.hermite=false; tOp.antihermite=true
    elseif Omega.hermite && O0.antihermite; tOp.hermite=true;tOp.antihermite=false
    elseif Omega.antihermite && O0.hermite; tOp.hermite=true;tOp.antihermite=false;
    else tOp.hermite = false;end
    Nested.hermite = tOp.hermite; Nested.antihermite = tOp.antihermite
    
    MS = HFobj.modelspace; p_sps = MS.p_sps; n_sps = MS.n_sps
    Chan2b = Chan2bD.Chan2b
    aOp1_p_bOp2!(O0,Nested,1.0,0.0) # B -> Nested
    aOp1_p_bOp2!(Nested,Hs,1.0,0.0) # Os = O0

    nx = getNorm(O0,p_sps,n_sps,Chan2b)
    ny = getNorm(Omega,p_sps,n_sps,Chan2b)
    epsilon = nx * exp(-2*ny) * tol / (2*ny)
    factor_kth = 1.0
    for iter=1:maxit
        aOp!(tOp,0.0)
        factor_kth /= iter
        @timeit to "comm." OpCommutator!(Omega,Nested,tOp,HFobj,Chan1b,Chan2bD,dictMono,dict6j,PandyaObj,to)
        ncomm[1] += 1
        aOp1_p_bOp2!(tOp,Nested,1.0,0.0)
        aOp1_p_bOp2!(Nested,Hs,factor_kth,1.0)
        epsilon *= iter+1
        normNest = getNorm(Nested,p_sps,n_sps,Chan2b)
        if normNest < epsilon; break;end
        if iter == maxit
            println("BCH transform didn't converge after $iter iteration")
        end
    end
    return nothing
end

function comm110ss!(X::Op,Y::Op,ret::Op,HFobj::HamiltonianNormalOrdered) where {Op<:Operator}
    if X.hermite && Y.hermite; return nothing;end
    if X.antihermite && Y.antihermite; return nothing;end
    Z0 = ret.zerobody; x1bs = X.onebody; y1bs = Y.onebody
    sps = HFobj.modelspace.sps; holes = HFobj.modelspace.holes
    for pn = 1:2        
        x1b = x1bs[pn]; y1b = y1bs[pn]
        xyyx = x1b*y1b - y1b*x1b
        for hole in holes[pn]
            ja = sps[hole].j
            hidx = div(hole,2) + hole%2
            Z0[1] += (ja+1) * sps[hole].occ * xyyx[hidx,hidx]
        end
    end
    return nothing
end 

"""
    comm220ss!(X::Op,Y::Op,Z::Op,HFobj::HamiltonianNormalOrdered,Chan2b::Vector{chan2b}) where Op<:Operator

``[X_2,Y_2]_0 = 2 \\sum_{J}(2J+1) \\mathrm{Tr}(X_{hh'pp'}  Y_{pp'hh'}) ``
"""
function comm220ss!(X::Op,Y::Op,Z::Op,HFobj::HamiltonianNormalOrdered,Chan2b::Vector{chan2b}) where Op<:Operator
    sps = HFobj.modelspace.sps
    x2bs = X.twobody; y2bs = Y.twobody
    for ch in eachindex(x2bs)
        x2b = x2bs[ch]; y2b = y2bs[ch]; tbc = Chan2b[ch]
        J = tbc.J; kets = tbc.kets
        NJ = 2*J+1
        zsum = 0.0
        for ib in eachindex(kets)
            a,b = kets[ib]
            na = sps[a].occ; nb = sps[b].occ
            if na * nb == 0.0; continue;end
            for ik in eachindex(kets)
                c,d = kets[ik]
                nc = sps[c].occ; nd = sps[d].occ
                zsum += (x2b[ib,ik]*y2b[ik,ib]-y2b[ib,ik]*x2b[ik,ib])*na*nb *(1-nc)*(1-nd)
            end 
        end
        Z.zerobody[1] += NJ * zsum
    end                                     
    return nothing
end

"""
    comm111ss!(X,Y,ret;inifac=1.0)

returns ``[X_1,Y_1]``
"""
function comm111ss!(X::Op,Y::Op,ret::Op;inifac=1.0) where Op<:Operator
    m1bs = ret.onebody; x1bs = X.onebody; y1bs = Y.onebody
    for pn = 1:2
        m1b = m1bs[pn]; x1b = x1bs[pn]; y1b = y1bs[pn]
        BLAS.gemm!('N','N',1.0,x1b,y1b,inifac,m1b) 
        BLAS.gemm!('N','N',-1.0,y1b,x1b,1.0,m1b)
    end
    return nothing
end 

"""
    comm121ss!(X,Y,ret,HFobj,Chan1b,Chan2b,dictMono,PandyaObj)

returns ``[X_1,Y_2] - [Y_1,X_2]``, whose elements are given as 

`` [X_1,Y_2]_{ij} = \\frac{1}{2j_i+1}\\sum_{ab} (n_a \\bar{n}_b) \\sum_{J} (2J+1) (X_{ab} Y^J_{biaj} - X_{ba} Y^J_{aibj}) ``
"""
function comm121ss!(X::Op,Y::Op,ret::Op,HFobj::HamiltonianNormalOrdered,Chan1b::chan1b,dictMono,PandyaObj::PandyaObject) where Op <:Operator
    x1bs = X.onebody; y1bs = Y.onebody; m1bs = ret.onebody
    x2bs = X.twobody; y2bs = Y.twobody; hZ = ifelse(ret.hermite,1.0,-1.0)
    sps = HFobj.modelspace.sps
    keys6j = PandyaObj.keys6j
    dim1b = size(x1bs[1])[1]
    for pn = 1:2
        y1b = y1bs[pn]; m1b = m1bs[pn]        
        chan1bs = Chan1b.chs1b[pn]
        imin = ifelse(pn==1,1,2)
        imax = ifelse(pn==1,2*dim1b-1,2*dim1b)
        for i = imin:2:imax
            tkey = keys6j[threadid()]
            tkey = @view tkey[1:4]
            idx_i = div(i,2) + i%2
            jdeno = 1.0/(sps[i].j+1.0)
            js = chan1bs[i]
            for j in js
                idx_j = div(j,2) + j%2
                x1y2 = y1x2 = 0.0
                for pn_ab = 1:2
                    Tz = 0
                    if pn == pn_ab ==1; Tz = -2;end
                    if pn == pn_ab ==2; Tz =  2;end
                    x1b = x1bs[pn_ab]; y1b = y1bs[pn_ab]
                    targetDict = dictMono[2+div(Tz,2)]
                    chan1bs_ab = Chan1b.chs1b[pn_ab]
                    amin = ifelse(pn_ab==1,1,2)
                    amax = ifelse(pn_ab==1,2*dim1b-1,2*dim1b)
                    for a = amin:2:amax
                        oa = sps[a]; na = oa.occ; idx_a = div(a,2)+a%2
                        nonzero_bs = chan1bs_ab[a]
                        for b in nonzero_bs
                            ob = sps[b]; nb = ob.occ; idx_b = div(b,2)+b%2
                            #if na + nb != 1; continue;end
                            if (na * nb != 0.0 || na + nb == 0.0); continue;end # for fractional
                            nab = na * (1-nb)
                            if !(x1b[idx_a,idx_b] ==0.0 && x1b[idx_b,idx_a] ==0.0)
                                x1y2 += jdeno * nab* single_121(a,b,i,j,x1b,y2bs,sps,tkey,targetDict)
                            end
                            if !(y1b[idx_a,idx_b] ==0.0 && y1b[idx_b,idx_a] ==0.0)
                                y1x2 += jdeno * nab* single_121(a,b,i,j,y1b,x2bs,sps,tkey,targetDict)
                            end
                        end
                    end
                end
                m1b[idx_i,idx_j] += x1y2 -y1x2
                if idx_i != idx_j; m1b[idx_j,idx_i] = hZ * m1b[idx_i,idx_j];end
            end
        end 
    end
    return nothing
end


"""
    single_121(a,b,i,j,o1b,o2bs,sps,key,targetDict;verbose=false)

calc 121part ``[X_1,Y_2]-[Y_1,X_2]`` for given `i`,`j` and `a`,`b`.

``\\sum_{J} [J]^2 (o_{1,ab}o_{2,biaj} - o_{1,ba}o_{2,aibj})``
"""
function single_121(a,b,i,j,o1b,o2bs,sps,key,targetDict;verbose=false)
    key1 = @view key[1:2]; key2 = @view key[3:4]
    idx_a = div(a,2) + a%2; idx_b = div(b,2) + b%2
    ja = sps[a].j; jb = sps[b].j; ji = sps[i].j; jj = sps[j].j
    ## ajbi term
    r_ajbi = 0.0
    flip_aj = flip_bi = false
    fac_aj = fac_bi = 1.0
    key1[1] = a; key1[2] = j
    if a==j; fac_aj *= sqrt(2.0);end
    if a > j; flip_aj = true;key1[1]=j;key1[2]=a;end
    key2[1] = b; key2[2] = i
    if b==i; fac_bi *= sqrt(2.0);end
    if b > i; flip_bi = true;key2[1]=i;key2[2]=b;end
    fac = fac_aj * fac_bi
    for tmp_ket in targetDict[key1].vals                           
        ch,i_aj,J = tmp_ket
        o2b = o2bs[ch]; NJ = (2*J+1)
        phase_aj = ifelse(flip_aj, (-1)^(div(ja+jj,2)+J+1), 1.0)
        for tmp_bra in targetDict[key2].vals  
            ch_bra,i_bi,Jbra = tmp_bra
            if ch_bra != ch;continue;end
            phase_bi = ifelse(flip_bi, (-1)^(div(jb+ji,2)+Jbra+1),1.0)
            tmp = fac * NJ * phase_aj * phase_bi * o1b[idx_a,idx_b] *o2b[i_bi,i_aj]
            r_ajbi += tmp
        end
    end
    ## aibj term
    r_aibj = 0.0
    flip_ai = flip_bj = false
    fac_ai = fac_bj = 1.0
    key1[1] = a; key1[2] = i
    if a==i; fac_ai *= sqrt(2.0);end
    if a>i; flip_ai = true;key1[1]=i;key1[2]=a;end
    key2[1] = b; key2[2] = j
    if b==j; fac_bj *= sqrt(2.0);end
    if b >j; flip_bj = true;key2[1]=j;key2[2]=b;end
    fac = fac_ai * fac_bj
    for tmp_ket in targetDict[key1].vals
        ch,i_ai,J = tmp_ket
        o2b = o2bs[ch]; NJ = (2*J+1)
        phase_ai = ifelse(flip_ai, (-1)^(div(ja+ji,2)+J+1), 1.0)
        for tmp_bra in targetDict[key2].vals
            ch_bra,i_bj,Jbra = tmp_bra
            if ch_bra != ch;continue;end
            phase_bj = ifelse(flip_bj, (-1)^(div(jb+jj,2)+Jbra+1),1.0)
            r_aibj += fac * NJ * phase_ai * phase_bj * o1b[idx_b,idx_a] *o2b[i_ai,i_bj]
        end
    end
    return r_ajbi - r_aibj 
end

"""
    comm221ss!(X::Op,Y::Op,ret::Op,HFobj::HamiltonianNormalOrdered,Chan1b::chan1b,Chan2bD::Chan2bD,PandyaObj::PandyaObject) where Op<:Operator 

returns ``[X_2,Y_2]_1 - [Y_2,X_2]_1``, whose elements are given as

``[X_2,Y_2]_{ij} = 1/(2[j_i]) \\sum_{abc}\\sum_{J}[J](n'_an'_bn_c-n_an_bn'_c)(X_{2,ciab}Y_{2,abcj}-Y_{2,ciab}X_{2,abcj})``
"""
function comm221ss!(X::Op,Y::Op,ret::Op,HFobj::HamiltonianNormalOrdered,Chan1b::chan1b,Chan2bD::Chan2bD,PandyaObj::PandyaObject) where Op<:Operator
    Chan2b = Chan2bD.Chan2b
    mats_nab = PandyaObj.mats_nab
    mats_nab_bar = PandyaObj.mats_nab_bar
    tMats = PandyaObj.tMat
    hZ = ifelse(ret.hermite,1.0,-1.0)
    x2bs = X.twobody; y2bs = Y.twobody; m1bs = ret.onebody
    sps = HFobj.modelspace.sps
    dim1b = size(m1bs[1])[1]; nchan = length(Chan2b)
    Dict_chidx_from_ketJ = Chan2bD.dict_ch_idx_from_ket
    key6js = PandyaObj.keys6j
    nthre = nthreads()
    for i=1:2*nthre; PandyaObj.copy_1bmat[i] .= 0.0;end
    @threads for ch = 1:nchan
        tbc = Chan2b[ch]; x2 = x2bs[ch]; y2 = y2bs[ch]
        J = tbc.J; Tz = tbc.Tz; tkets = tbc.kets; nket = length(tkets)
        NJ = 2*J+1
        tkey = key6js[threadid()]
        tkey = @view tkey[1:2]
        mat_nab = mats_nab[ch]
        mat_nab_bar = mats_nab_bar[ch]

        tMat =tMats[threadid()]
        XY  = @view tMat[1:nket,1:nket]
        XYbar = @view tMat[nket+1:2*nket,nket+1:2*nket]
        XnY = @view tMat[1:nket,nket+1:2*nket]
        YnX = @view tMat[nket+1:2*nket,1:nket]    

        BLAS.gemm!('N','N',1.0,mat_nab,y2,0.0,XnY);  BLAS.gemm!('N','N',1.0,x2,XnY,0.0,XY)
        BLAS.gemm!('N','N',1.0,mat_nab,x2,0.0,YnX);  BLAS.gemm!('N','N',-1.0,y2,YnX,1.0,XY)

        BLAS.gemm!('N','N',1.0,mat_nab_bar,y2,0.0,XnY); BLAS.gemm!('N','N',1.0,x2,XnY,0.0,XYbar)
        BLAS.gemm!('N','N',1.0,mat_nab_bar,x2,0.0,YnX); BLAS.gemm!('N','N',-1.0,y2,YnX,1.0,XYbar)

        for c = 1:2*dim1b
            oc = sps[c]; jc = oc.j; nc = oc.occ; nbar_c = 1-nc
            pn_ij = 0
            if Tz==0 
                pn_ij = ifelse(oc.tz==-1,2,1)
            else 
                pn_ij = ifelse(oc.tz==-1,1,2)
            end
            tmat = PandyaObj.copy_1bmat[(pn_ij-1)*nthre+threadid()]
            chan1b = Chan1b.chs1b[pn_ij]
            for idx_i = 1:dim1b
                i = 2*(idx_i-1) + pn_ij
                oi = sps[i]; ji = oi.j
                if Tz != sps[i].tz + sps[c].tz; continue;end
                if !tri_check(jc,ji,2*J);continue;end
                jdeno = 1.0/(ji+1)
                for j in chan1b[i]
                    idx_j = div(j,2) + j%2 
                    jj = sps[j].j
                    if Tz != sps[j].tz + sps[c].tz; continue;end            
                    if !tri_check(jc,jj,2*J);continue;end
                    if abs(Tz) ==2 && J%2 == 1 && (c==i||c==j);continue;end
                    tdict = Dict_chidx_from_ketJ[2+div(Tz,2)][J+1]
                    fflip_ci = (-1)^(div(ji+jc,2)+J+1)
                    fflip_cj = (-1)^(div(jj+jc,2)+J+1)
                    phase_ci = ifelse(c>i,fflip_ci,1.0)
                    phase_cj = ifelse(c>j,fflip_cj,1.0)
                    i_ci = i_cj = 0
                    if c<=i
                        tkey[1] = c; tkey[2] = i
                    else    
                        tkey[1] = i; tkey[2] = c
                    end
                    for tmp in tdict[tkey]
                        if tmp[1] == ch
                            i_ci =  tmp[2]
                            break
                        end
                    end
                    if i_ci == 0; continue; end                   
                    if c<=j
                        tkey[1] = c; tkey[2] = j
                    else
                        tkey[1] = j; tkey[2] = c
                    end
                    for tmp in tdict[tkey]
                        if tmp[1] == ch
                            i_cj =  tmp[2]
                            break
                        end
                    end
                    if i_cj == 0; continue; end   
                    sqfac_ci = ifelse(c==i,sqrt(2.0),1.0)
                    sqfac_cj = ifelse(c==j,sqrt(2.0),1.0)
                    sqfac = sqfac_ci*sqfac_cj
                    v_1 = nbar_c * XY[i_ci,i_cj]
                    v_2 = nc * XYbar[i_ci,i_cj]
                    tmat[idx_i,idx_j] += phase_ci * phase_cj *sqfac* jdeno * NJ * (v_1 + v_2)
                    if idx_i != idx_j; tmat[idx_j,idx_i] = hZ * tmat[idx_i,idx_j];end
                end
            end
        end
    end
    for i = 1:2*nthre
        pn_ij = ifelse(i<=nthre,1,2)
        t1b = m1bs[pn_ij]        
        axpy!(1.0,PandyaObj.copy_1bmat[i],t1b)
    end
    return nothing
end

function comm222ph_ss!(X::Op,Y::Op,ret::Op,HFobj::HamiltonianNormalOrdered,Chan2bD::Chan2bD,dict6j,PandyaObj::PandyaObject,to) where Op<:Operator
    hy = ifelse(Y.hermite,1.0,-1.0)
    hz = ifelse(ret.hermite,1.0,-1.0)
    MS = HFobj.modelspace; sps = MS.sps
    numbers = PandyaObj.numbers    
    Chan2b_Pandya = PandyaObj.Chan2b
    nchPandya = length(Chan2b_Pandya)    
    phkets = PandyaObj.phkets
    PhaseMats = PandyaObj.PhaseMats
    Zbars = PandyaObj.Zbars; XYbars = PandyaObj.XYbars; keys6j = PandyaObj.keys6j
    @threads for ich=1:nchPandya
        tmpMat = PandyaObj.tMat[threadid()]
        key6j = keys6j[threadid()]
        Zbar = Zbars[ich]; Zbar .= 0.0
        tbc_cc = Chan2b_Pandya[ich]        
        tnumbers = numbers[ich]
        ch_cc,nKets_cc,nhh,nph = tnumbers        
        nph_kets = nhh + nph 
        Xbar_ph = @view XYbars[threadid()][1][1:nKets_cc,1:2*nph_kets]
        Ybar_ph = @view XYbars[threadid()][2][1:2*nph_kets,1:nKets_cc]        
        if nKets_cc * nph_kets !=0
            Xbar_ph .= 0.0; Ybar_ph .= 0.0
            DoPandyaTransformation(X, Xbar_ph, tbc_cc, Chan2bD, HFobj, tnumbers,dict6j,key6j,"T") 
            DoPandyaTransformation(Y, Ybar_ph, tbc_cc, Chan2bD, HFobj, tnumbers,dict6j,key6j,"N") 
            PhaseMat = PhaseMats[ich]; PhaseMat .= 1.0
            tkets = tbc_cc.kets
            for iket = 1:nKets_cc
                p,q = tkets[iket]
                jp = sps[p].j; jq = sps[q].j
                if div(jp+jq,2) % 2 != 0;continue;end
                tv = @view PhaseMat[iket,:]; tv .*= -1.0
                tv = @view PhaseMat[:,iket]; tv .*= -1.0
            end
            PhaseMatY = @view PhaseMat[phkets[ich],:]
            Zleft  = @view Zbar[:,1:nKets_cc]
            Zright = @view Zbar[:,nKets_cc+1:2*nKets_cc]
            calcZbar!(Xbar_ph,Ybar_ph,PhaseMat,PhaseMatY,tmpMat,hy,nph_kets,nKets_cc,Zleft,Zright,hz)
        end
    end
    @timeit to "AddInv" AddInvPandya!(Zbars,ret,Chan2bD,dict6j,PandyaObj,sps,to) 
    return nothing
end 

"""
    calcZbar!(Xbar,Ybar,PhaseMat,PhaseMatY,tmpMat,hy,nph_kets,nKets_cc,Zlefthalf,Zrighthalf,hz)

- `Xbar`: (`nKets_cc`, 2*`nph_kets`) matrix or SubArray
- `Ybar`: (2*`nph_kets`,`nKets_cc`) matrix or SubArray
- `PhaseMatY`: (`nph_kets`,`nKets_cc`) matrix or SubArray
- `Zbar`: (`nKets_cc`, 2*`nKets_cc`) matrix or SubArray
"""
function calcZbar!(Xbar,Ybar,PhaseMat,PhaseMatY,tmpMat,hy,nph_kets,nKets_cc,Zlefthalf,Zrighthalf,hz) 
    BLAS.gemm!('N','N',1.0,Xbar,Ybar,0.0,Zlefthalf)
    ru = @view Ybar[nph_kets+1:2*nph_kets,:] 
    rl = @view Ybar[1:nph_kets,:]
    upper = @view tmpMat[1:nph_kets,1:nKets_cc]; upper .= (ru .* PhaseMatY)
    lower = @view tmpMat[nph_kets+1:2*nph_kets,1:nKets_cc]; lower .= (rl .* PhaseMatY)
    tmpMatr = @view tmpMat[1:2*nph_kets,1:nKets_cc]
    BLAS.gemm!('N','N',hy,Xbar,tmpMatr,0.0,Zrighthalf)

    tM = @view tmpMat[1:nKets_cc,1:nKets_cc];
    tM .= Zrighthalf' .* PhaseMat; Zrighthalf .+= tM

    tM .= 0.0; axpy!(hz,Zlefthalf',tM); Zlefthalf .+= tM
    return nothing
end 

function AddInvPandya!(Zbars::Vector{Matrix{Float64}},ret::Operator,Chan2bD::Chan2bD,dict6j,PandyaObj::PandyaObject,sps,to)
    Chan2b =Chan2bD.Chan2b
    Chan2b_Pandya = PandyaObj.Chan2b
    dict_2b_ch = Chan2bD.dict_ch_JPT
    dict_ch2ich = PandyaObj.dict_ch2ich
    numbers_addinv = PandyaObj.numbers_addinv
    tdict_ichidx = PandyaObj.dict_ich_idx_from_ketcc
    hermite = ret.hermite
    ofst1 = 1000
    @inbounds @threads for ich in eachindex(numbers_addinv)
        ch = numbers_addinv[ich][1]
        tbc = Chan2b[ch]; tkets = tbc.kets; J = tbc.J
        nKets = length(tkets)
        if nKets == 0; continue;end
        key6j = PandyaObj.keys6j[threadid()]
        Z2b = ret.twobody[ch]
        tdict = dict6j[J+1]
        @inbounds for ibra = 1:nKets
            i,j = tkets[ibra]
            li = sps[i].l; tzi = sps[i].tz
            ji = sps[i].j; jj = sps[j].j
            ketmin = ifelse(hermite,ibra,ibra+1)            
            @inbounds for iket = ketmin:nKets
                k,l = tkets[iket]
                lk = sps[k].l; jk = sps[k].j; tzk = sps[k].tz
                ll = sps[l].l; jl = sps[l].j; tzl = sps[l].tz
                commij = commji = 0.0 
                prty_cc = (-1)^((li+ll)%2)
                aTz_cc = abs(tzi+tzl)
                Jpmin = max(abs(div(ji-jl,2)),abs(div(jk-jj,2)))
                Jpmax = min(div(ji+jl,2), div(jk+jj,2))
                for Jpr = Jpmin:Jpmax
                    nkey = get_nkey_from_key6j(ji,jj,jk,jl,Jpr;ofst_unit=ofst1)
                    sixj = tdict[nkey]
                    if abs(sixj) < 1.e-8;continue;end
                    key_JPT = @view key6j[1:3]
                    key_JPT[1] = aTz_cc; key_JPT[2] = prty_cc; key_JPT[3] = Jpr
                    ich_cc = dict_ch2ich[ dict_2b_ch[key_JPT].Vch ]
                    nketcc = length(Chan2b_Pandya[ich_cc].kets)
                    flipil = ifelse(i>l,true,false); flipkj = ifelse(k>j,true,false)
                    tkey_il = @view key6j[1:2]; tkey_kj = @view key6j[3:4] 
                    update_key!(i,l,flipil,tkey_il)
                    update_key!(k,j,flipkj,tkey_kj)
                    
                    idx_il = tdict_ichidx[ich_cc][tkey_il[1] * ofst1 + tkey_il[2]]
                    idx_kj = tdict_ichidx[ich_cc][tkey_kj[1] * ofst1 + tkey_kj[2]] + ifelse(k>j,nketcc,0)

                    phaseil = (-1)^( div(ji+jl,2)+Jpr+1 )
                    phase = ifelse(flipil,phaseil,1.0)
                    me1 = Zbars[ich_cc][idx_il,idx_kj] 
                    commij -= (2*Jpr+1) * sixj * me1 *phase
                end
            
                if k==l
                    commji = commij
                elseif i==j
                    commji = (-1)^((div(ji+jj,2)+div(jk+jl,2))%2) *commij
                else
                    prty_cc = (-1)^((li+lk)%2)
                    aTz_cc = abs(tzi+tzk)
                    Jpmin = max( abs(div(jj-jl,2)), abs(div(jk-ji,2)))
                    Jpmax = min( div(jj+jl,2), div(jk+ji,2))
                    for Jpr = Jpmin:Jpmax                            
                        nkey = get_nkey_from_key6j(jj,ji,jk,jl,Jpr;ofst_unit=ofst1)
                        sixj = tdict[nkey]

                        if abs(sixj) < 1.e-8;continue;end
                      
                        key_JPT = @view key6j[1:3]
                        key_JPT[1] = aTz_cc; key_JPT[2] = prty_cc; key_JPT[3] = Jpr
                        ch_cc = dict_2b_ch[key_JPT].Vch
                        ich_cc = dict_ch2ich[ch_cc]   
                        nketcc = length(Chan2b_Pandya[ich_cc].kets)
                        flipik = ifelse(i>k,true,false); fliplj = ifelse(l>j,true,false)
                        tkey_ik = @view key6j[1:2]; tkey_lj = @view key6j[3:4] 
                        update_key!(i,k,flipik,tkey_ik)
                        update_key!(l,j,fliplj,tkey_lj)

                        idx_ik = tdict_ichidx[ich_cc][tkey_ik[1]*ofst1+tkey_ik[2]] 
                        idx_lj = tdict_ichidx[ich_cc][tkey_lj[1]*ofst1+tkey_lj[2]]+ ifelse(l>j,nketcc,0)                        

                        phaseik = (-1)^( div(ji+jk,2)+Jpr+1)
                        phase = ifelse(flipik,phaseik,1.0)
                        me1 = Zbars[ich_cc][idx_ik,idx_lj] * phase 
                        commji -= (2*Jpr+1) * sixj * me1
                    end                
                end
                d_ij = delta(i,j);d_kl = delta(k,l)
                tnorm = ifelse(d_ij==d_kl, 1.0+delta(i,j),sqrt(2.0))
                phase = (-1)^(div(jk+jl,2)-J)
                Z2b[ibra,iket] -= (commij - phase * commji) /tnorm
                if hermite 
                    Z2b[iket,ibra] = Z2b[ibra,iket]
                elseif iket != ibra 
                    Z2b[iket,ibra] = -Z2b[ibra,iket]
                end
            end 
        end
    end
    return nothing
end
function update_key!(i,j,flip::Bool,tkey)
    if flip
        tkey[1] = j; tkey[2] = i
    else
        tkey[1] = i; tkey[2] = j
    end
    return nothing
end

function getlocalidx(i,j,tkets_cc;ofst=true)
    flip = false  
    ti = i; tj = j
    if i > j 
        ti = j; tj = i
        flip = true
    end
    tidx = 0
    for (idx,ket) in enumerate(tkets_cc)
        if ket[1] == ti && ket[2] == tj
            tidx = idx
            break
        end  
    end
    if ofst 
        return tidx + ifelse(i>j,length(tkets_cc),0),flip
    else
        return tidx,flip
    end
end

function prep122(HFobj::HamiltonianNormalOrdered,Chan1b::chan1b,Chan2bD::Chan2bD)
    Chan2b = Chan2bD.Chan2b
    idx_dict = Chan2bD.dict_idx_from_chket
    MS = HFobj.modelspace; sps = MS.sps
    nch = length(Chan2b)
    tkey = zeros(Int64,2)
    util122 = [ single_util122[] for ch=1:nch]
    for ch = 1:nch
        J = Chan2b[ch].J; tkets = Chan2b[ch].kets
        dim = length(tkets)
        targetDict = idx_dict[ch]
        @inbounds for idx_ij = 1:dim                       
            i,j = tkets[idx_ij]
            pn_i = 2-i%2; pn_j = 2-j%2
            ji = sps[i].j; jj = sps[j].j
            phaseij = (-1)^(div(ji+jj,2)+J+1)
            ind1_ia = Int64[]; ind1_ja = Int64[ ]
            ind2_ai = Int64[]; ind2_aj = Int64[ ]
            factor_ia = Float64[]; factor_ja = Float64[]
            for a in Chan1b.chs1b_redundant[pn_i][i]
                if a <= j;tkey[1] = a; tkey[2] = j
                else;tkey[1] = j; tkey[2] = a; end
                ind2 = get(targetDict,tkey,0)
                if ind2 == 0;continue;end
                idx_a = div(a,2)+a%2
                push!(ind1_ia,idx_a); push!(ind2_aj,ind2)
                ###
                phaseaj = (-1)^(div(sps[a].j+jj,2)+J+1)
                if phaseaj != phaseij; println("warn!"); exit();end
                #### just in case...
                tfac = ifelse(a>j,phaseij,1.0) * ifelse(a==j,sqrt(2.0),1.0)
                push!(factor_ia, tfac)
            end
            if i!=j
                for a in Chan1b.chs1b_redundant[pn_j][j]
                    if a <= i
                        tkey[1] = a; tkey[2] = i
                    else
                        tkey[1] = i; tkey[2] = a
                    end
                    ind2 = get(targetDict,tkey,0)
                    if ind2==0;continue;end
                    idx_a = div(a,2)+a%2
                    push!(ind1_ja,idx_a)
                    push!(ind2_ai,ind2)
                    ### just in case...
                    phaseai = (-1)^(div(sps[a].j+ji,2)+J+1)
                    if phaseai != phaseij; println("warn!"); exit();end
                    ###
                    tfac = ifelse(i>a,phaseij,1.0) * ifelse(i==a,sqrt(2.0),1.0)
                    push!(factor_ja,tfac)
                end
            end    
            if i==j
                factor_ia ./= sqrt(2.0)
                factor_ja ./= sqrt(2.0)
            end
            ind2s = vcat(ind2_aj,ind2_ai) 
            dimia = length(factor_ia); dimja = length(factor_ja)
            mia = zeros(Float64,dimia,dimia); for i=1:dimia; mia[i,i] = factor_ia[i];end
            mja = zeros(Float64,dimja,dimja); for i=1:dimja; mja[i,i] = factor_ja[i];end
            push!(util122[ch],single_util122(ind1_ia,ind1_ja,ind2s,mia,mja))
        end
    end
    return util122
end
function comm122ss!(X::Op,Y::Op,ret::Op,Chan2b::Vector{chan2b},PandyaObj::PandyaObject,to) where Op<:Operator
    hZ = ifelse(ret.hermite,1.0,-1.0)
    x1 = X.onebody; x2bs = X.twobody
    y1 = Y.onebody; y2bs = Y.twobody; m2bs = ret.twobody
    nch = length(Chan2b)
    util122 = PandyaObj.util122
    tMats = PandyaObj.tMat    
    XYbars = PandyaObj.XYbars
    @threads for ch = 1:nch
        x2 = x2bs[ch]; y2 = y2bs[ch]; z2 = m2bs[ch]
        tkets = Chan2b[ch].kets; dim = length(tkets)
        tmp_util = util122[ch]     
        tM = tMats[threadid()]   
        w2  = @view tM[1:dim,1:dim]
        M_L = XYbars[threadid()][1]
        M_R = XYbars[threadid()][2]
        tx1_i = @view tM[1:dim,dim+1]
        tx1_j = @view tM[1:dim,dim+2]
        ty1_i = @view tM[1:dim,dim+3]
        ty1_j = @view tM[1:dim,dim+4]
        @inbounds for idx_ij = 1:dim
            i,j = tkets[idx_ij]
            idx_i = div(i,2)+i%2; pn_i = 2-i%2
            idx_j = div(j,2)+j%2; pn_j = 2-j%2

            util = tmp_util[idx_ij]
            ind1_ia = util.ind1_ia
            ind1_ja = util.ind1_ja
            idxs = util.ind2s
            factor_ia = util.factor_ia
            factor_ja = util.factor_ja

            lidxs = length(idxs)
            y2col = @view y2[:,idxs]
            x2col = @view x2[:,idxs]
            w2col = @view w2[:,idx_ij]
            x1_i = @view tx1_i[1:length(ind1_ia)]; x1_i .= @view x1[pn_i][ind1_ia,idx_i]
            x1_j = @view tx1_j[1:length(ind1_ja)]; x1_j .= @view x1[pn_j][ind1_ja,idx_j]
            y1_i = @view ty1_i[1:length(ind1_ia)]; y1_i .= @view y1[pn_i][ind1_ia,idx_i]
            y1_j = @view ty1_j[1:length(ind1_ja)]; y1_j .= @view y1[pn_j][ind1_ja,idx_j]
            left = @view M_L[1:dim,1:2*lidxs]
            lL = @view left[1:dim,1:lidxs]; lL .= y2col
            lR = @view left[1:dim,lidxs+1:2*lidxs]; lR .= x2col

            right = @view M_R[1:2*lidxs,1:1]
            lxi = length(x1_i); lxj = length(x1_j); lyi = length(y1_i); lyj = length(y1_j)
            rtmp = @view right[1:lxi,1:1]
            BLAS.gemm!('N','N',1.0,factor_ia,x1_i,0.0,rtmp) 
            rtmp = @view right[lxi+1:lxi+lxj,1:1]
            BLAS.gemm!('N','N',1.0,factor_ja,x1_j,0.0,rtmp)
            rtmp = @view right[lxi+lxj+1:lxi+lxj+lyi,1:1]
            BLAS.gemm!('N','N',-1.0,factor_ia,y1_i,0.0,rtmp)
            rtmp = @view right[lxi+lxj+lyi+1:lxi+lxj+lyi+lyj,1:1]
            BLAS.gemm!('N','N',-1.0,factor_ja,y1_j,0.0,rtmp)
            BLAS.gemm!('N','N',1.0,left,right,0.0,w2col)
            if i==j; w2col .*= 2.0;end
        end
        z2 .-= w2 
        axpy!(-hZ,w2',z2) 
    end
    return nothing
end 

function comm222_pphh_ss!(X::Op,Y::Op,ret::Op,HFobj::HamiltonianNormalOrdered,Chan2bD::Chan2bD,PandyaObj::PandyaObject,to) where Op<:Operator
    Chan2b = Chan2bD.Chan2b
    m2bs = ret.twobody; x2bs = X.twobody; y2bs = Y.twobody
    nch = length(Chan2b)
    Mats_hh = PandyaObj.Mats_hh; Mats_pp = PandyaObj.Mats_pp; Mats_ph = PandyaObj.Mats_ph
    @threads for ch = 1:nch
        tmpMat = PandyaObj.tMat[threadid()]
        x2 = x2bs[ch]; y2 = y2bs[ch]; z2 = m2bs[ch]
        tbc = Chan2b[ch]; kets = tbc.kets; nkets = length(kets)
        ppidx = get(HFobj.modelspace.spaces.pp,ch,Int64[]); npp = length(ppidx) 
        hhidx = get(HFobj.modelspace.spaces.hh,ch,Int64[]); nhh = length(hhidx) 
        phidx = get(HFobj.modelspace.spaces.ph,ch,Int64[]); nph = length(phidx)
        Matpp = Mats_pp[ch]; Mathh = Mats_hh[ch]; Matph = Mats_ph[ch]
        # pp term, a&b are particle
        if npp != 0
            XppL = @view x2[:,ppidx]; XppR = @view x2[ppidx,:]    
            YppL = @view y2[:,ppidx]; YppR = @view y2[ppidx,:]
            Mpp_1 = @view tmpMat[1:nkets,1:npp]
            Mpp_2 = @view tmpMat[1:nkets,npp+1:2*npp]
            tM1 = @view tmpMat[nkets+1:2*nkets,1:nhh]
            Mpp_1 .= XppL; Mpp_2 .= YppR'; BLAS.gemm!('N','T', 1.0,Mpp_1,Mpp_2,1.0,z2)
            Mpp_1 .= YppL; Mpp_2 .= XppR'; BLAS.gemm!('N','T',-1.0,Mpp_1,Mpp_2,1.0,z2) 
        end        
        if nph != 0
            XphL = @view x2[:,phidx]; XphR = @view x2[phidx,:]    
            YphL = @view y2[:,phidx]; YphR = @view y2[phidx,:]
            Mph_1 = @view tmpMat[1:nkets,1:nph]
            Mph_2 = @view tmpMat[1:nkets,nph+1:2*nph]
            tM1 = @view tmpMat[nkets+1:2*nkets,1:nph]
            Mph_1 .= XphL; Mph_2 .= YphR' ; BLAS.gemm!('N','N',1.0,Mph_1,Matph,0.0,tM1); BLAS.gemm!('N','T', 1.0,tM1,Mph_2,1.0,z2)
            Mph_1 .= YphL; Mph_2 .= XphR' ; BLAS.gemm!('N','N',-1.0,Mph_1,Matph,0.0,tM1); BLAS.gemm!('N','T',1.0,tM1,Mph_2,1.0,z2)
        end        
        ##hh term, a&b are hole
        if nhh != 0 
            XhhL = @view x2[:,hhidx]; XhhR = @view x2[hhidx,:]
            YhhL = @view y2[:,hhidx]; YhhR = @view y2[hhidx,:]
            Mhh_1 = @view tmpMat[1:nkets,1:nhh]       
            Mhh_2 = @view tmpMat[1:nkets,nhh+1:2*nhh] 
            tM = @view tmpMat[nkets+1:2*nkets,1:nhh]
            Mhh_1 .= XhhL; Mhh_2 .= YhhR' ; BLAS.gemm!('N','N',-1.0,Mhh_1,Mathh-Matpp,0.0,tM); BLAS.gemm!('N','T',1.0,tM,Mhh_2,1.0,z2)
            Mhh_1 .= YhhL; Mhh_2 .= XhhR' ; BLAS.gemm!('N','N', 1.0,Mhh_1,Mathh-Matpp,0.0,tM); BLAS.gemm!('N','T',1.0,tM,Mhh_2,1.0,z2)    
        end
    end
    return nothing
end
