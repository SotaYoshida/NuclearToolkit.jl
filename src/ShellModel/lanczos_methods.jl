function TRL(vks,uks,Tmat,k,
             pbits,nbits,jocc_p,jocc_n,
             SPEs,ppinfo,nninfo,bis,bfs,block_tasks,
             p_NiNfs,n_NiNfs,Mps,delMs,Vpn,
              eval_jj,oPP,oNN,oPNu,oPNd,Jidxs,
             tdims,num_ev,num_history,lm,ls,en,tol,to,doubleLanczos=false,checkorth=false)
    mdim = tdims[end]
    TF=[false]

    lnJ = ls
    Jvs = [zeros(Float64,mdim) for i=1:lnJ]
    Jvret = [zeros(Float64,mdim)]
    Jmat = zeros(Float64,lnJ,lnJ)
    Jvs[1] .= vks[1]
    Jtol = 1.e-6
    JTF = [false];beta_J = 1.0
    if doubleLanczos
        Jlanczos(Jvs,Jmat,TF,JTF,Jtol,Jvret,pbits,nbits,tdims,eval_jj,
                 Jidxs,oPP,oNN,oPNu,oPNd,beta_J,to)
        vks[1] .= Jvs[1]
        for i=1:lnJ;Jvs[i] .=0.0;end
    end
    elit = 1
    while TF[1]==false
        for it = k:lm-1
            vk =vks[it]; vkp1 =vks[it+1]
            @timeit to "operate H" begin
                operate_H!(vk,vkp1,
                           pbits,nbits,
                           jocc_p,jocc_n,SPEs,
                           ppinfo,nninfo,
                           tdims,bis,bfs,block_tasks,
                           p_NiNfs,n_NiNfs,Vpn,Mps,delMs,to)
            end
            talpha =  dot(vk,vkp1)
            Tmat[it,it] = talpha
            diagonalize_T!(it,num_ev,Tmat,en,num_history,TF,tol)            
            if it == mdim; TF[1]=true;end
            if TF[1];elit=it;break;end
            svks = @views vks[1:it]
            @timeit to "ReORTH" ReORTH(it,vkp1,svks)
            tbeta = sqrt(dot(vkp1,vkp1))
            vkp1 .*= 1.0/tbeta
            if checkorth; Check_Orthogonality(it,vks,en); end
            Tmat[it+1,it] = tbeta; Tmat[it,it+1] = tbeta
            if doubleLanczos
                if tbeta < Jtol;TF[1]=true;elit=it;break;end
                @timeit to "JJ lanczos" begin
                    for i =2:length(Jvs);Jvs[i] .=0.0;end
                    Jmat .= 0.0;JTF[1] = false
                    Jvs[1] .= vkp1
                    Jlanczos(Jvs,Jmat,TF,JTF,Jtol,Jvret,
                             pbits,nbits,tdims,eval_jj,
                             Jidxs,oPP,oNN,oPNu,oPNd,beta_J,to)
                    vkp1 .= Jvs[1]                            
                    if TF[1];elit=it;break;end
                end
            end
        end
        if TF[1] == false
            @timeit to "Restart" begin
                ThickRestart(vks,uks,Tmat,lm,ls)
            end
        end
        k = ls+1
    end
    return elit
end

function Jlanczos(Jvs,Jmat,TF,JTF,Jtol,Jvret,
                  pbits,nbits,tdims,eval_jj,
                  Jidxs,oPP,oNN,oPNu,oPNd,beta_J,to)
    mdim = tdims[end]                 
    lnJ = length(Jvs)
    eljit=1; k = 1; inow=k
    while JTF[1] == false
        for it = k:lnJ-1
            inow = it
            vk = Jvs[it]; vkp1 = Jvs[it+1]
            operate_J!(vk,vkp1,pbits,nbits,tdims,
                       Jidxs,oPP,oNN,oPNu,oPNd,beta_J)
            axpy!(eval_jj,vk,vkp1)## vkp1 .+= eval_jj .* vk  
            Jmat[it,it] = dot(vk,vkp1)
            teval = eigvals(@views Jmat[1:it,1:it])[1]
            if abs(teval-eval_jj) < Jtol
                eljit=it;JTF[1] = true;break                
            end
            if JTF[1];eljit=it;break;end
            ReORTH(it,vkp1,Jvs)
            beta = sqrt(dot(vkp1,vkp1))
            if beta < 1.e-4;eljit=it;TF[1]=true;JTF[1]=true;break;end
            vkp1 .*= 1.0/beta
            Jmat[it+1,it] = beta; Jmat[it,it+1] = beta
            eljit = it
        end
        if JTF[1]==false
            ThickRestart_J(Jvs,Jvret,Jmat,eljit+1,1,eval_jj,Jtol)
            k=2
        end
    end
    if inow > k
        ThickRestart_J(Jvs,Jvret,Jmat,eljit+1,1,eval_jj,Jtol)
    end
    return nothing
end

function operate_H!(wf,twf,
                    pbits,nbits,
                    jocc_p,jocc_n,
                    SPEs,ppinfo,nninfo,
                    tdims,bis,bfs,block_tasks,
                    p_NiNfs,n_NiNfs,Vpn,
                    Mps,delMs,
                    to=nothing) 
    #lblock = length(pbits)
    @inbounds @threads for bi in block_tasks
        if bi==0; continue;end #empty job
        ret = [0,0,0]; ret2 = [0,0,0]
        idim = tdims[bi]
        l_Np = length(pbits[bi])
        l_Nn = length(nbits[bi])
        offset = idim -l_Nn

        #@timeit to "pn1" begin
        @inbounds for (bfidx,bf) in enumerate(bfs[bi])
            bisearch!(delMs,Mps[bf]-Mps[bi],ret) #!! bf = bfs[bi][j]
            Vs = Vpn[ret[1]]
            fdim = tdims[bf]
            l_Nn_f = length(nbits[bf])
            p_NiNf = p_NiNfs[bi][bfidx]
            n_NiNf = n_NiNfs[bi][bfidx]
            off_f = fdim-l_Nn_f
            @inbounds for (nth,V) in enumerate(Vs)
                Npifs = p_NiNf[nth]; Nnifs = n_NiNf[nth]
                @inbounds for Npif in Npifs
                    tMi = offset+ Npif.i .* l_Nn
                    tMf = off_f + Npif.f .* l_Nn_f
                    phase_p = Npif.phase
                    @inbounds for Nnif in Nnifs
                        Mi = tMi + Nnif.i; Mf = tMf + Nnif.f
                        phase_n = Nnif.phase
                        twf[Mi] += ifelse(phase_p!=phase_n,-V,V) .* wf[Mf]
                    end
                end
            end
        end
        #end
        #@timeit to "pn2" begin
        @inbounds for j = 1:length(bis[bi])-1    #### tbf = bi !!!
            tbi = bis[bi][j]
            bisearch!(delMs,Mps[bi]-Mps[tbi],ret)
            bisearch_ord!(bfs[tbi],bi,ret2)
            fdim=tdims[tbi]
            l_Nn_i=length(nbits[tbi])
            p_NiNf = p_NiNfs[tbi][ret2[1]]
            n_NiNf = n_NiNfs[tbi][ret2[1]]
            off_f = fdim - l_Nn_i
            @inbounds for (nth,V) in enumerate(Vpn[ret[1]])
                Npifs = p_NiNf[nth]
                Nnifs = n_NiNf[nth]
                @inbounds for Npif in Npifs
                    tMi = off_f  + Npif.i *l_Nn_i #idim <-> fdim
                    tMf = offset + Npif.f*l_Nn
                    phase_p = Npif.phase
                    @inbounds for Nnif in Nnifs
                        Mi = tMi + Nnif.i; Mf = tMf + Nnif.f
                        phase_n = Nnif.phase
                        twf[Mf] += ifelse(phase_p!=phase_n,-V,V) .* wf[Mi]
                    end
                end
            end
        end
        #end
        #@timeit to "pp/nn" begin
        ### pp/nn interaction
        @inbounds for (Npi,tinfo) in enumerate(ppinfo[bi])
            tMi = offset + Npi*l_Nn
            @inbounds for (jj,tmp) in enumerate(tinfo)
                tMf = offset + l_Nn*tmp.f
                fac = tmp.coef
                @inbounds for nidx = 1:l_Nn
                    Mi = tMi + nidx; Mf = tMf + nidx
                    twf[Mf] += fac .* wf[Mi]
                    twf[Mi] += fac .* wf[Mf]
                end
            end
        end
        @inbounds for (Nni,tinfo) in enumerate(nninfo[bi])
            tMi = offset + Nni
            @inbounds for (jj,tmp) in enumerate(tinfo)
                tMf = offset + tmp.f
                fac = tmp.coef
                @inbounds for pidx = 1:l_Np
                    Mi = tMi + pidx*l_Nn; Mf = tMf + pidx*l_Nn
                    twf[Mf] += fac .* wf[Mi]
                    twf[Mi] += fac .* wf[Mf]
                end
            end
        end
        #end
        #@timeit to "1b" begin
        ### one-body operator
        @inbounds for pidx = 1:l_Np
            tMi = offset + pidx*l_Nn
            @inbounds for nidx =1:l_Nn
                Mi = tMi + nidx
                twf[Mi] += (dot(SPEs[1],jocc_p[bi][pidx])+
                            dot(SPEs[2],jocc_n[bi][nidx])) .* wf[Mi]
            end
        end
        #end
    end
    return nothing
end

function diagonalize_T!(k::Int64,num_ev::Int64,
                        Tmat,
                        en::Array{Array{Float64,1}},
                        num_history::Int64,
                        TF::Array{Bool,1},
                        tol::Float64) 
    for ith = num_history:-1:2
        en[ith] .= en[ith-1]
    end
    n = minimum([num_ev,k])
    @views en[1][1:n] .= eigvals(@views Tmat[1:k,1:k] )[1:n]
    if all( en[2]-en[1] .< tol) ; TF[1]=true; end
    return nothing
end

function TRBL(q,vks,uks,Tmat,Beta_H,pbits,nbits,jocc_p,jocc_n,SPEs,
               ppinfo,nninfo,bis,bfs,block_tasks,
               p_NiNfs,n_NiNfs,Mps,delMs,Vpn,tdims,
               eval_jj,oPP,oNN,oPNu,oPNd,Jidxs,
               num_ev,num_history,lm,ls_sub,en,tol,to,
               doubleLanczos=false)    
    ls = q*ls_sub   
    mdim = tdims[end]
    TF=[false]
    elit = 1
    inow = 1; itmin = 1; itmax = div(lm,q)-1
    rescount = 0
    Vt = zeros(Float64,q,mdim)
    R  = zeros(Float64,q,q)
    Beta_J = zeros(Float64,q,q)

    lnJ = 20 # maximum([q*4,10])
    Jvs = [zeros(Float64,q,mdim) for i=1:lnJ]
    Jvret = [zeros(Float64,q,mdim) ]
    Jmat = zeros(Float64,lnJ*q,lnJ*q)
    JTF = [false];Jtol=1.e-6
    U = zeros(Float64,mdim,lnJ*q)
    Mat = zeros(Float64,mdim,q)
    if doubleLanczos
        @timeit to "JJ lanczos" begin
            Jvs[1] .= vks[1]
            bl_JJ_Lanczos(q,Jvs,Jmat,Vt,R,Beta_J,JTF,Jtol,Jvret,
                          pbits,nbits,tdims,eval_jj,
                          Jidxs,oPP,oNN,oPNu,oPNd,to,Beta_H,U,Mat)
            vks[1] .= Jvs[1];Beta_J .= 0.0;JTF[1]=false
            for i=1:lnJ;Jvs[i] .= 0.0;end
            V = vks[1]
        end
    end
    while TF[1]==false
        for it = itmin:itmax
            inow = it
            V  = vks[it]
            HV = vks[it+1]
            @timeit to "bl_operateH" begin
                bl_operate_H!(q,V,HV,pbits,nbits,
                              jocc_p,jocc_n,SPEs,ppinfo,nninfo,
                              tdims,bis,bfs,block_tasks,
                              p_NiNfs,n_NiNfs,Vpn,Mps,delMs,to)
            end            
            BLAS.gemm!('N','T',1.0,V,HV,0.0,R) ##mul!(R,V,HV')
            if issymmetric(R) == false
                @inbounds for i=1:q;for j=i:q
                    R[i,j] = R[j,i]
                end;end
            end
            @views Tmat[q*it-q+1:q*it,q*it-q+1:q*it] .= R
            diagonalize_T!(it*q,num_ev,Tmat,en,num_history,TF,tol)
            #print_vec("En TR(d)BL $it ",en[1])          
            BLAS.gemm!('N','N',-1.0,R,V,1.0,HV)#mul!(HV,R,V,-1.0,1.0)
            s_vks = @views vks[1:it-1]
            @timeit to "ReORTH" bl_ReORTH(q,it,mdim,s_vks,HV,Vt,R)            
            bl_QR!(HV',Beta_H,mdim,q)#myQR!(HV',Beta_H,mdim,q)
            if doubleLanczos
                tnorm = norm(Diagonal(Beta_H),Inf)
                if tnorm < Jtol
                    println("Hbn norm $tnorm");TF[1]=true;elit=it;break
                end            
                @timeit to "JJ lanczos" begin
                    for i=2:lnJ; Jvs[i] .= 0.0; end
                    JTF[1] = false; Jmat .= 0.0;Jvs[1] .= HV
                    Vt .= 0.0; R .= 0.0
                    bl_JJ_Lanczos(q,Jvs,Jmat,Vt,R,Beta_J,JTF,Jtol,Jvret,pbits,nbits,
                                  tdims,eval_jj,Jidxs,oPP,oNN,oPNu,oPNd,to,Beta_H,U,Mat)
                    HV .= Jvs[1]; Beta_J .=0.0
                end
            end
            add_bl_T!(q,it,Tmat,Beta_H)
            if TF[1];elit = it;break;end            
        end
        if TF[1] == false
            @timeit to "bl_Restart" begin
                bl_ThickRestart(q,vks,uks,Beta_H,Tmat,inow,ls_sub,mdim,Vt)
            end
            itmin = ls_sub + 1
        end        
    end
    return elit
end

function bl_ReORTH(q,it,mdim,vks,HV,Vt,R)
    #LinearAlgebra.jl:  mul!(C, A, B, α, β) -> C: ABalpha + Cbeta
    #BLAS.gemm!(tA, tB, alpha, A, B,beta,C)   C := alpha*A*B + beta*C
    # :tA,tB 'N': normal 'T': transpose
    @inbounds for l = it-1:-1:1
        Vt .= vks[l]
        BLAS.gemm!('N','T', 1.0,HV,Vt,0.0,R) #mul!(R,HV,Vt')
        BLAS.gemm!('N','N',-1.0,R,Vt,1.0,HV) #mul!(HV,R,Vt,-1.0,1.0) 
    end
    R .= 0.0
    return nothing
end

# wf:V, twf:HV
function bl_operate_H!(q,wf,twf,pbits,nbits,jocc_p,jocc_n,
                       SPEs,ppinfo,nninfo,tdims,bis,bfs,block_tasks,
                       p_NiNfs,n_NiNfs,Vpn, Mps,delMs,to=nothing) 
    #lblock = length(pbits)
    @inbounds @threads for bi in block_tasks
        if bi==0; continue;end #empty job
        ret = [0,0,0]; ret2 = [0,0,0]
        idim = tdims[bi]
        l_Np = length(pbits[bi])
        l_Nn = length(nbits[bi])
        offset = idim -l_Nn
        #@timeit to "pn1" begin
        @inbounds for (bfidx,bf) in enumerate(bfs[bi])
            bisearch!(delMs,Mps[bf]-Mps[bi],ret) #!! bf = bfs[bi][j]
            Vs = Vpn[ret[1]]
            fdim = tdims[bf]
            l_Nn_f = length(nbits[bf])
            p_NiNf = p_NiNfs[bi][bfidx]
            n_NiNf = n_NiNfs[bi][bfidx]
            off_f = fdim-l_Nn_f
            @inbounds for (nth,V) in enumerate(Vs)
                Npifs = p_NiNf[nth]; Nnifs = n_NiNf[nth]
                @inbounds for Npif in Npifs
                    tMi = offset+ Npif.i .* l_Nn
                    tMf = off_f + Npif.f .* l_Nn_f
                    phase_p = Npif.phase
                    @inbounds for Nnif in Nnifs
                        Mi = tMi + Nnif.i; Mf = tMf + Nnif.f
                        phase_n = Nnif.phase
                        coeff = ifelse(phase_p!=phase_n,-V,V)
                        w_f = @views twf[:,Mi]
                        w_i = @views  wf[:,Mf]
                        @inbounds for b=1:q
                            w_f[b] += coeff .* w_i[b]
                        end
                    end
                end
            end
        end
        #end
        #@timeit to "pn2" begin
        @inbounds for j = 1:length(bis[bi])-1    #### tbf = bi !!!
            tbi = bis[bi][j]
            bisearch!(delMs,Mps[bi]-Mps[tbi],ret)
            bisearch_ord!(bfs[tbi],bi,ret2)
            fdim=tdims[tbi]
            l_Nn_i=length(nbits[tbi])
            p_NiNf = p_NiNfs[tbi][ret2[1]]
            n_NiNf = n_NiNfs[tbi][ret2[1]]
            off_f = fdim - l_Nn_i
            @inbounds for (nth,V) in enumerate(Vpn[ret[1]])
                Npifs = p_NiNf[nth]
                Nnifs = n_NiNf[nth]
                @inbounds for Npif in Npifs
                    tMi = off_f  + Npif.i *l_Nn_i #idim <-> fdim
                    tMf = offset + Npif.f*l_Nn
                    phase_p = Npif.phase
                    @inbounds for Nnif in Nnifs
                        Mi = tMi + Nnif.i; Mf = tMf + Nnif.f
                        phase_n = Nnif.phase
                        coeff = ifelse(phase_p!=phase_n,-V,V)
                        w_f = @views twf[:,Mf]
                        w_i = @views  wf[:,Mi]
                        @inbounds for b=1:q
                            w_f[b] += coeff .* w_i[b]
                        end
                    end
                end
            end
        end
        #end
        #@timeit to "pp/nn" begin
        ### pp/nn interaction
        @inbounds for (Npi,tinfo) in enumerate(ppinfo[bi])
            tMi = offset + Npi*l_Nn
            @inbounds for (jj,tmp) in enumerate(tinfo)
                tMf = offset + l_Nn*tmp.f
                fac = tmp.coef
                @inbounds for nidx = 1:l_Nn
                    Mi = tMi + nidx; Mf = tMf + nidx
                    w_f1 = @views twf[:,Mf]
                    w_i1 = @views  wf[:,Mi]
                    w_f2 = @views twf[:,Mi]                    
                    w_i2 = @views  wf[:,Mf]                        
                    @inbounds for b=1:q
                        w_f1[b] += fac .* w_i1[b]
                        w_f2[b] += fac .* w_i2[b]
                    end
                end
            end
        end
        @inbounds for (Nni,tinfo) in enumerate(nninfo[bi])
            tMi = offset + Nni
            @inbounds for (jj,tmp) in enumerate(tinfo)
                tMf = offset + tmp.f
                fac = tmp.coef
                @inbounds for pidx = 1:l_Np
                    Mi = tMi + pidx*l_Nn; Mf = tMf + pidx*l_Nn
                    w_f1 = @views twf[:,Mf]
                    w_i1 = @views  wf[:,Mi]
                    w_f2 = @views twf[:,Mi]                    
                    w_i2 = @views  wf[:,Mf]                        
                    @inbounds for b=1:q
                        w_f1[b] += fac .* w_i1[b]
                        w_f2[b] += fac .* w_i2[b]
                    end
                end
            end
        end
        #end
        #@timeit to "1b" begin
        ### one-body operator
        @inbounds for pidx = 1:l_Np
            tMi = offset + pidx*l_Nn
            @inbounds for nidx =1:l_Nn
                Mi = tMi + nidx
                coeff = (dot(SPEs[1],jocc_p[bi][pidx])+
                         dot(SPEs[2],jocc_n[bi][nidx]))
                w_f = @views twf[:,Mi]
                w_i = @views  wf[:,Mi]
                @inbounds for b=1:q
                    w_f[b] += coeff .* w_i[b]
                end

            end
        end
        #end
    end
    return nothing
end

function Jcompress(q,vks,Tmat,inow,ls_sub,mdim,R,bnout,U,Mat,Beta_J;use_LAPACK=false)
    lm = q*inow; ls = q*ls_sub
    vals = [0.0]
    vecs = [0.0;0.0]
    if use_LAPACK # used for only for debug (you need "l_diag.F90" by S.Yoshida)
        M = Tmat[1:lm,1:lm]        
        ccall((:diagonalize_double_,"l_diag.so"),Nothing,
              (Ref{Int64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Int64}),
              lm,M,vals,vecs,ls)
    else
        M = @views Tmat[1:lm,1:lm]       
        vals,vecs = eigen(M) # eigen(Symmetric(M,:L))
    end
    # x = maximum([ 1.0-sum(vecs[1:lm,i].^2) for i=1:ls])
    # if x > 1.e-6
    #     println("WARNING: JJ block norm diff. in jj_refine ", x)
    # end
    
    tv = @views vecs[1:ls,1:ls]
    BLAS.gemm!('T','N',1.0,tv,bnout,0.0,R)# mul!(R,tv,bnout)  
    bnout .= R
    
    Tmat .= 0.0; U .= 0.0
    for i=1:inow
        bv = vks[i]
        for b=1:q
            j = q*(i-1)+b
            u = @views U[:,j]
            u .= @views bv[b,:]
        end
    end
    tU = @views U[:,1:lm]
    tv = @views vecs[1:lm,1:ls]
    BLAS.gemm!('T','T',1.0,tv,tU,0.0,vks[1])
    Beta_J .= @views vecs[1:ls,1:ls]   
    
    #println("V1 compressed ",vks[1][1,:])
    #println("V2 compressed ",vks[1][2,:])
    return nothing
end

function bl_ThickRestart(q,vks,uks,R,Tmat,inow,ls_sub,mdim,Vt)
    lm = q*inow
    ls = q*ls_sub
    r = zeros(Float64,q,ls)
    vals,vecs = eigen(@views Tmat[1:lm,1:lm])
    Tmat .= 0.0    
    for k = 1:ls; Tmat[k,k] = vals[k]; end
    tv = @views vecs[lm-q+1:lm,1:ls]    
    mul!(r,R,tv) #BLAS.gemm!('N','N',1.0,R,tv,0.0,r)
    Tmat[ls+1:ls+q,1:ls] .= r
    Tmat[1:ls,ls+1:ls+q] .= r'
    for bi = 1:ls_sub
        uk = uks[bi] #.= 0.0
        uk .= 0.0
        for b=1:q
            k = q*(bi-1) + b
            for j=1:inow
                vk = vks[j]
                for bj = 1:q
                    v = @views vk[bj,:]
                    idx = q*(j-1) +bj
                    fac = vecs[idx,k]
                    @inbounds for m=1:mdim
                        uk[b,m] += fac .* v[m]
                    end
                end
            end
        end
    end
    for i = 1:ls_sub
        vks[i] .= uks[i]
    end
    vks[ls_sub+1] .= vks[inow+1]
    for i= ls_sub+2:length(vks)
        vks[i] .= 0.0
    end
    for i=1:ls_sub+1
        V = vks[i]
        if i>1
            for j = i-1:-1:1
                Vt .= vks[j]
                mul!(R,V,Vt')
                mul!(V,R,Vt,-1.0,1.0) 
            end
            bl_QR!(V',R,mdim,q)
        end
    end
    return nothing
end

function bl_JJ_Lanczos(q,Jvs,Jmat,Vt,R,Beta_J,JTF,Jtol,Jvret,
                       pbits,nbits,tdims,eval_jj,
                       Jidxs,oPP,oNN,oPNu,oPNd,to,bnout,
                       U,Mat;verbose=false)
    mdim = tdims[end]                 
    lnJ = length(Jvs)
    itmin = 1;itmax = lnJ-1; inow = 0
    rescount = 0
    while JTF[1] == false
        for it = itmin:itmax
            inow = it            
            V = Jvs[it]; JV = Jvs[it+1]
            bl_operate_J!(q,V,JV,pbits,nbits,tdims,
                          Jidxs,oPP,oNN,oPNu,oPNd)
            axpy!(eval_jj,V,JV)##JV .+= eval_jj .* V            
            BLAS.gemm!('N','T',1.0,V,JV,0.0,R)
            @inbounds for j=1:q;for i=1:j
                R[j,i] = R[i,j]
            end;end
            @views Jmat[q*it-q+1:q*it,q*it-q+1:q*it] .= R            
            tJmat = @views Jmat[1:it*q,1:it*q]

            jeval = eigvals(tJmat)[1:q]
            if verbose; print_vec("$it jeval=> ",jeval);end
            if jeval[1] - eval_jj < -Jtol
                if verbose;println("neg JJ @$it jeval");end
                JTF[1] = true;break
            end
            if all( abs.(jeval.-eval_jj) .< Jtol)
                if verbose;println("J converged @$it");end
                JTF[1] = true;break
            end
            BLAS.gemm!('N','N',-1.0,R,V,1.0,JV) #mul!(JV,R,V,-1.0,1.0)
            s_vs = @views Jvs[1:it-1]
            bl_ReORTH(q,it,mdim,s_vs,JV,Vt,R)
            bl_QR!(JV',Beta_J,mdim,q)

            tnorm = norm(Diagonal(Beta_J),Inf)
            add_bl_T!(q,it,Jmat,Beta_J)
            if tnorm < Jtol
                println("Jbn norm $tnorm")
                JTF[1]=true;break
            end
        end
        if rescount == 20;println("JJ not converged");return nothing;end
        if JTF[1]==false
            tJmat = @views Jmat[1:inow*q,1:inow*q]
            bl_ThickRestart(q,Jvs,Jvret,Beta_J,Jmat,
                            inow,1,mdim,Vt)
            itmin = 2;rescount += 1
        end
    end
    if inow > itmin
        #@timeit to "Restart"
        Jcompress(q,Jvs,Jmat,inow,1,mdim,R,bnout,U,Mat,Beta_J)
    end
    return nothing
end

function bl_operate_J!(q,Rvec,Jv,
                       pbits,nbits,tdims,
                       Jidxs,
                       oPP,oNN,oPNu,oPNd,
                       beta_J=1.0)
    lblock=length(pbits)
    @inbounds @threads for bi in Jidxs
        if bi==0;continue;end
        idim = tdims[bi]
        lNn = length(nbits[bi])
        lNp = length(pbits[bi])
        opPP = oPP[bi]
        opNN = oNN[bi]
        offset = idim-lNn
        @inbounds for tmp in opPP
            Npi =tmp.Mi; Npf=tmp.Mf; fac=tmp.fac .* beta_J
            tMi = offset + Npi*lNn
            tMf = offset + Npf*lNn
            @inbounds for nidx = 1:lNn
                Mi = tMi+nidx; Mf = tMf+nidx
                w_f1 = @views Jv[:,Mf]
                w_i1 = @views Rvec[:,Mi]
                w_f2 = @views Jv[:,Mi]
                w_i2 = @views Rvec[:,Mf]
                @inbounds for b=1:q
                    w_f1[b] += fac .* w_i1[b]
                    w_f2[b] += fac .* w_i2[b]
                end
            end
        end
        @inbounds for tmp in opNN #nn
            Nni =tmp.Mi; Nnf=tmp.Mf; fac=tmp.fac .* beta_J
            tMi = offset + Nni
            tMf = offset + Nnf
            @inbounds for pidx = 1:lNp
                Mi = tMi+pidx*lNn; Mf = tMf+pidx*lNn
                w_f1 = @views Jv[:,Mf]
                w_i1 = @views Rvec[:,Mi]
                w_f2 = @views Jv[:,Mi]
                w_i2 = @views Rvec[:,Mf]
                @inbounds for b=1:q
                    w_f1[b] += fac .* w_i1[b]
                    w_f2[b] += fac .* w_i2[b]
                end
            end
        end
    end
    @inbounds @threads for bi = 2:lblock
        operator = oPNd[bi]
        bf = bi-1
        l_Nn_i = length(nbits[bi])
        l_Nn_f = length(nbits[bf])
        off_i = tdims[bi] - l_Nn_i
        off_f = tdims[bf] - l_Nn_f
        @inbounds for top in operator
            pj = top.pjump
            nj = top.njump
            fac =top.fac .* beta_J
            @inbounds for tmp_p in pj
                phase_p=tmp_p.phase
                tMi = off_i + tmp_p.i * l_Nn_i
                tMf = off_f + tmp_p.f * l_Nn_f
                @inbounds for tmp_n in nj
                    phase_n=tmp_n.phase
                    Mi = tMi + tmp_n.i
                    Mf = tMf + tmp_n.f
                    coeff = ifelse(phase_p!=phase_n,-fac,fac)
                    w_f = @views Jv[:,Mf]; w_i = @views Rvec[:,Mi]
                    @inbounds for b=1:q
                        w_f[b] += coeff .* w_i[b]
                    end
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
        @inbounds for top in operator
            pj  = top.pjump
            nj  = top.njump
            fac = top.fac .* beta_J
            @inbounds for tmp_p in pj
                phase_p=tmp_p.phase
                tMi = off_i + tmp_p.i * l_Nn_i
                tMf = off_f + tmp_p.f * l_Nn_f
                @inbounds for tmp_n in nj
                    Nni = tmp_n.i; Nnf = tmp_n.f; phase_n=tmp_n.phase
                    Mi = tMi + tmp_n.i
                    Mf = tMf + tmp_n.f
                    coeff = ifelse(phase_p!=phase_n,-fac,fac)
                    w_f = @views Jv[:,Mf]; w_i = @views Rvec[:,Mi]
                    @inbounds for b=1:q
                        w_f[b] += coeff .* w_i[b]
                    end
                end
            end
        end
    end
    return nothing
end
