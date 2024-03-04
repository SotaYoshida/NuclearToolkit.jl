"""
The main API for VMC (w/ pairwise ansatz).
"""
function main_vmc(psi_exact,eigenvals,tdims,msps_p,msps_n,
                  pbits,nbits,jocc_p,jocc_n,SPEs,pp_2bjump,nn_2bjump,bis,bfs,block_tasks,p_NiNfs,n_NiNfs,Mps,delMs,Vpn,
                  Jidxs,oPP,oNN,oPNu,oPNd)
    state_idx = 1
    psi_gs = psi_exact[state_idx]
    mdim = tdims[end]
    @assert mdim <= 300 "vmc function is now valid for mdim <= 300"

    # Make block w.f. for the ground state
    psi_block = [ zeros(Float64,length(pbits[bi]),length(nbits[bi])) for bi = 1:length(pbits)]
    make_block_wf!(psi_gs,psi_block,pbits,nbits)
    #svd_block_wf(psi_block)

    # Flatten Hamiltonian operator to (mdim,mdim) matrix => Hflat
    Hflat = calc_Hflat(tdims,pbits,nbits,jocc_p,jocc_n,SPEs,pp_2bjump,nn_2bjump,bis,bfs,block_tasks,p_NiNfs,n_NiNfs,Mps,delMs,Vpn)
    print_vec("eigvals(Hflat)",eigvals(Hflat))

    ### low-rank approximation of "many-body wavefunction"

    pdim = sum( [size(psi_block[i])[1] for i = 1:length(pbits)] )
    ndim = sum( [size(psi_block[i])[2] for i = 1:length(pbits)] )

    psi_config = zeros(Float64,pdim,ndim)
    p_ofst = n_ofst = 0 
    flat_idx = 1
    for bi = 1:length(pbits) 
        for pidx = 1:length(pbits[bi])
            for nidx = 1:length(nbits[bi])
                psi_config[p_ofst+pidx,n_ofst+nidx] = psi_gs[flat_idx]
                flat_idx += 1
            end
        end
        p_ofst += length(pbits[bi])
        n_ofst += length(nbits[bi])
    end

    # construct product state 
    dict_idxs = prep_product_state(pbits,nbits,msps_p,msps_n) 
 
    U,S,V = svd(psi_config)
    print_vec("S:\n",S)

    pbits_flat = vcat(pbits...)
    nbits_flat = vcat(nbits...)
    println("n_msps $msps_n")
    println("neutron configs. $nbits_flat")


    for trank in 1:rank(psi_config)    
        idxs_rank = collect(1:trank)

        # # # for 8Be remove index from Array
        # idxs_rank = Int64[ ]
        # for i = 1:size(V)[2]
        #     #println("abs ", abs.(V'[i,[3,4,5,11,12,13]]), " sum(abs)", sum(abs.(V'[i,[3,4,5,11,12,13]])))
        #     if sum(abs.(V'[i,[3,4,5,11,12,13]])) > 1.e-8; continue;end
        #     push!(idxs_rank,i)
        # end
        # idxs_rank = idxs_rank[1:min(length(idxs_rank),trank)]
        # println("idxs_rank $idxs_rank")

        subC = U[:,idxs_rank]*Diagonal(S[idxs_rank])*V'[idxs_rank,:]


        subpsi = subC  
        subpsi ./= sqrt(sum(subpsi.^2))

        #print_mat("U",U)
        #print_mat("V^T",V')

        vec_trial = zeros(Float64,mdim)
        p_ofst = n_ofst = 0; flat_idx = 1
        for bi = 1:length(pbits)
            for pidx = 1:length(pbits[bi])
                for nidx = 1:length(nbits[bi])
                    vec_trial[flat_idx] = subpsi[p_ofst+pidx,n_ofst+nidx]
                    flat_idx += 1
                end
            end
            p_ofst += length(pbits[bi])
            n_ofst += length(nbits[bi])
            # println("block $bi\n",psi_block[bi])
            # println("pbits $(pbits[bi])")
            #println("nbits $(nbits[bi])")
        end
        Etri = dot(vec_trial,Hflat*vec_trial)
        num_non0 = sum(abs.(vec_trial) .> 1.e-12)
        println("trank ",@sprintf("%3i",trank)," <H>_trial = ", @sprintf("%10.5f", Etri), "\t # nonzero $num_non0 \t rel. error (%) ", @sprintf("%9.3e", 100*abs((Etri-eigenvals[state_idx])/eigenvals[state_idx])))

        if trank == rank(psi_config)
            print_mat("V^T",V')
            print_mat("subpsi",subpsi)
            for i = 1:size(V)[2]
                println(V'[i,:])
            end
        end
       
    end        
    return nothing
    
    psi_svd = zeros(Float64,pdim,ndim)
    psi_svd .= subC  
    psi_svd ./= sqrt(sum(psi_svd.^2))

    psi_flat_new = zeros(Float64,mdim)
    p_ofst = n_ofst = 0 
    flat_idx = 1
    for bi = 1:length(pbits) 
        for pidx = 1:length(pbits[bi])
            for nidx = 1:length(nbits[bi])
                psi_flat_new[flat_idx] = psi_config[p_ofst+pidx,n_ofst+nidx] 
                flat_idx += 1
            end
        end
        p_ofst += length(pbits[bi])
        n_ofst += length(nbits[bi])
    end
    Egs = dot(psi_flat_new,Hflat*psi_flat_new)
    println("trank = $trank, Egs = $Egs error ",abs(Egs-eigenvals[1]))
    for i = 1:pdim
        print_vec("",psi_svd[i,:])
    end
   
    # SVD of Hflat
    #M0idxs = get_M0idxs(pbits,nbits,msps_p,msps_n)
    #svd_mat(Hflat,M0idxs)
    
    # Check the energies by direct diagonalization of Hflat
    num_ev = length(eigenvals)
    vals = real.(eigsolve(Hflat,num_ev,:SR)[1])[1:num_ev]
    @assert norm(vals-eigenvals, Inf) < 1.e-6 "Hflat not giving the same eigenvalues as H"
   
    vtmp = zeros(Float64,mdim)
    operate_J!(psi_gs,vtmp,pbits,nbits,tdims,Jidxs,oPP,oNN,oPNu,oPNd)
    Jexpec = dot(psi_gs,vtmp)
    totalJ = J_from_JJ1(Jexpec)
    println("<Hflat>_0 = ", dot(psi_gs,Hflat * psi_gs), " <J>_0 = $totalJ")

    # Specifying pairwise truncations for VMC
    # mask_idxs, Midxs = truncated_Midxs_pnpair(pbits,nbits,msps_p,msps_n)
    # println("pn-pair truncated dim. => $(length(Midxs))")
    mask_idxs, Midxs = truncated_Midxs_absM(pbits,nbits,msps_p,msps_n)
    println("pn-pair + jpartner truncated dim. => $(length(Midxs))")
    # mask_idxs, Midxs = truncated_Midxs_MpMn0(pbits,nbits,msps_p,msps_n)
    # println("Mp=Mn=0 configs. dim. => $(length(Midxs))")

    # Make trial wavefunction
    psi_trial = copy(psi_gs); psi_trial[mask_idxs] .= 0.0

    subH = @view Hflat[Midxs,Midxs]
    svals, svecs = eigsolve(subH,1,:SR)
    psi_trial[Midxs] .= svecs[1]
    psi2 = dot(psi_trial,psi_trial)
    #println("psi_trial $psi_trial psi2 $psi2")
    psi_trial ./= (psi2)

    # check total angular momentum J 
    vtmp .= 0.0
    operate_J!(psi_trial,vtmp,pbits,nbits,tdims,Jidxs,oPP,oNN,oPNu,oPNd)
    Jexpec = dot(psi_trial,vtmp)
    totalJ = J_from_JJ1(Jexpec)
    println("<H>_tri = ", dot(psi_trial,Hflat * psi_trial), " <J>_tri = $totalJ Jexpec $Jexpec \n")

    return nothing
end

"""
For each proton or neutron configuration have 'jz-mirror' configurations. They should have the same amplitude.
"""
function prep_product_state(pbits,nbits,msps_p,msps_n)
    num_pbit = length(msps_p)
    num_nbit = length(msps_n)
    pbits_flat = vcat(pbits...)
    nbits_flat = vcat(nbits...)

    partner_dict_p = Dict{Int64,Int64}()
    for (idx,tmp) in enumerate(msps_p)
        n,l,j,tz,m = tmp
        #println("idx $idx nljtzm $n $l $j $tz $m")
        for (idx_,tmp_) in enumerate(msps_p)
            n_,l_,j_,tz_,m_ = tmp_
            if n == n_ && l == l_ && j == j_ && tz == tz_ && m == -m_
                #println("partner => idx_ $idx_ nljtzm $n_ $l_ $j_ $tz_ $m_")
                partner_dict_p[idx] = idx_
            end
        end
    end

    dict_idxs = Dict{Int64,Int64}()
    idx = 1
    for (conf_idx_1,pbit_1) in enumerate(pbits_flat)
        occ_vec_p = get_bitarry_from_int(pbit_1, num_pbit)
        hole_idxs = findall(occ_vec_p .==1)
        partner_idxs = [partner_dict_p[hole_idx] for hole_idx in hole_idxs]
        for (conf_idx_2,pbit_2) in enumerate(pbits_flat)
            occ_vec_p_2 = get_bitarry_from_int(pbit_2, num_pbit)
            hole_idxs_2 = findall(occ_vec_p_2 .==1)
            if Set(partner_idxs) != Set(hole_idxs_2); continue; end
            println("conf_idx $conf_idx_1 $pbit_1 partner => $conf_idx_2 $pbit_2")
            dict_idxs[conf_idx_1] = conf_idx_2
        end            
    end
    return dict_idxs
end

function print_mat(text,mat)
    println("$text:")
    for i = 1:size(mat)[1]
        print_vec("",mat[i,:])
    end
    return nothing
end

function Midxs_random(mdim,num_use)
    tmp = randperm(mdim)
    Midxs = tmp[1:num_use]
    mask_idxs = tmp[num_use+1:end]
    return mask_idxs,Midxs
end

"""
Function to return M-scheme indices for truncated configurations
"""
function truncated_Midxs_pnpair(pbits,nbits,msps_p,msps_n)
    Midxs = Int64[ ]
    mask_idxs = Int64[ ]
    idx = 0 
    num_pbit = length(msps_p)
    num_nbit = length(msps_n)
    for block_i = 1:length(pbits)
        l_Np = length(pbits[block_i])
        l_Nn = length(nbits[block_i])
        for pidx = 1:l_Np
            pbit = pbits[block_i][pidx]
            occ_vec_p = get_bitarry_from_int(pbit, num_pbit)
            for nidx = 1:l_Nn
                idx += 1
                nbit = nbits[block_i][nidx]
                occ_vec_n = get_bitarry_from_int(nbit, num_nbit)
                valid = true
                for (idx_nbit,occ_n) in enumerate(occ_vec_n)
                    n_n,l_n,j_n,tz_n,m_n = msps_n[idx_nbit]
                    if occ_n == 1 
                        for (idx_pbit,occ_p) in enumerate(occ_vec_p)
                            n_p,l_p,j_p,tz_p,m_p = msps_p[idx_pbit]
                            if n_p != n_n || l_p != l_n || j_p != j_n || tz_p == tz_n || m_p != -m_n
                                continue
                            end
                            valid *= ifelse(occ_p==occ_n,true,false)
                        end
                    end
                end
                if valid
                    push!(Midxs,idx)
                else
                    push!(mask_idxs,idx)
                end
            end
        end
    end
    return mask_idxs,Midxs
end

function truncated_Midxs_absM(pbits,nbits,msps_p,msps_n;jpartner=false)
    Midxs = Int64[ ]
    mask_idxs = Int64[ ]
    idx = 0 
    num_pbit = length(msps_p)
    num_nbit = length(msps_n)
    for block_i = 1:length(pbits)
        l_Np = length(pbits[block_i])
        l_Nn = length(nbits[block_i])
        for pidx = 1:l_Np
            pbit = pbits[block_i][pidx]
            occ_vec_p = get_bitarry_from_int(pbit, num_pbit)
            for nidx = 1:l_Nn
                idx += 1
                nbit = nbits[block_i][nidx]
                occ_vec_n = get_bitarry_from_int(nbit, num_nbit)
                absMn = absMp = 0
                for (idx_nbit,occ_n) in enumerate(occ_vec_n)
                    n_n,l_n,j_n,tz_n,m_n = msps_n[idx_nbit]
                    if occ_n == 1; absMn += abs(m_n);end
                end
                for (idx_pbit,occ_p) in enumerate(occ_vec_p)
                    n_p,l_p,j_p,tz_p,m_p = msps_p[idx_pbit]
                    if occ_p == 1; absMp += abs(m_p);end
                end
                if absMn==absMp
                    push!(Midxs,idx)
                else
                    push!(mask_idxs,idx)
                end
            end
        end
    end
    return mask_idxs,Midxs
end

function truncated_Midxs_MpMn0(pbits,nbits,msps_p,msps_n)
    Midxs = Int64[ ]
    mask_idxs = Int64[ ]
    idx = 0 
    num_pbit = length(msps_p)
    num_nbit = length(msps_n)
    for block_i = 1:length(pbits)
        l_Np = length(pbits[block_i])
        l_Nn = length(nbits[block_i])
        for pidx = 1:l_Np
            pbit = pbits[block_i][pidx]
            occ_vec_p = get_bitarry_from_int(pbit, num_pbit)
            for nidx = 1:l_Nn
                idx += 1
                nbit = nbits[block_i][nidx]
                occ_vec_n = get_bitarry_from_int(nbit, num_nbit)
                Mn = Mp = 0
                for (idx_nbit,occ_n) in enumerate(occ_vec_n)
                    n_n,l_n,j_n,tz_n,m_n = msps_n[idx_nbit]
                    if occ_n == 1; Mn += m_n;end
                end
                for (idx_pbit,occ_p) in enumerate(occ_vec_p)
                    n_p,l_p,j_p,tz_p,m_p = msps_p[idx_pbit]
                    if occ_p == 1; Mp += m_p;end
                end
                if  Mp == Mn == 0
                    push!(Midxs,idx)
                else
                    push!(mask_idxs,idx)
                end
            end
        end
    end
    return mask_idxs,Midxs
end

function calc_Hflat(tdims,pbits,nbits,jocc_p,jocc_n,SPEs,pp_2bjump,nn_2bjump,bis,bfs,block_tasks,p_NiNfs,n_NiNfs,Mps,delMs,Vpn)
    mdim = tdims[end]
    Hflat = zeros(Float64,mdim,mdim)
    for bi in block_tasks
        if bi==0; continue;end #empty job
        ret = [0,0,0]; ret2 = [0,0,0]
        idim = tdims[bi]
        l_Np = length(pbits[bi])
        l_Nn = length(nbits[bi])
        offset = idim -l_Nn

        # 1. one-body term 
        for pidx = 1:l_Np
            tMi = offset + pidx*l_Nn
            for nidx =1:l_Nn
                Mi = tMi + nidx
                Hflat[Mi,Mi] += (dot(SPEs[1],jocc_p[bi][pidx])+dot(SPEs[2],jocc_n[bi][nidx]))
            end
        end

        # 2. two-body term (pp/nn)
        for (pidx,tinfo) in enumerate(pp_2bjump[bi])
            tMi = offset + pidx*l_Nn
            for (jj,tmp) in enumerate(tinfo)
                tMf = offset + l_Nn*tmp.f
                fac = tmp.coef
                for nidx = 1:l_Nn
                    Mi = tMi + nidx; Mf = tMf + nidx
                    @assert Mi <= Mf "Mi > Mf happened in pp2b"
                    Hflat[Mi,Mf] += fac
                    Hflat[Mf,Mi] += fac
                end
            end
        end
        for (Nni,tinfo) in enumerate(nn_2bjump[bi])
            tMi = offset + Nni
            for (jj,tmp) in enumerate(tinfo)
                tMf = offset + tmp.f
                fac = tmp.coef
                for pidx = 1:l_Np
                    Mi = tMi + pidx*l_Nn; Mf = tMf + pidx*l_Nn
                    @assert Mi <= Mf "Mi > Mf happened in nn2b"
                    Hflat[Mi,Mf] += fac
                    Hflat[Mf,Mi] += fac
                end
            end
        end
        
        # 3. two-body term (pn)
        for (bfidx,bf) in enumerate(bfs[bi])
            bisearch!(delMs,Mps[bf]-Mps[bi],ret) #!! bf = bfs[bi][j]
            Vs = Vpn[ret[1]]
            fdim = tdims[bf]
            l_Nn_f = length(nbits[bf])
            p_NiNf = p_NiNfs[bi][bfidx]
            n_NiNf = n_NiNfs[bi][bfidx]
            off_f = fdim-l_Nn_f
            for (nth,V) in enumerate(Vs)
                Npifs = p_NiNf[nth]; Nnifs = n_NiNf[nth]
                for Npif in Npifs
                    tMi = offset+ Npif.i * l_Nn
                    tMf = off_f + Npif.f * l_Nn_f
                    phase_p = Npif.phase
                    for Nnif in Nnifs
                        Mi = tMi + Nnif.i; Mf = tMf + Nnif.f
                        phase_n = Nnif.phase
                        Hflat[Mi,Mf] += ifelse(phase_p!=phase_n,-V,V)                       
                    end
                end
            end
        end
        
        for j = 1:length(bis[bi])-1    #### tbf = bi !!!
            tbi = bis[bi][j]
            bisearch!(delMs,Mps[bi]-Mps[tbi],ret)
            bisearch_ord!(bfs[tbi],bi,ret2)
            fdim=tdims[tbi]
            l_Nn_i=length(nbits[tbi])
            p_NiNf = p_NiNfs[tbi][ret2[1]]
            n_NiNf = n_NiNfs[tbi][ret2[1]]
            off_f = fdim - l_Nn_i
            for (nth,V) in enumerate(Vpn[ret[1]])
                Npifs = p_NiNf[nth]
                Nnifs = n_NiNf[nth]
                for Npif in Npifs
                    tMi = off_f  + Npif.i *l_Nn_i #idim <-> fdim
                    tMf = offset + Npif.f*l_Nn
                    phase_p = Npif.phase
                    for Nnif in Nnifs
                        Mi = tMi + Nnif.i; Mf = tMf + Nnif.f
                        phase_n = Nnif.phase
                        @assert Mi <= Mf "Mi > Mf happened in pn2b-2"
                        Hflat[Mf,Mi] += ifelse(phase_p!=phase_n,-V,V)
                    end
                end
            end
        end
    end
    tnorm = norm(Hflat-Hflat',Inf)
    @assert tnorm < 1e-10 "Hflat must be symmetric!! norm(,Inf) = $tnorm"
    return Hflat
end

function svd_block_wf(wf)
    for block = 1:length(wf)
        block_wf = wf[block]
        println("block $block")
        for i = 1:size(block_wf)[1]
            print_vec("",block_wf[i,:])
        end
        U,S,V = svd(block_wf)
        println("U:")
        for i = 1:size(U)[1]
            print_vec("",U[i,:])
        end
        println("S:")
        print_vec("",S)
        println("V:")
        for i = 1:size(V)[1]
            print_vec("",V[i,:])
        end
        println("|psi'>=V^T|psi>")
        nwf = V'*block_wf
        for i = 1:size(nwf)[1]
            print_vec("",nwf[i,:])
        end
        println("norm(SVD-H,2)",norm(U*Diagonal(S)*V'-block_wf,2))
    end
    return nothing
end

function get_M0idxs(pbits,nbits,msps_p,msps_n)
    M0idxs = [Int64[ ],Int64[ ]]
    nblock = length(pbits)
    idx = 0
    for bi = 1:nblock
        t_pbits = pbits[bi]
        t_nbits = nbits[bi]
        for pidx = 1:length(t_pbits)
            pbit = t_pbits[pidx]
            Mp = 0 
            for (idx_pbit,bit_p) in enumerate(pbit)
                n_p, l_p, j_p, tz_p, m_p = msps_p[idx_pbit]
                if bit_p == 1
                    Mp += m_p
                end
            end
            for nidx = 1:length(t_nbits)
                nbit = t_nbits[nidx]
                Mn = 0
                for (idx_nbit,bit_n) in enumerate(nbit)
                    n_n, l_n, j_n, tz_n, m_n = msps_n[idx_nbit]
                    idx += 1
                    if bit_n == 1
                        Mn += m_n
                    end                    
                end
                if Mp + Mn == 0
                    push!(M0idxs[1],pidx)
                    push!(M0idxs[2],nidx)
                end
            end
        end
    end
    println("idx check $idx")
    return M0idxs
end

function svd_mat(Hflat,M0idxs)
    U,S,V = svd(Hflat)
    print_vec("S:",S)

    # for trank = 1:length(S)
    #     subH = U[:,1:trank]*Diagonal(S[1:trank])*V[:,1:trank]'
    #     vals,vecs = eigsolve(subH,1,:SR)
    #     Egs = real.(vals)[1]
    #     diffnorm = norm(Hflat-subH,2)
    #     diffnorm_M0 = norm(Hflat[M0idxs[1],M0idxs[2]]-subH[M0idxs[1],M0idxs[2]],2)
    #     println("trank = $trank, Egs = $Egs diffnorm $diffnorm diffnorm(M0) $diffnorm_M0 ")
    # end
    return nothing
end

function make_block_wf!(psi_flat,psi_block,pbits,nbits)
    flat_idx = 1
    for block = 1:length(pbits)
        pbits_in_b = pbits[block]
        nbits_in_b = nbits[block]
        target = psi_block[block]
        ln_nbits = length(nbits_in_b)
        for pidx = 1:length(pbits_in_b)
            for nidx = 1:length(nbits_in_b)
                target[pidx,nidx] = psi_flat[flat_idx]
                flat_idx += 1
            end
        end
    end   

    mdim = length(psi_flat)
    Psi_block = zeros(Float64,mdim,mdim)
    p_ofst = n_ofst = 0 
    flat_idx = 1
    for block = 1:length(pbits)
        pbits_in_b = pbits[block]
        nbits_in_b = nbits[block] 
        #println("block $block pbits $(length(pbits_in_b)) nbits $(length(nbits_in_b))")      
        for pidx = 1:length(pbits_in_b)
            for nidx = 1:length(nbits_in_b)
                #println("p_ofst $p_ofst + pidx $pidx,  n_ofst $n_ofst + nidx $nidx")
                Psi_block[p_ofst+pidx,n_ofst+nidx] = psi_flat[flat_idx]
                flat_idx += 1
            end
        end
        p_ofst += length(pbits_in_b)
        n_ofst += length(nbits_in_b)
    end   

    # U,S,V = svd(Psi_block)
    # Snon0 = Float64[ ]
    # tol = 1.e-8
    # for tmp in S
    #     if tmp > tol
    #         push!(Snon0,tmp)
    #     end
    # end
    # print_vec("singular values (#=$(length(Snon0))):",Snon0;ine=true)
    # println("size $(size(Psi_block))")
    # for trank = 1:length(Snon0)
    #     subH = U[:,1:trank]*Diagonal(S[1:trank])*V[:,1:trank]'
    #     Egs = dot(psi_flat,subH*psi_flat)
    #     println("trank = $trank, Egs = $Egs")
    # end

    return nothing
end
