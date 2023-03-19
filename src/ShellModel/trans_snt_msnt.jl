function make_dict_convention(similar_to_prevwork::Bool)
    dict_sps = Dict{Int64,Int64}()    
    if similar_to_prevwork 
        p_sps,n_sps
    else
        dict_sps[i]
    end
    return dict_sps
end

function write_msnt_sps_spe(io,p_sps,n_sps,mstates_p,mstates_n,SPEs,similar_to_prevwork)
    println(io, "!single particle space and SPEs")
    println(io, "!num   n   l   j  tz   mz       SPE(MeV)")
    ii = 0
    for pn = 1:2
        sps  = ifelse(pn==1,p_sps,n_sps)
        msps = ifelse(pn==1,mstates_p,mstates_n)
        SPE = SPEs[pn]
        for i = 1:length(msps)
            n,l,j,tz,mz,idx = msps[i]
            ii += 1
            qbit_idx = ifelse(similar_to_prevwork,ii-1,ii)
            println(io,@sprintf("%4i",qbit_idx),@sprintf("%4i",n), @sprintf("%4i",l), @sprintf("%4i",j), @sprintf("%4i",tz), @sprintf("%5i",mz),@sprintf("%15.6f",SPE[idx]))
        end
    end
    return nothing
end

function trans_snt_msnt(fn,Anum,Mtot,p_sps,n_sps,mstates_p,mstates_n,SPEs,olabels,oTBMEs,similar_to_prevwork)
    io = open(replace(fn,".snt"=>".msnt"),"w")
    println(io,"!snt:$(fn)\n!")
    write_msnt_sps_spe(io,p_sps,n_sps,mstates_p,mstates_n,SPEs,similar_to_prevwork)
    write_msnt_tbmes(io,Mtot,p_sps,n_sps,mstates_p,mstates_n,SPEs,olabels,oTBMEs,similar_to_prevwork)
    close(io)
    return nothing
end

struct ket_abJ
    a::Int64
    b::Int64
    J::Int64
end

function write_msnt_tbmes(io,Mtot,p_sps,n_sps,mstates_p,mstates_n,SPEs,olabels,oTBMEs,similar_to_prevwork)
    lp = length(p_sps)
    lpm = length(mstates_p)
    m_sps = vcat(mstates_p,mstates_n)
    allSPEs = vcat(SPEs[1],SPEs[2])
    dict_msps2sps = Dict{Int64,Int64}()
    for (idx,tmp) in enumerate(mstates_p)
        dict_msps2sps[idx] = tmp[end]
    end
    for (idx,tmp) in enumerate(mstates_n)
        dict_msps2sps[idx+lpm] = tmp[end] + lp
    end

    dictTBMEs = Dict{Vector{Int64},Float64}()
    for (i,tkey) in enumerate(olabels)
        a,b,c,d,J,oidx = tkey
        dictTBMEs[[a,b,c,d,J]] = oTBMEs[i]
        dictTBMEs[[c,d,a,b,J]] = oTBMEs[i]
    end

    kets = [ ket_abJ[ ] for pnrank=1:3]
    for pnrank = 1:3
        ctxt = "!Vpp:"
        if pnrank==2; ctxt = "!Vpn:";end
        if pnrank==3; ctxt = "!Vnn:";end
        println(io,ctxt)
        msps_a = ifelse(pnrank<=2,mstates_p,mstates_n)
        msps_b = ifelse(pnrank>=2,mstates_n,mstates_p)
        ofst_a = ifelse(pnrank<=2,0,lpm)
        ofst_b = ifelse(pnrank>=2,lpm,0)
        for (ia,a) in enumerate(msps_a)
            na,la,ja,tza,mza,idx_a = a
            ia += ofst_a
            for (ib,b) in enumerate(msps_b)
                ib += ofst_b
                if ia >= ib; continue;end
                nb,lb,jb,tzb,mzb,idx_b = b
                Tz = tza + tzb
                for J = abs(div(ja-jb,2)):div(ja+jb,2)
                    if pnrank % 2 == 1 && J % 2 == 1; continue; end
                    push!(kets[pnrank],ket_abJ(ia,ib,J))
                end
            end
        end
        dim = length(kets[pnrank])
        mat = zeros(Float64,dim,dim)
        idxs = [ Int64[ ] for J =0:10]
        for (idx_bra,bra) in enumerate(kets[pnrank])
            a = bra.a
            b = bra.b            
            qbit_a = a + ifelse(similar_to_prevwork,-1,0)
            qbit_b = b + ifelse(similar_to_prevwork,-1,0)    
            sps_a = dict_msps2sps[a]
            sps_b = dict_msps2sps[b]
            ja,ma = m_sps[a][3:2:5]
            jb,mb = m_sps[b][3:2:5]
            Jbra = bra.J
            if ma + mb == Mtot; push!(idxs[Jbra+1],idx_bra);end
            for (idx_ket,ket) in enumerate(kets[pnrank])
                if idx_bra > idx_ket; continue;end
                c = ket.a
                d = ket.b
                qbit_c = c + ifelse(similar_to_prevwork,-1,0)
                qbit_d = d + ifelse(similar_to_prevwork,-1,0)    
                sps_c = dict_msps2sps[c]
                sps_d = dict_msps2sps[d]
                jc,mc = m_sps[c][3:2:5]
                jd,md = m_sps[d][3:2:5]
                Jket = ket.J
                if ma + mb - mc - md != 0; continue; end                
                if Jbra != Jket; continue;end
                if abs(ma+mb) > 2*Jket;continue;end
                tbme = dictTBMEs[ [sps_a,sps_b,sps_c,sps_d,Jket] ]
                CG1 = clebschgordan(Float64,ja//2,ma//2,jb//2,mb//2, Jket, div(ma+mb,2))
                CG2 = clebschgordan(Float64,jc//2,mc//2,jd//2,md//2, Jket, div(mc+md,2))
                tbme_M = tbme * sqrt( (1.0+deltaf(sps_a,sps_b)) *(1.0+deltaf(sps_c,sps_d)) ) * CG1 * CG2
                if abs(tbme_M) < 1.e-9; continue;end
                SPE_a = allSPEs[sps_a]
                SPE_b = allSPEs[sps_b]
                e1b = ifelse(a==c && b==d,SPE_a+SPE_b,0.0)
                mat[idx_bra,idx_ket] = mat[idx_ket,idx_bra] = e1b + tbme_M
                println(io, @sprintf("%5i",qbit_a),@sprintf("%5i",qbit_b),
                        @sprintf("%5i",qbit_c),@sprintf("%5i",qbit_d), @sprintf("%5i",Jket), @sprintf("%15.6f",tbme_M))               
            end
        end
        if true && pnrank == 2 # Li6
            #for J = 0:10
            for J = 2:3
                idx = idxs[J+1]
                if length(idx) == 0; continue; end
                submat = mat[idx,idx]
                evals,evecs = eigen(submat)
                println("J $J idx $idx")
                for nn = 1:length(idx)
                    ket = kets[pnrank][idx[nn]]
                    M = m_sps[ket.a][5] +  m_sps[ket.b][5] 
                    print("|$(ket.a) $(ket.b);M=$M> ")
                end
                print("\n")
                for nn = 1:length(idx)
                    print_vec("", submat[nn,:])
                end 
                for nn = 1:length(evals)
                    if abs(evals[nn]) < 1.e-12; continue;end
                    vec = evecs[:,nn]
                    print_vec(" E "*@sprintf("%15.6f",evals[nn]), vec)
                end
            end
        end
    end
    return nothing
end

"""
if `similar_to_prevwork == true`, indices for qbit starts from 0 and protons have tz = 1, which is opposite to the convensiton used in many-body methods in NuclearToolkit.jl.
"""
function main_trans_msnt(fn,target_nuc,target_Js=[];similar_to_prevwork=false) 
    Anum = parse(Int64, match(reg,target_nuc).match)
    Mtot = 0;if Anum % 2 != 0; Mtot = 1;end
    if Anum % 2 != Mtot % 2
        @error "invalid target_Js = $target_Js"
    end
    if length(target_Js) > 0
        Mtot = minimum(target_Js)
        tJ=target_Js[1]
        eval_jj = 0.5*tJ*(tJ/2+1)
    end
    lp,ln,cp,cn,massop,Aref,p,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsmsnt(fn,Anum)

    target_el = replace.(target_nuc, string(Anum)=>"")
    Z,N,vp,vn = getZNA(target_el,Anum,cp,cn)
    mstates_p,mstates_n,mz_p,mz_n = def_mstates(p_sps,n_sps)
    trans_snt_msnt(fn,Anum,Mtot,p_sps,n_sps,mstates_p,mstates_n,SPEs,olabels,oTBMEs,similar_to_prevwork)
end

function svd_li6(Rvecs)
    println("svd_li6 function is called. This is only used for Li6 on the p-shell space.")
    dict_sps = Dict(1=>[3,3],2=>[1,1],3=>[1,4],4=>[4,1],5=>[4,4],
                    6=>[2,2],7=>[2,5],8=>[5,2],9=>[5,5],10=>[6,6])
    nsps = 6
    Mat = zeros(Float64,nsps,nsps)
    for n = 1:length(Rvecs)
        Rvec = Rvecs[n]
        println("n $n")
        for i = 1:length(Rvec)
            idxs = dict_sps[i]
            Mat[idxs[1],idxs[2]] = Rvec[i]
        end
        for j = 1:nsps
            print_vec("",Mat[j,:];long=true)
        end
        SVD = LinearAlgebra.svd(Mat)
        U = SVD.U; Sig = Diagonal(SVD.S); Vt = SVD.Vt; V = Vt'
        println("U")
        for j = 1:nsps
            print_vec("",U[j,:];long=true)
        end 
        println("Sig")
        for j = 1:nsps
            print_vec("",Sig[j,:];long=true)
        end 
        println("V")
        for j = 1:nsps
            print_vec("",V[j,:];long=true)
        end 
        println("")
        # println("Vt")
        # for j = 1:nsps
        #     print_vec("",Vt[j,:])
        # end 
        Mat .= 0.0
    end
end