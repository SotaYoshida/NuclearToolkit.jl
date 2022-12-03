function write_msnt_sps_spe(io,p_sps,n_sps,mstates_p,mstates_n,SPEs)
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
            println(io,@sprintf("%4i",ii),@sprintf("%4i",n), @sprintf("%4i",l), @sprintf("%4i",j), @sprintf("%4i",tz), @sprintf("%5i",mz),@sprintf("%15.6f",SPE[idx]))
        end
    end
    return nothing
end

function trans_snt_msnt(fn,Anum,Mtot,p_sps,n_sps,mstates_p,mstates_n,SPEs,olabels,oTBMEs)
    io = open(replace(fn,".snt"=>".msnt"),"w")
    println(io,"!snt:$fn")
    println(io,"!totalM = $Mtot")
    write_msnt_sps_spe(io,p_sps,n_sps,mstates_p, mstates_n,SPEs)
    write_msnt_tbmes(io,Mtot,p_sps,n_sps,mstates_p,mstates_n,SPEs,olabels,oTBMEs)
    close(io)
    return nothing
end

struct ket_abJ
    a::Int64
    b::Int64
    J::Int64
end

function write_msnt_tbmes(io,Mtot,p_sps,n_sps,mstates_p,mstates_n,SPEs,olabels,oTBMEs)
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
        ctxt = "pp:"
        if pnrank==2; ctxt = "pn:";end
        if pnrank==3; ctxt = "nn:";end
        msps_a = ifelse(pnrank<=2,mstates_p,mstates_n)
        msps_b = ifelse(pnrank>=2,mstates_n,mstates_p)
        ofst_a = ifelse(pnrank<=2,0,lpm)
        ofst_b = ifelse(pnrank>=2,lpm,0)
        for (ia,a) in enumerate(msps_a)
            na,la,ja,tza,mza,idx_a = a
            ia += ofst_a
            for (ib,b) in enumerate(msps_b)
                ib += ofst_b
                if ia > ib; continue;end
                nb,lb,jb,tzb,mzb,idx_b = b
                if Mtot != mza + mzb; continue;end
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
            sps_a = dict_msps2sps[a]
            sps_b = dict_msps2sps[b]
            ja,ma = m_sps[a][3:2:5]
            jb,mb = m_sps[b][3:2:5]
            Jbra = bra.J
            push!(idxs[Jbra+1],idx_bra)
            for (idx_ket,ket) in enumerate(kets[pnrank])
                if idx_bra > idx_ket; continue;end
                c = ket.a
                d = ket.b
                sps_c = dict_msps2sps[c]
                sps_d = dict_msps2sps[d]
                jc,mc = m_sps[c][3:2:5]
                jd,md = m_sps[d][3:2:5]
                Jket = ket.J
                if Jbra != Jket; continue;end
                tbme = dictTBMEs[ [sps_a,sps_b,sps_c,sps_d,Jket] ]
                CG1 = clebschgordan(Float64,ja//2,ma//2,jb//2,mb//2, Jket, Mtot)
                CG2 = clebschgordan(Float64,jc//2,mc//2,jd//2,md//2, Jket, Mtot)
                tbme_M = tbme * sqrt( (1.0+deltaf(sps_a,sps_b)) *(1.0+deltaf(sps_c,sps_d)) ) * CG1 * CG2
                println(io, ctxt, @sprintf("%5i",a),@sprintf("%5i",b),@sprintf("%5i",c),@sprintf("%5i",d), @sprintf("%5i",Jket), @sprintf("%15.6f",tbme_M))                
                mat[idx_bra,idx_ket] = mat[idx_ket,idx_bra] =  tbme_M
            end
        end
        if pnrank == 2
            for J = 0:10
                idx = idxs[J+1]
                if length(idx) == 0; continue; end
                submat = mat[idx,idx]
                evals,evecs = eigen(submat)
                println("J $J ")
                for nn = 1:length(evals)
                    if abs(evals[nn]) < 1.e-12; continue;end
                    println(@sprintf("%15.6f", evals[nn]) )
                #     vec = evecs[:,nn]
                #     E1b = 0.0
                #     for tmp in idx
                #         tket =  kets[pnrank][tmp]
                #         a = dict_msps2sps[tket.a]
                #         b = dict_msps2sps[tket.b]                            
                #         E1b += (allSPEs[a] + allSPEs[b])
                #     end
                #     E = E1b + evals[nn]
                #     println("E1b ", @sprintf("%15.6f",E1b), " E ",@sprintf("%15.6f",E))
                end
            end
        end
    end
    return nothing
end

function main_trans_msnt(fn,target_nuc,target_Js=[]) 
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
    trans_snt_msnt(fn,Anum,Mtot,p_sps,n_sps,mstates_p,mstates_n,SPEs,olabels,oTBMEs)
end