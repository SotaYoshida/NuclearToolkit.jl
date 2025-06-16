function Calc_3NF_in_Jacobi_coordinate(params3N,LECs,to)
    Rnls_r,Rnls_p = prep_Rnls(params3N)
    @timeit to "prep. ch&Mat" begin
        Jacobi_idxs = prep_Jacobi_idxs(params3N)        
        chan3b = estimate_required(params3N)
        vec_freg_nl = ret_fregvec_nonlocal(params3N;npow=2)
        ## Calculate 3NF 
        pqdim = length(params3N.meshpoints.pqs)
        Mats_pq = [ zeros(Float64,pqdim,pqdim) for i = 1:nthreads()]
        @timeit to "Contact" Jacobi_3NF_Contact(LECs,chan3b,Jacobi_idxs,params3N,Mats_pq,vec_freg_nl)        
    end
end

function prep_Jacobi_idxs(params)
    pqs = params.meshpoints.pqs
    dim = size(pqs)[1]
    pq_idxs = Vector{Int64}[]
    for idx_ket = 1:dim
        for idx_bra = 1:idx_ket
            push!(pq_idxs,[idx_bra,idx_ket])
        end
    end
    return pq_idxs
end

"""

calculate |pq;alpha>
This reproduces Table.3 in K.Hebeler Physics Reports 890 (2021) 1–116.
"""
function estimate_required(params;verbose=true)
    J12max = params.J12max
    Np = length(params.meshpoints.ps)
    Nq = length(params.meshpoints.qs)
    j3max = params.j3max
    chs = channel3b_Jj[]
    ich = 0 
    for P123 = 1:-2:-1  
        for dT123 = 1:2:3
            for dJ123 = 1:2:j3max
                @assert abs(P123)==1 "P123 must be +1 or -1"
                @assert (dT123 == 1 || dT123 == 3) "dT123 must be 1 or 3"
                T12min = ifelse(dT123==3,1,0)
                Nalpha = 0
                for J12 = 0:J12max
                    for S12 = 0:1
                        for L12 = abs(J12-S12):J12+S12 
                            for T12 = T12min:1
                                if (-1)^(L12+S12+T12) != -1; continue;end                    
                                for dj3 = max(1,abs(dJ123-2*J12)):2:dJ123+2*J12
                                    for l3 = div(max(0,2*dj3-1),2):div(2*dj3+1,2)                                    
                                        if (-1)^(L12+l3) != P123; continue; end
                                        Nalpha += 1
                                        ich += 1
                                        push!(chs,channel3b_Jj(dJ123,dT123,P123,L12,S12,J12,T12,dj3,l3))
                                    end
                                end
                            end
                        end
                    end
                end
                dim = (Nalpha * Np * Nq)^2
                if verbose
                    println("J12max $(@sprintf("%2i",J12max)): J=$(@sprintf("%2i",dJ123))/2  ",
                        "T=$(@sprintf("%2i",dT123))/2  P=$(@sprintf("%2i",P123))",
                        "  Nα $(@sprintf("%4i",Nalpha))  Dim.=> ",@sprintf("%6.e",dim))
                end
            end
        end
    end
    return chs
end

function make_pq(ps,wps,qs,wqs)
    dim = length(ps) * length(qs)
    pqs = zeros(Float64,dim,2)
    wpqs = zeros(Float64,dim)
    idx = 0
    for (ip,p) in enumerate(ps)
        for (iq,q) in enumerate(qs)
            idx += 1
            pqs[idx,1] = p
            pqs[idx,2] = q
            wpqs[idx] = wps[ip] * wqs[iq]
        end
    end
    return pqs,wpqs,dim
end

function Jacobi_3NF_Contact(LECs,chan3b,ch_idxs,params,Mats_pq,vec_freg_nl)
    cE = LECs.dLECs["cE"] 
    fac_cont = hc * (hc/Fpi)^4 * (hc/params.LambdaChi) / (12.0 *sqrt(3.0) *pi^4)
    E = cE * fac_cont
    pqs = params.meshpoints.pqs
    pqdim = size(pqs)[1]
    alphadim = length(ch_idxs)
    println("dim. of |pq;alpha>: $(pqdim*alphadim)")

    @threads for idx in eachindex(ch_idxs)
        ch_bra,ch_ket = ch_idxs[idx]
        tid = threadid()
        mat = Mats_pq[tid]
        bra = chan3b[ch_bra]
        ket = chan3b[ch_ket]

        dT123 = bra.dT123
        L12 = bra.L12
        S12 = bra.S12
        J12 = bra.J12
        T12 = bra.T12
        dj3 = bra.dj3
        l3 = bra.l3
        if L12 != 0; continue; end
        if l3 != 0; continue;end

        dT123_k = ket.dT123
        L12_k = ket.L12
        S12_k = ket.S12
        J12_k = ket.J12 
        T12_k = ket.T12 
        dj3_k = ket.dj3
        l3_k = ket.l3
 
        if dT123_k != dT123; continue;end
        if L12_k != 0; continue; end
        if l3_k != 0; continue;end
        if dj3 != dj3_k; continue;end
        if S12 != S12_k; continue; end
        if J12 != J12_k; continue; end
    
        tdot,rnav = IsospinDot(T12,T12_k,dT123)
        me = tdot * E
        if tdot == 0.0;continue;end

        #println("idx $idx ch_bra $ch_bra ch_ket $ch_ket T12 $T12 T12_k $T12_k tdot $tdot rnav $rnav E $E")

        for pqidx_bra = 1:pqdim
            p_bra,q_bra = pqs[pqidx_bra,:]
            freg_bra = vec_freg_nl[pqidx_bra]
            for pqidx_ket = 1:pqdim
                p_ket,q_ket = pqs[pqidx_ket,:]
                freg_ket = vec_freg_nl[pqidx_ket]
                mat[pqidx_bra,pqidx_ket] = freg_bra * me * freg_ket
            end
        end
    end
    return nothing
end
