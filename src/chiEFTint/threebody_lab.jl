struct sps_nljtz
    n::Int64
    l::Int64
    j::Int64
    e::Int64
    tz::Int64
    isocc::Bool
end

struct sps_nlj
    n::Int64
    l::Int64
    j::Int64
    e::Int64
end

struct ket_labHO
    p::sps_nljtz
    q::sps_nljtz
    r::sps_nljtz
    Tz::Int64
    J::Int64
    T::Int64
end

function prep_lab_space(params)
    e3max = params.e3max
    N3max = params.N3max
    J3max = params.j3max
    J12max = params.J12max
    emax = params.e2max
    e2max = emax * 2
    println("emax $emax e2max $e2max")
    J12max = J3max #ad hoc

    nsps = div(emax*(emax+1),2)
    spsISO = sps_nlj[ ]
    for e = 0:emax
        for l = 0:e
            if (e-l)%2 != 0; continue;end
            n = div(e-l,2)
            for j = abs(2*l-1):2:2*l+1
                push!(spsISO,sps_nlj(n,l,j,e))
            end
        end
    end
    dict_lab_JPT = Dict{Int64,Vector{Int64}}()
    for Jtot = 1:2:max(J3max,2*e3max+3)
        for Ptot = 1:-2:-1
            for Ttot = 1:2:3
                Ncnt = Dict{Int64,Int64}(); for N = 0:N3max; Ncnt[N] = 0; end
                for i1 = 1:nsps
                    o1 = spsISO[i1]; l1 = o1.l; j1 = o1.j
                    for i2 = 1:i1
                        o2 = spsISO[i2]; l2 = o2.l; j2 = o2.j
                        for i3 = 1:i2
                            o3 = spsISO[i3]; l3 = o3.l; j3 = o3.j
                            if o1.e+o2.e > e2max || o1.e+o3.e > e2max || o2.e+o3.e > e2max;continue;end
                            if o1.e+o2.e+o3.e > e3max; continue;end                            
                            if (-1)^(l1+l2+l3) != Ptot; continue; end
                            N = o1.e+o2.e+o3.e 
                            for j12 = div(abs(j1-j2),2):div(j1+j2,2)
                                if !tri_check(j12,j3/2,Jtot/2);continue;end
                                for t12 = 0:1
                                    if !tri_check(t12,1/2,Ttot/2);continue;end
                                    if i1 == i2 && (j12 +t12)%2==0;continue;end
                                    Ncnt[N] += 1
                                end
                            end
                        end
                    end       
                end
                for N = 0:N3max
                    tkey = get_nkey4(Jtot,Ptot,Ttot,N)
                    dict_lab_JPT[tkey] = [ Ncnt[N],0,0 ]
                end
            end
        end
    end
    return dict_lab_JPT
end

function const_labHO(params,JPTs,dict_lab_JPT,dict_JPTN_dim,idx_dict_JPT)
    #j12      L12      lam     lcm
    #     Jpq      lr       l3      Jrel
    # J        jr       sr      j3 
    for (idx_JPT,JPT) in enumerate(JPTs)
        J,P,T = JPT
        ch_kets_labHO(params,J,P,T,dict_lab_JPT,dict_JPTN_dim,idx_dict_JPT)
    end
end

function ch_kets_labHO(params,J,P,T,dict_lab_JPT,dict_JPTN_dim,idx_dict_JPT)
    e3max = params.e3max
    J3max = params.j3max
    labHO_mat = zeros(Int64,1,1)
    labHO_mat123 = zeros(Int64,1,1)
    for loop = 1:2
        ndim = ndim123 = 0
        for n123 = 0:e3max
            for E = 0:n123
                ecm = n123 - E           
                for ncm = 0:div(ecm,2)
                    lcm = ecm - 2*ncm                 
                    for jrel = abs(J-2*lcm):2:min(J3max,2*e3max+3,J+2*lcm)
                        for prel = -1:2:1
                            if P != prel * (-1)^(lcm); continue; end
                            idx_rel = idx_dict_JPT[[jrel,prel,T]]
                            no = dict_JPTN_dim["orth"][idx_rel][E][1]                                  
                            if no == 0; continue; end
                            if loop == 2
                                for j = 1:no
                                    ndim += 1
                                    tv = @view labHO_mat[:,ndim]
                                    tv[1] = E; tv[2] = prel; tv[3] = jrel
                                    tv[4] = j; tv[5] = ncm;  tv[6] = lcm
                                end
                            else
                                ndim += no
                            end
                        end
                    end
                end
            end                
            for e12 = 0:n123
                for e3 = 0:n123-e12
                    if e12 + e3 > e3max; continue;end
                    ecm = n123-e12-e3
                    for n12 = 0:div(e12,2)
                        l12 = e12 - 2*n12
                        for n3 = 0:div(e3,2)
                            l3 = e3 - 2*n3
                            for ncm = 0:div(ecm,2)
                                lcm = ecm - 2*ncm
                                for s12 = 0:1
                                    for j12 = abs(l12-s12):l12+s12
                                        for j3 = abs(2*l3-1):2:2*l3+1
                                            for t12 = 0:1
                                                if !tri_check(t12,1/2,T/2); continue;end
                                                if (-1)^(l12+s12+t12) != -1;continue;end
                                                if (-1)^(l12+l3+lcm) != P; continue;end
                                                for jrel = abs(2*j12-j3):2:min(2*j12+j3,J3max)
                                                    if !tri_check(lcm,jrel/2,J/2);continue;end
                                                    ndim123 += 1
                                                    if loop == 2
                                                        tv = @view labHO_mat123[:,ndim123]
                                                        tv .= [n12,l12,s12,j12,t12,n3,l3,j3,ncm,lcm,jrel,e12+e3]
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
            end
        end
        if loop == 1
            labHO_mat = zeros(Int64,6,ndim)
            labHO_mat123 = zeros(Int64,12,ndim123)
            nlab = 0 
            for N = 0:params.N3max
                tkey = get_nkey4(J,P,T,N)
                nlab += dict_lab_JPT[tkey][1]
            end
            tkey = get_nkey4(J,P,T,params.N3max)                
            dict_lab_JPT[tkey][2] = ndim
            dict_lab_JPT[tkey][3] = ndim123
            println("J",@sprintf("%3i",J)," P",@sprintf("%3i",P)," T ",@sprintf("%3i",T)," nlab $nlab ndim $ndim ndim123 $ndim123")
        end
    end
    return nothing
end
