struct single_sps
    e::Int
    n::Int
    l::Int
    j2::Int
    tz::Int
end

struct sps_3Blab
    e1max::Int
    e1max_file::Int
    e2max_file::Int
    e3max::Int
    e3max_file::Int
    norbits::Int
    norbits_file::Int
    sps::Dict{Int64,single_sps}
    sps_file::Dict{Int64,single_sps}
end

struct Obj_3BME
    use3BME::Bool
    sps_3b::sps_3Blab
    dict_3b_idx::Dict{UInt64,Int64}
    dict_idx_to_snt::Dict{Int,Int}
    dict_idx_to_me3j::Dict{Int,Int}
    v3bme::Vector{Float64}
    v3monopole::Dict{UInt64,Float64}
end

function main_read_me3j(fn_3nf, e1max, e1max_file, e2max_file, e3max, e3max_file, sps_snt, dWS, to) 
    sps_3b = get_modelspace(e1max, e1max_file, e2max_file, e3max, e3max_file)
    if fn_3nf == "" 
        return Obj_3BME(false, sps_3b, Dict{UInt64,Int64}(),Dict{Int,Int}(), Dict{Int,Int}(), [0.0], Dict{UInt64,Float64}())
    end      
    println("ModelSpace: {e1max $e1max e3max $e3max}  File: {e1max $e1max_file e2max $e2max_file e3max $e3max_file}")

    @timeit to "Read" begin
        println("Reading 3nf... \n$fn_3nf")
        @timeit to "count_nreads_File" dict_idxThBME = count_nreads(sps_3b,"File",to)
        @timeit to "count_me3jgz" count_ME_file = count_me3jgz(sps_3b)
        println("count_ME (File) $count_ME_file")
        @timeit to "read_me3jgz" ThBME = read_me3jgz(fn_3nf, count_ME_file, to)
    end

    @timeit to "alloc/store" begin
        println("Allocating v3bme...")
        @timeit to "count_nreads_MS" nreads_v3bme = count_nreads(sps_3b,"ModelSpace",to)
        v3bme,dict_3b_idx= allocate_3bme(sps_3b)
        println("Storing 3BME to v3bme...")
        if e1max_file == e1max && e2max_file == 2*e1max && e3max_file == e3max
            v3bme .= ThBME
        else
            @timeit to "store" store_me3jgz!(sps_3b, ThBME,v3bme, nreads_v3bme, dWS, dict_idxThBME)
        end
    end

    println("Building V3monopole...")
    @timeit to "Build V3monopole" V3mono = monopole_V3(e3max, sps_3b,dict_3b_idx,v3bme,dWS)
    dict_idx_to_snt, dict_idx_to_me3j = get_dict_idx_to_snt(sps_snt, sps_3b)
    return Obj_3BME(true, sps_3b, dict_3b_idx, dict_idx_to_snt, dict_idx_to_me3j, v3bme, V3mono)
end

"""
$(SIGNATURES)

Function to read me3j.gz using GZip. The values are stored in a vector ThBME.
The ordering of `ThBME` is not considered here.
"""
function read_me3jgz(fn,count_ME_file, to; verbose=false)
    isfile(fn) || error("File not found: $fn")
    size_ME = count_ME_file*8/1024^3
    @assert size_ME < 0.9 * (Sys.total_memory() / 2^20 /1024) "# of ThBME=$(size_ME) is beyond available memory"
    ThBME = zeros(Float64,count_ME_file)    
    stream = GzipDecompressorStream(open(fn))
    @inbounds for (idx,line) in enumerate(eachline(stream))
        idx_i = 1 + (idx-2)*10
        idx_f = idx_i + 9
        subsize = 10
        if idx == 1;continue; end
        if idx_i > count_ME_file
            break
        end
        if idx_f > count_ME_file
            subsize = count_ME_file - idx_i + 1
        end
        for i = 1:subsize
            tl = @view line[16*(i-1)+1:16*i]
            ThBME[idx_i+i-1] = Parsers.parse(Float64, tl)
        end
    end
    close(stream)
    return ThBME
end

function get_modelspace(e1max, e1max_file, e2max_file, e3max, e3max_file;verbose=false, l_descending_order=false)
    sps = Dict{Int64,single_sps}()
    sps_file = Dict{Int64,single_sps}()
    norbits_file = norbits_ms = 0
    for mode in ["File", "ModelSpace"]
        norbits = 0
        target = ifelse(mode=="File", sps_file, sps)
        for temax = 0:ifelse(mode=="File", e1max_file, e1max)
            lmin, lstep, lmax = temax%2, 2, temax
            if l_descending_order
                lmin, lstep, lmax = temax, -2, temax%2
            end
            for l = lmin:lstep:lmax
                n = div(temax-l,2)
                for j2 = abs(2*l-1):2:2*l+1
                    for tz = -1:2:1
                        norbits += 1
                        target[norbits] = single_sps(temax, n, l, j2, tz)
                        if verbose; println("norbits $norbits e = $temax  n $n l $l j2 $j2 tz $tz");end
                    end
                end
            end
        end
        if mode == "File"
            norbits_file = norbits
        else
            norbits_ms = norbits
        end
    end
    return sps_3Blab(e1max, e1max_file, e2max_file, e3max, e3max_file, norbits_ms, norbits_file, sps, sps_file)
end

"""
Allocating 3BME vector and dictionary for the indices.
Note that the dimension of `v3bme` is the one that is used for the 3BME (for modelspace) instead of the number of elements in the input ThBME file.
"""
function allocate_3bme(sps_3b, ME_is_double=true)
    norbits = sps_3b.norbits
    sps = sps_3b.sps
    e1max = sps_3b.e1max
    e3max = sps_3b.e3max
    dict_3b_idx = Dict{UInt64,Int64}()
    total_dim = 0
    for a = 1:2:norbits
        ea = sps[a].e
        la = sps[a].l
        if ea > e1max; continue; end
        if ea > e3max; continue; end
        for b = 1:2:a
            eb = sps[b].e
            lb = sps[b].l
            if ea+eb > e3max; continue; end
            Jab_min = div(abs(sps[a].j2 - sps[b].j2),2)
            Jab_max = div(abs(sps[a].j2 + sps[b].j2),2)
            for c = 1:2:b
                ec = sps[c].e
                lc = sps[c].l
                if ea+eb+ec > e3max; continue; end
                for d = 1:2:a
                    ed = sps[d].e
                    ld = sps[d].l
                    for e = 1:2:ifelse(d==a,b,d)
                        ee = sps[e].e
                        le = sps[e].l
                        for f = 1:2:ifelse((a==d && b==e),c,e)
                            ef = sps[f].e
                            lf = sps[f].l
                            if ed+ee+ef > e3max; continue; end
                            if (la + lb + lc + ld + le + lf) % 2 != 0; continue; end
                            orbit_hash = get_nkey6(a,b,c,d,e,f)
                            dict_3b_idx[orbit_hash] = total_dim
                            Jde_min = div(abs(sps[d].j2 - sps[e].j2),2)
                            Jde_max = div(abs(sps[d].j2 + sps[e].j2),2)
                            for Jab = Jab_min:Jab_max
                                for Jde = Jde_min:Jde_max
                                    J2_min = max( abs(2*Jab-sps[c].j2), abs(2*Jde-sps[f].j2) )
                                    J2_max = min( 2*Jab+sps[c].j2, 2*Jde+sps[f].j2 )
                                    for J2 = J2_min:2:J2_max
                                        total_dim += 5
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    size_3bme = total_dim * ifelse(ME_is_double, 8, 4) / 1024.0^3
    println("# of 3BME: ", @sprintf("%12i", total_dim), " Mem. ", @sprintf("%12.5e", size_3bme), " GB")
    @assert size_3bme < 0.9 * (Sys.total_memory() / 2^20 /1024) "size(3BME) $(size_3bme) is beyond your environment"
    v3bme = zeros(Float64,total_dim)    
    return v3bme, dict_3b_idx
end

"""

While the snt format is given ascending order for ``\\ell`` and ``j``
the 3BME is stored in descending order for ``\\ell`` and ``j``.
Hence, we need to reorder the indices for the 3BME.
Note that this is ad hoc and not guaranteed to work for all cases.
"""
function get_dict_idx_to_snt(sps_snt, sps_3b)
    dict_idx_me3j_to_snt = Dict{Int,Int}()
    dict_idx_snt_to_me3j = Dict{Int,Int}()
    count_d = 0
    emax = sps_3b.e1max
    for temax = 0:emax
        lmin, lstep, lmax = temax, -2, temax%2
        for l = lmin:lstep:lmax
            n = div(temax-l,2)
            for j2 = 2*l+1:-2:abs(2*l-1)
                for tz = -1:2:1
                    count_d += 1
                    for idx_s = 1:length(sps_snt)
                        nn, ll, jj, tzz = sps_snt[idx_s].n, sps_snt[idx_s].l, sps_snt[idx_s].j, sps_snt[idx_s].tz
                        if nn == n && ll == l && jj == j2 && tzz == tz                            
                            dict_idx_me3j_to_snt[count_d] = idx_s
                            dict_idx_snt_to_me3j[idx_s] = count_d
                            break
                        end
                    end
                end
            end
        end
    end
    return dict_idx_me3j_to_snt, dict_idx_snt_to_me3j
end

"""
Since we gonna store in isospin formalism, we use only proton part.
"""
function count_nreads(sps_3b, mode, to;verbose=false, save_verbose=false)
    io = stdout 
    norbits = ifelse(mode=="File", sps_3b.norbits_file, sps_3b.norbits)
    e1max = ifelse(mode=="File", sps_3b.e1max_file, sps_3b.e1max)
    e2max = ifelse(mode=="File", sps_3b.e2max_file, sps_3b.e1max*2)
    e3max = ifelse(mode=="File", sps_3b.e3max_file, sps_3b.e3max)
    sps = ifelse(mode=="File", sps_3b.sps_file, sps_3b.sps)
    dict_idx_ThBME = Dict{UInt64,Int64}()
    nread = 0    
    nreads = zeros(Int64,div(norbits,2))
    @timeit to "loop" for idx_a = 1:2:norbits
        nreads[div(idx_a,2)+1] = nread
        oa = sps[idx_a]
        ea = oa.e
        if ea > e1max; continue; end
        for idx_b = 1:2:idx_a 
            ob = sps[idx_b]
            eb = ob.e
            if ea+eb > e2max; continue; end
            for idx_c = 1:2:idx_b
                oc = sps[idx_c]
                ec = oc.e
                if ea+eb+ec > e3max; continue; end

                # J_bra < J_ket
                JabMax = div(oa.j2 + ob.j2,2)
                JabMin = div(abs(oa.j2 - ob.j2),2)
                twoJCMindownbra = 0
                if abs(oa.j2 -ob.j2) > oc.j2
                    twoJCMindownbra = abs(oa.j2 -ob.j2) - oc.j2
                elseif oc.j2 < (oa.j2+ob.j2)
                    twoJCMindownbra = 1
                else
                    twoJCMindownbra = oc.j2 - oa.j2 - ob.j2
                end
                twoJCMaxupbra = oa.j2 + ob.j2 + oc.j2 

                # loop for ket 
                for idx_d = 1:2:idx_a
                    od = sps[idx_d]
                    ed = od.e
                    for idx_e =1:2:ifelse(idx_a==idx_d,idx_b,idx_d)
                        oe = sps[idx_e]
                        ee = oe.e
                        idx_f_max = ifelse( (idx_a==idx_d && idx_b==idx_e), idx_c, idx_e)
                        for idx_f = 1:2:idx_f_max
                            of = sps[idx_f]
                            ef = of.e
                            if ed+ee+ef > e3max; continue; end
                            if (oa.l + ob.l + oc.l + od.l + oe.l + of.l) % 2 != 0; continue; end

                            JdeMax = div(od.j2+oe.j2,2)
                            JdeMin = div(abs(od.j2-oe.j2),2)
                            twoJCMindownket = 0
                            if abs(od.j2 -oe.j2) > of.j2
                                twoJCMindownket = abs(od.j2 -oe.j2) - of.j2
                            elseif of.j2 < (od.j2+oe.j2)
                                twoJCMindownket = 1
                            else
                                twoJCMindownket = of.j2 - od.j2 - oe.j2
                            end
                            twoJCMaxupket = od.j2 + oe.j2 + of.j2

                            twoJCMindown = max(twoJCMindownbra, twoJCMindownket)
                            twoJCMaxup = min(twoJCMaxupbra, twoJCMaxupket)
                            if twoJCMindown > twoJCMaxup; continue;end
                            if mode == "File"
                                nkey = get_nkey6(idx_a,idx_b,idx_c,idx_d,idx_e,idx_f) 
                                dict_idx_ThBME[ nkey ] = nread
                            end
                            for Jab = JabMin:JabMax
                                for Jde = JdeMin:JdeMax
                                    twoJCMin = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2) )
                                    twoJCMax = min( 2*Jab+oc.j2, 2*Jde+of.j2 )
                                    if twoJCMin > twoJCMax; continue; end                                   
                                    blocksize = (div(twoJCMax - twoJCMin,2) + 1)*5
                                    nread += blocksize
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if mode == "File"
        nkeys = length(keys(dict_idx_ThBME))
        println("size of dict_idx_ThBME: $nkeys ",nkeys*2*8/1024^3," GB ")
        return dict_idx_ThBME
    else
        return nreads
    end
end

"""
Counting the offset for the 3BME indices.
"""
function count_me3jgz(sps_3b::sps_3Blab;mode="File")    
    e1max = sps_3b.e1max
    e1max_file = sps_3b.e1max_file 
    e2max_file = sps_3b.e2max_file 
    e3max_file = sps_3b.e3max_file 
    e3max = sps_3b.e3max
    e1max_check = ifelse(mode=="File", e1max_file, e1max)
    e3max_check = ifelse(mode=="File", e3max_file, e3max)
    sps = ifelse(mode=="File", sps_3b.sps_file, sps_3b.sps)
    norbits = ifelse(mode=="File", sps_3b.norbits_file, sps_3b.norbits)
    count_ME_file = 0
    for idx_a = 1:2:norbits
        oa = sps[idx_a]
        ea = oa.e
        if ea > e1max_check; continue; end
        for idx_b = 1:2:idx_a 
            ob = sps[idx_b]
            eb = ob.e
            if ea+eb > e2max_file; continue; end
            for idx_c = 1:2:idx_b
                oc = sps[idx_c]
                ec = oc.e
                if ea+eb+ec > e3max_check; continue; end
                JabMax = div(oa.j2 + ob.j2,2)
                JabMin = div(abs(oa.j2 - ob.j2),2)
                twoJCMindownbra = 0
                if abs(oa.j2 -ob.j2) > oc.j2
                    twoJCMindownbra = abs(oa.j2 -ob.j2) - oc.j2
                elseif oc.j2 < (oa.j2+ob.j2)
                    twoJCMindownbra = 1
                else
                    twoJCMindownbra = oc.j2 - oa.j2 - ob.j2
                end
                twoJCMaxupbra = oa.j2 + ob.j2 + oc.j2 
               
                for idx_d = 1:2:idx_a
                    od = sps[idx_d]
                    ed = od.e
                    if ed > e1max_check; continue;end
                    for idx_e =1:2:ifelse(idx_a==idx_d,idx_b,idx_d)
                        oe = sps[idx_e]
                        ee = oe.e
                        if ee > e1max_check; continue;end
                        idx_f_max = ifelse( (idx_a==idx_d && idx_b==idx_e), idx_c, idx_e)
                        for idx_f = 1:2:idx_f_max
                            of = sps[idx_f]
                            ef = of.e
                            if ef > e1max_check; continue;end
                            if ed+ee+ef > e3max_check; continue; end                            
                            if (oa.l + ob.l + oc.l + od.l + oe.l + of.l) % 2 != 0; continue; end
                    
                            JdeMax = div(od.j2 + oe.j2,2)
                            JdeMin = div(abs(od.j2 - oe.j2),2)
                            twoJCMindownket = 0
                            if abs(od.j2 -oe.j2) > of.j2
                                twoJCMindownket = abs(od.j2 -oe.j2) - of.j2
                            elseif of.j2 < (od.j2+oe.j2)
                                twoJCMindownket = 1
                            else
                                twoJCMindownket = of.j2 - od.j2 - oe.j2
                            end
                            twoJCMaxupket = od.j2 + oe.j2 + of.j2

                            twoJCMindown = max(twoJCMindownbra, twoJCMindownket)
                            twoJCMaxup = min(twoJCMaxupbra, twoJCMaxupket)
                            if twoJCMindown > twoJCMaxup; continue;end

                            for Jab = JabMin:JabMax
                                for Jde = JdeMin:JdeMax
                                    twoJCMin = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2) )
                                    twoJCMax = min( 2*Jab+oc.j2, 2*Jde+of.j2 )
                                    if twoJCMin > twoJCMax; continue; end
                                   
                                    blocksize = (div(twoJCMax - twoJCMin,2) + 1)*5
                                    count_ME_file += blocksize
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return count_ME_file
end

function valid_check(ea,eb,ec,ed,ee,ef,e1max,e2max,e3max)
    if ea > e1max; return false; end
    if eb > e1max; return false; end
    if ec > e1max; return false; end
    if ed > e1max; return false; end
    if ee > e1max; return false; end
    if ef > e1max; return false; end
    if ea+eb > e2max; return false; end
    if ea+ec > e2max; return false; end
    if ed+ee > e2max; return false; end
    if ed+ef > e2max; return false; end
    if ee+ef > e2max; return false; end
    if ea+eb+ec > e3max; return false; end
    if ed+ee+ef > e3max; return false; end
    return true
end

function store_me3jgz!(sps_3b::sps_3Blab, ThBME,  v3bme, nreads_v3bme, dWS, dict_idxThBME)
    e1max = sps_3b.e1max
    e2max = e1max*2
    e1max_file = sps_3b.e1max_file 
    e2max_file = sps_3b.e2max_file 
    e3max_file = sps_3b.e3max_file 
    e3max = sps_3b.e3max
    l3max = e1max
    sps = sps_3b.sps
    norbits = sps_3b.norbits 
    count_ME_file = 0
    for idx_a = 1:2:norbits
        oa = sps[idx_a]
        ea = oa.e
        nread_v3bme = nreads_v3bme[div(idx_a,2)+1]
        #@assert nread_v3bme == nread_ThBME "nread_v3bme $(nread_v3bme) nread_ThBME $(nread_ThBME)"
        if ea > e1max; continue; end
        for idx_b = 1:2:idx_a 
            ob = sps[idx_b]
            eb = ob.e
            if ea+eb > e2max; continue; end
            for idx_c = 1:2:idx_b
                oc = sps[idx_c]
                ec = oc.e
                if ea+eb+ec > e3max; continue; end
                JabMax = div(oa.j2 + ob.j2,2)
                JabMin = div(abs(oa.j2 - ob.j2),2)
                twoJCMindownbra = 0
                if abs(oa.j2 -ob.j2) > oc.j2
                    twoJCMindownbra = abs(oa.j2 -ob.j2) - oc.j2
                elseif oc.j2 < (oa.j2+ob.j2)
                    twoJCMindownbra = 1
                else
                    twoJCMindownbra = oc.j2 - oa.j2 - ob.j2
                end
                twoJCMaxupbra = oa.j2 + ob.j2 + oc.j2 
               
                # loop for ket 
                for idx_d = 1:2:idx_a
                    od = sps[idx_d]
                    ed = od.e
                    if ed > e1max; continue;end
                    for idx_e =1:2:ifelse(idx_a==idx_d,idx_b,idx_d)
                        oe = sps[idx_e]
                        ee = oe.e
                        if ee > e1max; continue;end
                        idx_f_max = ifelse( (idx_a==idx_d && idx_b==idx_e), idx_c, idx_e)
                        for idx_f = 1:2:idx_f_max
                            of = sps[idx_f]
                            ef = of.e
                            if ef > e1max; continue;end
                            if ed+ee+ef > e3max; continue; end                            
                            if (oa.l + ob.l + oc.l + od.l + oe.l + of.l) % 2 != 0; continue; end                  
                           
                            valid_for_v3bme = valid_check(ea,eb,ec,ed,ee,ef,e1max,e2max,e3max)
                            if !valid_for_v3bme; continue; end
                            valid_for_ThBME = valid_check(ea,eb,ec,ed,ee,ef,e1max_file,e2max_file,e3max_file)
                            if !valid_for_ThBME; continue; end

                            JdeMax = div(od.j2 + oe.j2,2)
                            JdeMin = div(abs(od.j2 - oe.j2),2)
                            twoJCMindownket = 0
                            if abs(od.j2 -oe.j2) > of.j2
                                twoJCMindownket = abs(od.j2 -oe.j2) - of.j2
                            elseif of.j2 < (od.j2+oe.j2)
                                twoJCMindownket = 1
                            else
                                twoJCMindownket = of.j2 - od.j2 - oe.j2
                            end
                            twoJCMaxupket = od.j2 + oe.j2 + of.j2

                            twoJCMindown = max(twoJCMindownbra, twoJCMindownket)
                            twoJCMaxup = min(twoJCMaxupbra, twoJCMaxupket)
                            if twoJCMindown > twoJCMaxup; continue;end

                            offset_ThBME = dict_idxThBME[ get_nkey6(idx_a,idx_b,idx_c,idx_d,idx_e,idx_f)]
                            idx_ThBME = offset_ThBME 
                            
                            for Jab = JabMin:JabMax
                                for Jde = JdeMin:JdeMax
                                    twoJCMin = max( abs(2*Jab-oc.j2), abs(2*Jde-of.j2) )
                                    twoJCMax = min( 2*Jab+oc.j2, 2*Jde+of.j2 )
                                    if twoJCMin > twoJCMax; continue; end
                                    blocksize = (div(twoJCMax - twoJCMin,2) + 1)*5
                                    for JTind = 0:twoJCMax-twoJCMin+1
                                        twoJC = twoJCMin + div(JTind,2)*2
                                        twoT = 1 + (JTind%2)*2
                                        for tab = 0:1
                                            for tde = 0:1
                                                if twoT > min(2*tab+1,2*tde+1); continue;end
                                                index_ab = div(5*(twoJC-twoJCMin),2)+2*tab+tde+div(twoT-1,2)
                                                v3idx = nread_v3bme + index_ab + 1
                                                idx_ThBME += 1
                                                ThBME_idx = idx_ThBME 
                                                V = 0.0
                                                autozero = false
                                                if (oa.l > l3max || ob.l > l3max || oc.l > l3max || od.l > l3max || oe.l > l3max || of.l > l3max )
                                                    V = 0.0
                                                end

                                                v3bme[v3idx] = ThBME[ThBME_idx] 
                                                if (idx_a==idx_b && (tab+Jab)%2==0) || (idx_d==idx_e && (tde+Jde)%2==0); autozero=true;; end
                                                if (idx_a==idx_b && idx_a==idx_c && twoT==3 && oa.j2<3); autozero = true; end
                                                if (idx_d==idx_e && idx_d==idx_f && twoT==3 && od.j2<3); autozero = true; end
                                                # if (autozero && V > 1.e-8) 
                                                #     @error "This should not happen V $V autozero $autozero"
                                                # end
                                            end
                                        end
                                    end                                   
                                    if valid_for_v3bme
                                        nread_v3bme += blocksize
                                    end
                                    count_ME_file += blocksize
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return nothing
end

function monopole_V3(E3max, sps_3b,dict_3b_idx,v3bme,dWS)
    n_orbits = sps_3b.norbits    
    sps = sps_3b.sps
    keys = UInt64[ ]
    for i = 1:n_orbits
        oi = sps[i]
        ei = oi.e 
        for j = i:n_orbits
            ej = sps[j].e
            oj = sps[j]
            if oi.l != oj.l || oi.j2 != oj.j2 || oi.tz != oj.tz; continue; end
            for a = 1:n_orbits
                oa = sps[a]
                ea = sps[a].e
                for b = 1:n_orbits
                    ob = sps[b]
                    eb = sps[b].e
                    if oa.l != ob.l || oa.j2 != ob.j2 || oa.tz != ob.tz; continue; end
                    for c = 1:n_orbits 
                        ec = sps[c].e
                        oc = sps[c]
                        if ea + ec + ei > E3max; continue; end
                        for d = 1:n_orbits
                            ed = sps[d].e
                            od = sps[d]
                            if oc.l != od.l || oc.j2 != od.j2 || oc.tz != od.tz; continue; end
                            if eb + ed + ej > E3max; continue; end
                            if (oi.l + oa.l + ob.l + oc.l + od.l + oj.l) % 2 != 0; continue; end
                            key = get_nkey6(a,c,i,b,d,j)
                            push!(keys ,key)
                        end
                    end
                end
            end
        end
    end
    nkeys = length(keys)
    println("Number of keys for V3mono: ", nkeys)
    Vmon3 = Dict{UInt64,Float64}()
    for idx = 1:nkeys
        key = keys[idx]
        Vmon3[key] = 0.0
    end
    @threads for idx = 1:nkeys
        key = keys[idx]
        Vmon3[key] = 0.0
        a,c,i,b,d,j = unhash_key6j(key)
        @assert key == get_nkey6(a,c,i,b,d,j) "key mismatch"
        ja = sps[a].j2
        jc = sps[c].j2
        ji = sps[i].j2
        jb = sps[b].j2
        jd = sps[d].j2
        jj = sps[j].j2
        j2min = div(max( abs(ja-jc), abs(jb-jd) ),2)
        j2max = div(min( ja+jc, jb+jd ),2)
        v = 0.0
        for j2 = j2min:j2max
            Jmin = max( abs(2*j2-ji), abs(2*j2-jj) )
            Jmax = 2 *j2 + min(ji,jj)
            for J2 = Jmin:2:Jmax
                vtmp = get_V3_pn(idx, E3max, v3bme,j2,j2,J2,a,c,i,b,d,j,sps_3b,dict_3b_idx,dWS) *(J2+1) 
                v += vtmp                 
            end
        end
        Vmon3[key] += v /(ji+1)
        #vallsum += Vmon3[key]
        #if 4 <= idx <= 5 && Vmon3[idx] != 0.0
        # if idx <= 10 && Vmon3[idx] != 0.0
        #     println("idx $idx  vmon3 ", Vmon3[idx])
        # end        
    end
    println("norm(Vmon3) = ", norm(values(Vmon3)),)
    return Vmon3
end

function get_V3_pn(indx, E3max, v3bme,Jab,Jde,J2,a,b,c,d,e,f,sps_3b,dict_3b_idx,dWS;V_in=0.0)
    tza = sps_3b.sps[a].tz
    tzb = sps_3b.sps[b].tz
    tzc = sps_3b.sps[c].tz
    tzd = sps_3b.sps[d].tz
    tze = sps_3b.sps[e].tz
    tzf = sps_3b.sps[f].tz
    dcg_spin = dWS.dcg_spin
    Vpn = 0.0
    Tmin = max( abs(tza+tzb+tzc), abs(tzd+tze+tzf) )
    for tab = div(abs(tza+tzb),2):1
        CG1 = dcg_spin[ get_nkey6_shift(1,tza,1,tzb,tab*2, tza+tzb) ]
        for tde = div(abs(tzd+tze),2):1
            CG2 = dcg_spin[ get_nkey6_shift(1,tzd,1,tze,tde*2, tzd+tze) ]
            if CG1*CG2 == 0; continue; end
            Tmax = min(1+2*tab, 1+2*tde)
            for T2 = Tmin:2:Tmax
                CG3 = dcg_spin[ get_nkey6_shift(tab*2, (tza+tzb), 1, tzc, T2, (tza+tzb+tzc)) ]
                CG4 = dcg_spin[ get_nkey6_shift(tde*2, (tzd+tze), 1, tzf, T2, (tzd+tze+tzf)) ]
                if CG3*CG4 == 0; continue; end
                tbme = Get3BME_ISO(indx, E3max, v3bme,dict_3b_idx,sps_3b,Jab,Jde,J2,tab,tde,T2,a,b,c,d,e,f,dWS,V_in)
                Vpn += (CG1*CG2*CG3*CG4) * tbme
            end
        end
    end 
    return Vpn
end

"""
The (re)order is done with only odd(tz=-1) indices.
abc=0, bca=1, cab=2, acb=3, bac=4, cba=5

"""
function sort_3_orbits(a_in,b_in,c_in)
    a_in = ifelse(a_in%2==0,a_in-1,a_in)
    b_in = ifelse(b_in%2==0,b_in-1,b_in)
    c_in = ifelse(c_in%2==0,c_in-1,c_in)
    a = a_in; b = b_in; c = c_in
    if a < b; a,b = b,a; end
    if b < c; b,c = c,b; end
    if a < b; a,b = b,a; end
    idx = -1    
    if a_in == a 
        idx = ifelse(b_in == b,0,3)
    elseif a_in == b
        idx = ifelse(b_in == a,4,1)
    else
        idx = ifelse(b_in == a,2,5)
    end
    return a, b, c, idx
end 

function RecouplingCG(idx_abc,ja2,jb2,jc2,Jab_in,Jab,J2, dWS)::Float64
    if div(abs(ja2-jb2),2) > Jab || div(ja2+jb2,2) < Jab; return 0.0; end
    if div(abs(jc2-J2),2) > Jab || div(jc2+J2,2) < Jab; return 0.0; end
    if idx_abc == 0
        return ifelse(Jab==Jab_in,1.0,0.0)
    elseif idx_abc == 1 # bca        
        phase = (-1)^( div(jb2+jc2,2)+Jab_in + 1 )
        t6j = dWS.d6j_lj[ get_key6j_sym(ja2,jb2,Jab*2,jc2,J2,Jab_in*2)]
        return phase * hat(Jab_in) * hat(Jab) * t6j
    elseif idx_abc == 2 # cab
        phase = (-1)^( div(ja2+jb2,2)-Jab+1)
        t6j = dWS.d6j_lj[ get_key6j_sym(jb2,ja2,Jab*2,jc2,J2,Jab_in*2)]
        return phase * hat(Jab_in) * hat(Jab) * t6j
    elseif idx_abc == 3 # acb
        phase = (-1)^( div(jb2+jc2,2)+Jab_in-Jab)
        t6j = dWS.d6j_lj[ get_key6j_sym(jb2,ja2,Jab*2,jc2,J2,Jab_in*2)]
        ret = phase * hat(Jab_in) * hat(Jab) * t6j
        return ret 
    elseif idx_abc == 4 # bac
        if Jab == Jab_in
            phase = (-1)^(div(ja2+jb2,2)-Jab)
            return 1.0 * phase
        else
            return 0.0
        end
    elseif idx_abc == 5 # cba
        t6j = dWS.d6j_lj[ get_key6j_sym(ja2,jb2,Jab*2,jc2,J2,Jab_in*2)]
        return - hat(Jab_in) * hat(Jab) * t6j
    else
        @assert false "This should not happen"
    end
end

""" 

Function to get the 3BME from the vector v3bme.
To this end, we need to calculate the index of the 3BME in the vector v3bme from given `Jab,Jde,J2,tab,tde,T,a,b,c,d,e,f`.

Since the dict for 3BME, the order if set so that `a>=b>=c, d>=e>=f`,
one needs to be careful to get the correct index for the 3BME.
"""
function Get3BME_ISO(ind, E3max, v3bme,dict_3b_idx,sps_3b,
                     Jab_in,Jde_in,J2,tab_in,tde_in,T2,
                     a_in,b_in,c_in,d_in,e_in,f_in,
                     dWS,V_in)
    sps = sps_3b.sps
    v = 0.0
    a,b,c,idx_abc = sort_3_orbits(a_in,b_in,c_in)
    d,e,f,idx_def = sort_3_orbits(d_in,e_in,f_in)

    if d > a  || (d==a && e>b) || (d==a && e==b && f>c)
        a,b,c,d,e,f = d,e,f,a,b,c
        Jab_in, Jde_in = Jde_in, Jab_in
        tab_in, tde_in = tde_in, tab_in
        idx_abc,idx_def = idx_def,idx_abc
    end

    tkey = get_nkey6(a,b,c,d,e,f)
    idx_3borbit = get(dict_3b_idx,tkey,-1) + 1
    idx_3borbit == 0 && return 0.0
   
    oa = sps[a]; ob = sps[b]; oc = sps[c]
    od = sps[d]; oe = sps[e]; of = sps[f]
    if oa.e + ob.e + oc.e >  E3max; return 0.0; end
    if od.e + oe.e + of.e >  E3max; return 0.0; end
  
    ja2 = oa.j2; jb2 = ob.j2; jc2 = oc.j2
    jd2 = od.j2; je2 = oe.j2; jf2 = of.j2
    Jab_min = div(abs(ja2-jb2),2); Jab_max = div(abs(ja2+jb2),2)
    Jde_min = div(abs(jd2-je2),2); Jde_max = div(abs(jd2+je2),2)

    tab_min = ifelse(T2==3,1,0); tab_max = 1
    tde_min = ifelse(T2==3,1,0); tde_max = 1

    J_index = count_inner = 0
    for Jab=Jab_min:Jab_max
        Cj_abc = RecouplingCG(idx_abc,ja2,jb2,jc2,Jab_in,Jab,J2,dWS)
        if div(idx_abc,3) != div(idx_def,3) # for odd permutation
            Cj_abc *= -1.0
        end
        for Jde=Jde_min:Jde_max
            Cj_def = RecouplingCG(idx_def,jd2,je2,jf2,Jde_in,Jde,J2,dWS)
            J2_min = max( abs(2*Jab-jc2), abs(2*Jde-jf2) )
            J2_max = min( 2*Jab+jc2, 2*Jde+jf2 )
            if J2_min > J2_max; continue; end

            J_index += div(J2-J2_min,2)*5 
            count_inner += 1

            if J2 >= J2_min && J2 <= J2_max && abs(Cj_abc*Cj_def) > 1.e-10
                for tab = tab_min:tab_max
                    Ct_abc = RecouplingCG(idx_abc,1,1,1,tab_in,tab,T2,dWS)        
                    for tde = tde_min:tde_max
                        Ct_def = RecouplingCG(idx_def,1,1,1,tde_in,tde,T2,dWS)
                        if abs(Ct_abc*Ct_def) < 1.e-10; continue; end
                        Tindex = 2*tab + tde + div(T2-1,2)
                        idx = idx_3borbit + J_index + Tindex                        
                        v += Cj_abc * Cj_def * Ct_abc * Ct_def * v3bme[idx]
                    end
                end
            end
            J_index += div(J2_max-J2+2,2)*5
        end
    end
    return v
end

