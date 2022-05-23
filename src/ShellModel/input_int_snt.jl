struct abcdJ
    a::Int64
    b::Int64
    c::Int64
    d::Int64
    tJ::Int64
end
function make_int(p_sps,n_sps)
    count = 0
    label_T1 = abcdJ[ ]
    label_T0 = abcdJ[ ]
    tV = [ ]
    for (i,tmp) in enumerate(p_sps)
        n1 = tmp[1]; l1=tmp[2]; j1=tmp[3]
        for (j,tmp2) in enumerate(p_sps)            
            if j<i;continue;end
            n2 = tmp2[1]; l2=tmp2[2]; j2=tmp2[3]
            for tJ in div(abs(j1-j2),2):1:div(j1+j2,2)
                if i==j && tJ%2==1;continue;end
                push!(tV,[i,j,tJ])
            end
        end
    end
    for i = 1:length(tV)
        for j = i:length(tV)
            if tV[i][3] != tV[j][3];continue;end
            a,b = tV[i][1:2]
            c,d = tV[j][1:2]
            tJ = tV[j][3]
            push!(label_T1,abcdJ(a,b,c,d,tJ))
            count += 1
        end
    end
    count = 0
    tV = [ ]
    ofst = length(p_sps)
    for (i,tmp) in enumerate(p_sps)
        n1 = tmp[1]; l1=tmp[2]; j1=tmp[3]
        i = i 
        for (j,tmp2) in enumerate(n_sps)            
            j = j 
            if j < i;continue;end
            n2 = tmp2[1]; l2=tmp2[2]; j2=tmp2[3]
            for tJ in div(abs(j1-j2),2):1:div(j1+j2,2)
                if i==j && tJ%2==0;continue;end
                push!(tV,[i,j,tJ])
            end
        end
    end
    for i = 1:length(tV)
        for j = i:length(tV)
            if tV[i][3] != tV[j][3];continue;end
            a,b = tV[i][1:2]
            c,d = tV[j][1:2]
            tJ = tV[j][3]
            count += 1
            push!(label_T0,abcdJ(a,b,c,d,tJ))
        end
    end
    return label_T1,label_T0
end

function int2snt(p_sps,n_sps,label_T1,label_T0,olabels)
    l_TBME = length(olabels)    
    ofst = length(p_sps)
    s_from_i = [ Int64[-1,-1] for i = 1:l_TBME]
    facs = [ Float64[0.0,0.0] for i = 1:l_TBME]    
    for (j,olabel) in enumerate(olabels)
        a,b,c,d,totJ = @views olabel[1:5]
        fac0 = 0.0; fac1 = 0.0;  rank = "--"

        if (a <= ofst && b <= ofst && c <= ofst && d <= ofst) || (
            a > ofst && b > ofst && c > ofst && d > ofst)
            fac1 = 1.0
            rank = "pp/nn"
            for (i,label) in enumerate(label_T1) 
                ta = label.a; tb = label.b; tc = label.c; td = label.d; tJ = label.tJ
                if tJ != totJ; continue;end
                if (a == ta && b == tb && c == tc && d == td ) || (
                    a == ta +ofst && b == tb+ofst && c == tc+ofst && d == td+ofst) 
                    s_from_i[j][1] = i
                    facs[j][1] = fac1
                    break
                end
            end
        else
            rank = "pn"
            for (i,label) in enumerate(label_T1) 
                ta = label.a; tb = label.b; tc = label.c; td = label.d; tJ = label.tJ
                if totJ != tJ; continue;end
                coef = 0.5 * sqrt( (1+deltaf(ta,tb) * (-1)^tJ) * (1+deltaf(tc,td) * (-1)^tJ))
                if (a == ta && b -ofst == tb && c == tc && d-ofst == td) || (
                    a == tc && b -ofst == td && c == ta && d-ofst == tb) #abcd or cdab
                    fac1 = coef
                elseif (a == tb && b -ofst == ta && c == tc && d-ofst == td) || (
                    a == tc && b -ofst == td && c == tb && d-ofst == ta) #bacd or cdba
                    ja = p_sps[tb][3]; jb = p_sps[ta][3]
                    fac1 = (-1)^(div(ja+jb,2)+tJ+1) .* coef
                elseif (a == ta && b -ofst == tb && c == td && d-ofst == tc) || (
                    a == td && b -ofst == tc && c == ta && d-ofst == tb) #abdc or dcab
                    jc = p_sps[td][3]; jd = p_sps[tc][3]
                    fac1 = (-1)^(div(jc+jd,2)+tJ+1) .* coef
                elseif (a == tb && b -ofst == ta && c == td && d-ofst == tc) || (
                    a == td && b -ofst == tc && c == tb && d-ofst == ta) #badc or dcba
                    ja = p_sps[tb][3]; jb = p_sps[ta][3]
                    jc = p_sps[td][3]; jd = p_sps[tc][3]
                    fac1 = (-1)^(div(ja+jb+jc+jd,2)) .* coef
                else
                    continue
                end
                s_from_i[j][1]=i
                facs[j][1] = fac1
                break
            end
            for (i,label) in enumerate(label_T0) 
                ta = label.a; tb = label.b; tc = label.c; td = label.d; tJ = label.tJ
                if totJ != tJ; continue;end
                coef = 0.5 * sqrt( (1-deltaf(ta,tb) * (-1)^tJ) * (1-deltaf(tc,td) * (-1)^tJ))
                if (a == ta && b-ofst == tb && c == tc && d-ofst == td) || (
                    a == tc && b-ofst == td && c == ta && d-ofst == tb)  #abcd or cdab
                    fac0 = coef
                elseif (a == tb && b-ofst == ta && c == tc && d-ofst == td) || (
                    a == tc && b-ofst == td && c == tb && d-ofst == ta) #bacd or cdba
                    ja = p_sps[tb][3]; jb = p_sps[ta][3]
                    fac0 = (-1)^(div(ja+jb,2)+tJ) * coef
                elseif (a == ta && b-ofst == tb && c == td && d-ofst == tc) || (
                    a == td && b-ofst == tc && c == ta && d-ofst == tb) #abdc or dcab
                    jc = p_sps[td][3]; jd = p_sps[tc][3]
                    fac0 = (-1)^(div(jc+jd,2)+tJ) * coef
                elseif (a == tb && b-ofst == ta && c == td && d-ofst == tc) || (
                    a == td && b-ofst == tc && c == tb && d-ofst == ta) #badc or dcba
                    ja = p_sps[tb][3]; jb = p_sps[ta][3]
                    jc = p_sps[td][3]; jd = p_sps[tc][3]
                    fac0 = (-1)^(div(ja+jb+jc+jd,2)) * coef
                else                    
                    continue
                end
                s_from_i[j][2]=i
                facs[j][2] = fac0
                break
            end
        end
    end
    return s_from_i,facs
end
function readVint(intf,label_T1,label_T0,TF=true) # ad hoc func. to read from **.int
    f = open(intf,"r");tlines = readlines(f);close(f)
    lines = rm_comment(tlines)
    SPEs = Float64[ ]
    ofst = 0
    V0s = zeros(Float64,length(label_T0))
    V1s = zeros(Float64,length(label_T1))

    if intf == ""
        SPEs = zeros(Float64,length(p_sps))
        for i = 1:length(label_T0);V0s[i] = randn(); end
        for i = 1:length(label_T1);V1s[i] = randn(); end
        return SPEs,V1s,V0s
    end
    for line in lines
        if length(split(line)) < 7
            push!(SPEs,parse(Float64,split(line)[2]))
            continue
        end       
        if TF;ofst = -length(SPEs);end
        tmp = split(line)
        a = parse(Int64,tmp[1])+ofst
        b = parse(Int64,tmp[2])+ofst
        c = parse(Int64,tmp[3])+ofst
        d = parse(Int64,tmp[4])+ofst
        J = parse(Int64,tmp[5])
        T = parse(Int64,tmp[6])
        V = parse(Float64,tmp[7])
        if T==0
            for (idx,ref) in enumerate(label_T0)
                ta = ref.a; tb = ref.b; tc = ref.c; td = ref.d; tJ = ref.tJ
                if a == ta && b==tb && c==tc && d==td && tJ== J
                    V0s[idx] = V                    
                    break
                end
            end
        else
            for (idx,ref) in enumerate(label_T1)
                ta = ref.a; tb = ref.b; tc = ref.c; td = ref.d; tJ = ref.tJ
                if a == ta && b==tb && c==tc && d==td && tJ== J
                    V1s[idx] = V                    
                    break
                end
            end
        end
    end
    return SPEs,V1s,V0s
end
function evalVsnt!(Vsnt,SPEs,V1s,V0s,idx_s_from_i,facs)
    ## adhoc: the same model space for proton/neutron 
    #Vsnt = zeros(Float64,length(SPEs)*2 + length(idx_s_from_i))
    ofst = length(SPEs)
    for i = 1:length(SPEs)
        Vsnt[i] = SPEs[i]
        Vsnt[i+ofst] = SPEs[i]
    end
    ofst = 2*length(SPEs)
    fac0 = 0.0; fac1 = 0.0
    for (n,idxs) in enumerate(idx_s_from_i)
        idx_T1,idx_T0 = idxs
        Vsnt[ofst+n] = 0.0
        if idx_T1 != -1
            fac1 = facs[n][1]  
            Vsnt[ofst+n] += V1s[idx_T1] .* fac1
        end
        if idx_T0 != -1 
            fac0 = facs[n][2]            
            Vsnt[ofst+n] += V0s[idx_T0] .* fac0
        end
    end
    return nothing
end

function evalVint!(Vint,SPEs,V1s,V0s)
    ln = length(SPEs)
    @inbounds for i = 1:ln
        Vint[i] = SPEs[i]
    end
    @inbounds for i = 1:length(V1s)
        Vint[ln+i] = V1s[i]
    end
    ln += length(V1s)
    @inbounds for i = 1:length(V0s)
        Vint[ln+i] = V0s[i]
    end
    return nothing
end

function proposal_Vint!(nSPEs,nV1s,nV0s,SPEs,V1s,V0s,varM)
    @inbounds for i = 1:length(SPEs)
        nSPEs[i] = SPEs[i] + varM * randn()
    end
    @inbounds for i = 1:length(V1s)
        nV1s[i] = V1s[i] + varM * randn()
    end
    @inbounds for i = 1:length(V0s)
        nV0s[i] = V0s[i] + varM * randn()
    end
    return nothing
end

# function test()
#     sntf = "../snts/usdb.snt"
#     Anum = 28
#     lp,ln,cp,cn,massop,Aref,pow,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsnt(sntf,Anum)
    
#     label_T1,label_T0 = make_int(p_sps,n_sps)
#     idx_s_from_i,facs = int2snt(p_sps,n_sps,label_T1,label_T0,olabels)

#     intf = "../snts/random_input/ints/sdshl/tmp_1.int"
#     SPEs,V1s,V0s = readVint(intf,label_T1,label_T0)
#     Vsnt = evalVsnt(SPEs,V1s,V0s,idx_s_from_i,facs)
# end
#test()
