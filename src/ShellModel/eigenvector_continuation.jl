function J2_to_J(J2)
    J=""
    if Int(J2) % 2 == 0
        J=string(Int(div(J2,2)))
    else
        J=string(Int(J2))*"/2"
    end
    return J
end

function myCholesky!(tmpA,ln::Int64,cLL::LowerTriangular{Float64,Array{Float64,2}})
    l11 = sqrt(tmpA[1,1]) 
    cLL[1,1] = l11
    cLL[2,1] = tmpA[2,1]/l11; cLL[2,2] = sqrt( tmpA[2,2]-cLL[2,1]^2)
    for i=3:ln
        for j=1:i-1
            cLL[i,j] = tmpA[i,j]
            for k = 1:j-1
                cLL[i,j] += - cLL[i,k]*cLL[j,k]                
            end
            cLL[i,j] = cLL[i,j] / cLL[j,j]            
        end
        cLL[i,i] = tmpA[i,i]
        for j=1:i-1
            cLL[i,i] += -cLL[i,j]^2
        end
        cLL[i,i] = sqrt(cLL[i,i])             
    end
    nothing
end

function read_wfinfos(inpf)
    f = open(inpf,"r");tlines = readlines(f);close(f)
    lines = rm_comment(tlines)
    ret = []
    for line in lines
        tl = rm_nan(split(line," "))
        push!(ret,[tl[1],parse(Int64,tl[3]),parse(Int64,tl[5])])
    end
    return ret
end

function read_tdmat(inpf)
    io=open(inpf,"r")
    Dim=read(io,Int)
    ln = read(io,Int)
    lTD = read(io,Int)
    TDs = [ Float64[] for i=1:ln]
    infos = [ Int64[] for i=1:ln]
    Nmat = zeros(Float64,Dim,Dim)
    for nth = 1:ln
        i = read(io,Int); j = read(io,Int); Nij = read(io,Int)
        infos[nth] = [i,j,Nij]
        tN = read(io,Float64)
        Nmat[i,j] = tN
        Nmat[j,i] = tN
        TDs[nth] = [ read(io,Float64) for k=1:lTD]
    end
    close(io)
    return Dim,ln,lTD,Nmat,infos,TDs
end

function read_wf_from_info(wfinfos,mdim,vs,wpath)
    ## e.g., ./wavsamples_sigma3/Mg24_sample_1.wav nth= 1 TDmatidx  1
    for tmp in wfinfos
        inpf,nth,TDidx = tmp
        if wpath != "./" && wpath != ""
            tn = split(inpf,"/")
            ln =length(tn)
            inpf = wpath*tn[ln-1]*"/"*tn[ln]
        end
        io = open(inpf,"r")
        num_ev = read(io,Int32)
        mtot = read(io,Int32)
        vals = [read(io,Float64) for i = 1:num_ev]
        ojs = [read(io,Int32) for i=1:num_ev]    
        Css = [ [read(io,Float64) for i=1:mdim] for j=1:num_ev]
        #println("TDidx $TDidx nth $nth")
        vs[TDidx] .= Css[nth]
        close(io)
    end
    return nothing
end

function readRvecs!(mdim,inpf,tJ,vecs,wfinfos;num=100)
    io = open(inpf,"r")
    num_ev = read(io,Int32)
    mtot = read(io,Int32)
    vals = [read(io,Float64) for i = 1:num_ev]
    ojs = [read(io,Int32) for i=1:num_ev]
    Css = [ [read(io,Float64) for i=1:mdim] for j=1:num_ev]
    n = minimum([num,num_ev])
    for i=1:n
        j = Int64(ojs[i])
        if j == tJ
            push!(vecs,Css[i])
            push!(wfinfos,[inpf,i,length(vecs)])
        end
    end
    close(io)
    return nothing
end

function prepEC(Hs,target_nuc,num_ev,num_ECsample,tJ,mode;
                num_history=3,lm=100,ls=15,tol=1.e-9,q=1,
                path_to_samplewav="",calc_moment=true,
                save_wav = false,tdmatdir="./tdmat/",
                gfactors = [1.0,0.0,5.586,-3.826],
                effcharge=[1.5,0.5])
    to = TimerOutput()
    sntf = Hs[1]    
    Anum = parse(Int64, match(reg,target_nuc).match)
    lp,ln,cp,cn,massop,Aref,pow,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsnt(sntf,Anum)
    hw, bpar = init_ho_by_mass(Anum,1) # mass formula 
    if 16<= Anum <= 40
        hw, bpar = init_ho_by_mass(Anum,2) # 2: mass formula for sd-shell
    end
    Mtot = 0;if Anum % 2 != 0; Mtot = 1;end
    if typeof(tJ) != Int64;println("targetJ must be integer!!");exit();end    
    Mtot = tJ
    eval_jj = 0.5*tJ*(tJ/2+1)
    if Anum % 2 != Mtot % 2; println("invalid target J");exit();end
    
    target_el = replace.(target_nuc, string(Anum)=>"")
    Z,N,vp,vn = getZNA(target_el,Anum,cp,cn)
    mstates_p,mstates_n,mz_p,mz_n = def_mstates(p_sps,n_sps)
    pbits,nbits,jocc_p,jocc_n,Mps,Mns,tdims = occ(p_sps,mstates_p,mz_p,vp,
                                                  n_sps,mstates_n,mz_n,vn,Mtot)
    oPP,oNN,oPNu,oPNd = prep_J(tdims,p_sps,n_sps,mstates_p,mstates_n,
                               pbits,nbits)
    mdim = tdims[end]
    lblock=length(pbits)
    mdim_print(target_nuc,Z,N,cp,cn,vp,vn,mdim,tJ)
    if mode == "sample"
        for (iter,sntf) in enumerate(Hs)
            @timeit to "make sample" begin
                println("sntf: $sntf")
                lp,ln,cp,cn,massop,Aref,pow,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsnt(sntf,Anum)
                bV1,V1 = HbitT1(p_sps,n_sps,mstates_p,mstates_n,labels,TBMEs)
                bVpn,Vpn,delMs = Hbitpn(p_sps,n_sps,mstates_p,mstates_n,
                                        labels[3],TBMEs[3])
                ppinfo = prep_pp(mstates_p,pbits,bV1[1],V1[1])
                nninfo = prep_nn(mstates_n,nbits,bV1[2],V1[2])
                bV1 = nothing
                l_pbit = length(mstates_p);l_nbit = length(mstates_n)
                bis,bfs,p_NiNfs,n_NiNfs,num_task = prep_pn(lblock,tdims,
                                                           l_pbit,l_nbit,
                                                           pbits,nbits,
                                                           Mps,delMs,bVpn,Vpn)
                bVpn=nothing
                block_tasks = make_distribute(num_task)

                Js = [ 0.5*Mtot*(0.5*Mtot+1) for i = 1:num_ev]
                Jtasks = zeros(Int64,lblock)
                for i = 1:lblock
                    Jtasks[i] = length(pbits[i])*length(nbits[i])
                end
                Jidxs = make_distribute_J(Jtasks)
            
                en =[ [1.e+4 for j=1:num_ev] for i = 1:num_history]
                #Rvecs = [ zeros(Float64,mdim) for i=1:num_ev]
                Tmat = zeros(Float64,lm,lm)
                itmin = 1; elit=1
                vks = [ zeros(Float64,mdim) for i=1:lm]
                uks = [ zeros(Float64,mdim) for i=1:ls]                
                initialize_wf(vks[1],"rand",tJ,mdim)
                @timeit to "Lanczos" elit = TRL(vks,uks,Tmat,itmin,
                                                pbits,nbits,jocc_p,jocc_n,SPEs,
                                                ppinfo,nninfo,bis,bfs,block_tasks,
                                                p_NiNfs,n_NiNfs,Mps,delMs,Vpn,
                                                eval_jj,oPP,oNN,oPNu,oPNd,Jidxs,
                                                tdims,num_ev,num_history,lm,ls,en,
                                                tol,to,true)
                Rvecs = @views uks[1:num_ev]               
                vals,vecs = eigen(@views Tmat[1:elit*q,1:elit*q])
                @inbounds for (nth,Rvec) in enumerate(Rvecs)
                    Rvec .= 0.0
                    @inbounds for k=1:length(vals)
                        Rvec .+= vecs[k,nth] .* vks[k]
                    end
                    Rvec .*= 1.0/sqrt(dot(Rvec,Rvec))    
                end
                Jv = zeros(Float64,mdim)
                for (nth,Rv) in enumerate(Rvecs)
                    Jv .= 0.0
                    operate_J!(Rv,Jv,pbits,nbits,tdims,
                               Jidxs,oPP,oNN,oPNu,oPNd)
                    Js[nth] += dot(Rv,Jv)
                end
                totJs = J_from_JJ1.(Js)
                println("J $totJs")#;println("Js $Js")
                print("En. ");map(x -> @printf("%24.15f ",x), en[1]);print("\n")
                csnt = split(split(sntf,"/")[end],".")[1]
                oupf= path_to_samplewav*target_nuc*"_"*csnt*".wav"
                if tJ != -1;oupf=path_to_samplewav*target_nuc*"_"*csnt*"_j"*string(tJ)*".wav";end
                if save_wav
                    writeRitzvecs(mdim,Mtot,en[1],totJs,Rvecs,oupf)
                end
                if calc_moment                    
                    @timeit to "moment" tx_mom = eval_moment(Mtot,Rvecs,totJs,p_sps,n_sps,
                                                          mstates_p,mstates_n,tdims,
                                                          jocc_p,jocc_n,pbits,nbits,bpar,
                                                          gfactors,effcharge)
                    if tx_mom != "";print(tx_mom);end
                end
                println("")
            end
        end
    elseif mode =="TD"
        try
            mkdir(tdmatdir)
        catch
            nothing
        end
        #@timeit to "read sample wavs." begin
        svecs = [ zeros(Float64,1) ];deleteat!(svecs,1)
        wfinfos = [ ]       
        for fidx = 0:num_ECsample-1
            inpf = path_to_samplewav*target_nuc*"_tmp_$fidx"*"_j"*string(tJ)*".wav"
            readRvecs!(mdim,inpf,tJ,svecs,wfinfos;num=num_ev) 
        end
        
        io=open(tdmatdir*"wfinfos_"*target_nuc*"_J"*string(tJ)*".dat","w")
        for tmp in wfinfos
            println(io,split(tmp[1],"/")[end], " nth= ",tmp[2], " TDmatidx  ",tmp[3])
        end
        close(io)
        @timeit to "misc" begin
            Dim = length(svecs)
            Hmat = zeros(Float64,Dim,Dim)               
            Nmat = zeros(Float64,Dim,Dim)
            println("Dim $Dim")
        end
        @timeit to "Nmat calc." begin
            @threads for j = 1:Dim
                v = svecs[j]
                for i = 1:j
                    tN = dot(svecs[i],v)
                    Nmat[i,j] = tN; Nmat[j,i]=tN[1]
                end
            end
        end
        @timeit to "OBTD" begin
            numv = length(svecs)
            idxs = Tridx[] 
            for i = 1:length(svecs)
                for j=i:length(svecs)
                    push!(idxs, Tridx(i,j,idx_from_ij(i,j,numv)))
                end
            end                    
            OBTDs = [ [ zeros(Float64,length(p_sps)),
                        zeros(Float64,length(n_sps))
                        ] for Nij = 1:length(idxs)  ]
            calcOBTD(OBTDs,idxs,p_sps,n_sps,
                     mstates_p,mstates_n,tdims,
                     jocc_p,jocc_n,pbits,nbits,svecs)
        end
        
        @timeit to "prepTBTD" begin
            opTBTD1,pjumps,njumps,bif_idxs = prepTBTD(tJ,idxs,p_sps,n_sps,
                                                      mstates_p,mstates_n,
                                                      pbits,nbits,
                                                      labels,TBMEs,svecs,
                                                      tdims,Mps)
            TBTDs= [zeros(Float64,length(oTBMEs)) for Nij = 1:length(idxs)]            
        end                
        @timeit to "calcTBTD" begin
            calcTBTD(TBTDs,opTBTD1,pjumps,njumps,
                     pbits,nbits,tdims,svecs,
                     idxs,bif_idxs,
                     olabels,oTBMEs,labels,to)
        end
        @timeit to "write TDs" begin
            oupf=tdmatdir*target_nuc*"_J"*string(tJ)*".dat"
            lTD = length(p_sps)+length(n_sps)+length(oTBMEs)
            io=open(oupf,"w")
            write(io,Dim)
            write(io,length(idxs))
            write(io,lTD)
            for i = 1:length(svecs)
                for j=i:length(svecs)
                    Nij = idx_from_ij(i,j,length(svecs))
                    write(io,i);write(io,j)
                    write(io,Nij)
                    write(io,Nmat[i,j])              
                    TDs = Float64[]
                    for nth=1:length(OBTDs[Nij][1])
                        push!(TDs,OBTDs[Nij][1][nth])
                    end
                    for nth=1:length(OBTDs[Nij][2])
                        push!(TDs,OBTDs[Nij][2][nth])
                    end
                    for nth=1:length(TBTDs[Nij])
                        push!(TDs,TBTDs[Nij][nth]) 
                    end
                    for TD in TDs
                        write(io,TD)
                    end
                end
            end                    
            close(io)
        end
    end
    show(to, allocations = true, compact = true)
    println("\n\n")
    return nothing
end

"""
    solveEC(Hs,target_nuc,tJNs;
            write_appwav=false,verbose=false,calc_moment=true,wpath="./",is_show=false,
            gfactors = [1.0,0.0,5.586,-3.826],effcharge=[1.5,0.5],exact_logf="")

To solve EC (generalized eigenvalue problem) to approximate the eigenpairs for a given interaction.
```math
H \\vec{v} = \\lambda N \\vec{v}
```
Transition densities and overlap matrix for H and N are read from "tdmat/" directory (To be changed to more flexible)

# Arguments:
- `Hs`:array of paths to interaction files (.snt fmt)
- `target_nuc`: target nucleus
- `tJNs`:array of target total J (doubled) and number of eigenstates to be evaluated
    e.g., [ [0,5], [2,5] ], when you want five lowest J=0 and J=1 states.

# Optional arguments:
- `write_appwav=false`:write out the approximate wavefunction
- `verbose=false`:to print (stdout) approx. energies for each interaction
- `calc_moment=true`: to calculate mu&Q moments
- `wpath="./"`: path to sample eigenvectors to construct approx. w.f.
- `is_show=false`: to show TimerOutput
- `gfactors=[1.0,0.0,5.586,-3.826]`: g-factors to evaluate moments
- `effcharge=[1.5,0.5]`:effective charges to evaluate moments

# Optional arguments for author's own purpose
- `exact_logf=""`:path to logfile for E(EC) vs E(Exact) plot

!!! note
    All the effective interactions must be in "the same order" and must be consistent with interaction file from which the transition density matrices were made.

"""
function solveEC(Hs,target_nuc,tJNs;
                 write_appwav=false,verbose=false,
                 calc_moment=true,wpath="./",tdmatpath="./",
                 is_show=false,
                 gfactors = [1.0,0.0,5.586,-3.826],
                 effcharge=[1.5,0.5],
                 exact_logf="")

    if length(tJNs)==0;println("specify targetJ and # of states");exit();end
    to = TimerOutput()
    sntf = Hs[1]    
    Anum = parse(Int64, match(reg,target_nuc).match)
    lp,ln,cp,cn,massop,Aref,pow,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsnt(sntf,Anum)
    hw, bpar = init_ho_by_mass(Anum,1) # mass formula 
    if 16<= Anum <= 40
        hw, bpar = init_ho_by_mass(Anum,2) # 2: mass formula for sd-shell
    end
    target_el = replace.(target_nuc, string(Anum)=>"")
    Z,N,vp,vn = getZNA(target_el,Anum,cp,cn)  
    mstates_p,mstates_n,mz_p,mz_n = def_mstates(p_sps,n_sps)
    #println("nuc: $target_el")
    exlines = ""
    try 
        f = open(exact_logf,"r");exlines = readlines(f);close(f)
    catch
        nothing
        #println("logfile of exact calc. is missing ")
    end

    Dims = Int64[]
    sumV=[]
    wfinfos=[]
    lTD = 0    
    for (jidx,tmp) in enumerate(tJNs)
        tJ,num_ev_target = tmp
        Mtot = tJ
        pbits,nbits,jocc_p,jocc_n,Mps,Mns,tdims = occ(p_sps,mstates_p,mz_p,vp,
                                                      n_sps,mstates_n,mz_n,vn,Mtot)
        mdim = tdims[end]        
        @timeit to "read TDmat" begin
            inpf = tdmatpath*target_nuc*"_J"*string(tJ)*".dat"
            Dim,ln,lTD,Nmat,infos,TDs = read_tdmat(inpf)
            push!(Dims,Dim)
            if Dim+1 <= num_ev_target;
                println("Error!: sample number must be larger than $num_ev_target")
                exit()
            end                
        end
        #println("for J=$tJ/2 states Dim.=$Dim")        
        @timeit to "Cholesky" begin
            tMat = zeros(Float64,Dim,Dim)
            tildH = zeros(Float64,Dim,Dim)
            L = LowerTriangular(zeros(Float64,Dim,Dim))
            try
                myCholesky!(Nmat,Dim,L)
            catch
                u = 1.1 * 1.e-16
                shift = tr(Nmat) * u * ( (Dim+2)/(1.0-((Dim+1)*(Dim+2))*u))
                #println("#epsilon prescription is used!! shift: $shift" )
                Nmat = Nmat + shift .* Matrix(1.0I, Dim,Dim)
                myCholesky!(Nmat,Dim,L)
            end
            Linv = inv(L)
        end
        Hmat = zeros(Float64,Dim,Dim)

        for sntidx = 1:length(Hs)
            sntf = Hs[sntidx]
            tMat .= 0.0
            @timeit to "Hmat" begin
                lp,ln,cp,cn,massop,Aref,pow,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsnt(sntf,Anum)
                MEs = [SPEs[1],SPEs[2],oTBMEs]
                MEs = [(MEs...)...]
                for (idx,TD) in enumerate(TDs)
                    i,j,Nij = infos[idx]
                    tmp = dot(TD,MEs)
                    Hmat[i,j] = tmp; Hmat[j,i] = tmp
                end
                @timeit to "Gen. Eigen" begin
                    mul!(tMat,Linv,Hmat)
                    mul!(tildH,Linv',tMat)                                    
                    vals,vecs = real.(Arpack.eigs(tildH,nev=num_ev_target,which=:SR,
                                                  tol=1.e-8, maxiter=300))
                    if verbose
                        print_vec("$target_nuc 2J=$tJ  En(EC)",vals)
                    end
                end                
                push!(sumV,[sntf,tJ,vals])
            end
        end
            
        ## write out approx. wav.
        if write_appwav
            wfinfos = read_wfinfos("./tdmat/wfinfos_"*target_nuc*"_J"*string(tJ)*".dat")
            @timeit to "read sample wavs." begin
                vs = [zeros(Float64,mdim) for i=1:Dim]
                read_wf_from_info(wfinfos,mdim,vs,wpath)
            end
            @timeit to "make .appwav" begin                        
                Rvecs = [ zeros(Float64,mdim) for i=1:length(vals)]
                for nth = 1:length(vals)
                    for k = 1:length(vs)
                        @. Rvecs[nth] += vecs[k,nth] * vs[k]
                    end
                    Rvecs[nth] ./= sqrt(sum(Rvecs[nth].^2))
                end
            end            
            @timeit to "write .appwav" begin
                csnt = split(split(sntf, "/")[end],".")[1]
                wavname = "./appwavs/"*target_nuc*"_"*csnt*"_j"*string(tJ)*".appwav"
                Jlist = [tJ for i=1:length(vals)]
                writeRitzvecs(mdim,Mtot,vals,Jlist,Rvecs,wavname)
            end
        end                
    end
    if length(Hs) > 1 && exact_logf != "" 
        @timeit to "scatter plot" plot_EC_scatter(target_nuc,Hs,sumV,tJNs,Dims,exlines)
    end
    if is_show
        show(to, allocations = true, compact = false)
    end
    #println("")
    return nothing
end


function solveEC_UQ(Hs,target_nuc,tJNs,Erefs,errors;
                    itnum_MCMC = 5000,burnin=500,var_M=2.e-2,
                    calc_moment=true,wpath="./",
                    is_show=false,
                    gfactors = [1.0,0.0,5.586,-3.826],
                    effcharge=[1.5,0.5],
                    exact_logf="",
                    intf="",
                    num_replica=10)
    UQ=true; write_appwav=false
    if length(tJNs)==0;println("specify targetJ and # of states");exit();end
    to = TimerOutput()
    
    sntf = Hs[1]    
    Anum = parse(Int64, match(reg,target_nuc).match)
    lp,ln,cp,cn,massop,Aref,pow,p_sps,n_sps,SPEs,olabels,oTBMEs,labels,TBMEs = readsnt(sntf,Anum)
    hw, bpar = init_ho_by_mass(Anum,1) # mass formula 
    if 16<= Anum <= 40
        hw, bpar = init_ho_by_mass(Anum,2) # 2: mass formula for sd-shell
    end
    target_el = replace.(target_nuc, string(Anum)=>"")
    Z,N,vp,vn = getZNA(target_el,Anum,cp,cn)  
    mstates_p,mstates_n,mz_p,mz_n = def_mstates(p_sps,n_sps)
    println("nuc: $target_nuc")

    Dims = Int64[]
    lTD = 0
    lnJ = length(tJNs)
    HNmats = [ [ zeros(Float64,2,2) ] for i =1:lnJ]
    for i=1:lnJ;deleteat!(HNmats[i],1);end
    AllTDs = [ [ zeros(Float64,2) ] ]; deleteat!(AllTDs,1)
    Allinfos = [ [ [0,0] ] ]; deleteat!(Allinfos,1)

    for (jidx,tmp) in enumerate(tJNs)
        tJ,num_ev_target = tmp
        Mtot = tJ
        @timeit to "read TDmat" begin
            inpf = "./tdmat/"*target_nuc*"_J"*string(tJ)*".dat"
            Dim,ln,lTD,Nmat,infos,TDs = read_tdmat(inpf)
            push!(Dims,Dim)
            if Dim+1 <= num_ev_target;
                println("Error!: sample number must be larger than $num_ev_target")
                exit()
            end                
        end
        @timeit to "Cholesky" begin
            tMat = zeros(Float64,Dim,Dim)
            tildH = zeros(Float64,Dim,Dim)
            L = LowerTriangular(zeros(Float64,Dim,Dim))
            try
                myCholesky!(Nmat,Dim,L)
            catch
                u = 1.1 * 1.e-16
                shift = tr(Nmat) * u * ( (Dim+2)/(1.0-((Dim+1)*(Dim+2))*u))
                println("#epsilon prescription is used!! shift: $shift" )
                Nmat = Nmat + shift .* Matrix(1.0I, Dim,Dim)
                myCholesky!(Nmat,Dim,L)
            end
            Linv = inv(L)
        end
        Hmat = zeros(Float64,Dim,Dim)
        push!(AllTDs,copy(TDs))
        push!(Allinfos,copy(infos))
        push!(HNmats[jidx],Hmat)
        push!(HNmats[jidx],Linv)
    end
    
    ## MCMC
    in_snt =false
    Theta = zeros(Float64,length(SPEs[1])+length(SPEs[2])+length(oTBMEs))
    lME = length(Theta)    
    tDim = Dims[1]
    tMat = zeros(Float64,tDim,tDim)
    tildH = zeros(Float64,tDim,tDim)
    label_T1,label_T0 = make_int(p_sps,n_sps)
    idx_s_from_i,facs = int2snt(p_sps,n_sps,label_T1,label_T0,olabels)
    
    ### random walk in .snt space (not used now)    
    #MEs = [SPEs[1],SPEs[2],oTBMEs]
    #Theta .= [(MEs...)...]
    #Thetas = [ zeros(Float64,lME) ]; deleteat!(Thetas,1)   

    ### random walk in .int space
    ### normal Metropolis-Hastings
    #c_Theta = zeros(Float64,lME); c_Theta .= Theta    
    # iSPEs,V1s,V0s = readVint(intf,label_T1,label_T0)
    # ln_int = length(iSPEs)+length(V1s)+length(V0s)
    # Vint = zeros(Float64,ln_int)
    # iThetas = [ zeros(Float64,ln_int)];deleteat!(iThetas,1)   
    # nSPEs = copy(iSPEs); nV1s = copy(V1s); nV0s = copy(V0s)
    # evalVsnt!(Theta,iSPEs,V1s,V0s,idx_s_from_i,facs)
    # tx, llhs, Ens = intMCMC(itnum_MCMC,burnin,var_M,iThetas,Theta,c_Theta,Vint,
    #                         nSPEs,nV1s,nV0s,iSPEs,V1s,V0s,
    #                         idx_s_from_i,facs,
    #                         lnJ,tJNs,HNmats,AllTDs,Allinfos,
    #                         tMat,tildH,Erefs,errors,to)
    # plot_MCMC(iThetas,llhs,tJNs,Ens,Erefs)

    ### parallel tempering
    oSPEs,oV1s,oV0s = readVint(intf,label_T1,label_T0)
    iSPEs = [ copy(oSPEs) for i = 1:num_replica]
    V1s = [ copy(oV1s) for i = 1:num_replica]
    V0s = [ copy(oV0s) for i = 1:num_replica]
    Ts = [ 1.0 + 0.25*(i-1) for i=1:num_replica]
    #Ts = [ 1.0 *i for i=1:num_replica]

    ln_int = length(iSPEs[1])+length(V1s[1])+length(V0s[1])
    Vint = zeros(Float64,ln_int)
    iThetas = [ zeros(Float64,ln_int)];deleteat!(iThetas,1)   
    #AllThetas = [ [ zeros(Float64,ln_int)] for i = 1:num_replica]
    #for i=1:num_replica; deleteat!(AllThetas[i],1); end
    nSPEs = deepcopy(iSPEs); nV1s = deepcopy(V1s); nV0s = deepcopy(V0s)
    evalVsnt!(Theta,iSPEs[1],V1s[1],V0s[1],idx_s_from_i,facs)
    Thetas = [ copy(Theta) for i=1:num_replica]
    c_Thetas = deepcopy(Thetas)
    tx, llhs, Ens = intMCMC_PT(itnum_MCMC,burnin,var_M,Ts,
                               iThetas,Thetas,c_Thetas,Vint,
                               nSPEs,nV1s,nV0s,iSPEs,V1s,V0s,
                               oSPEs,oV1s,oV0s,
                               idx_s_from_i,facs,
                               lnJ,tJNs,HNmats,AllTDs,Allinfos,
                               tMat,tildH,Erefs,errors,to)
    
    Ns = length(iThetas)
    #for n =1:ln_int
    #    println("i= ", @sprintf("%4i", n),
    #            "mean: ",@sprintf("%10.4f ", mean([iThetas[i][n] for i =1:Ns])),
    #            "std: ", @sprintf("%10.4f ", std([iThetas[i][n] for i =1:Ns])))
    #end
    println(tx)    
    #plot_MCMC_PT(iThetas,Ts,llhs,tJNs,Ens,Erefs)   
    write_history(iThetas,Ts,llhs,tJNs,Ens)   
    show(to, allocations = true, compact = false)
    println("")
    return nothing
end

function write_history(iThetas,Ts,llhs,tJNs,Ens)
    try
        mkdir("history")
    catch
        nothing
    end
    nr = length(Ts)
    Ns = length(llhs[1])
    num_states = sum([tmp[2] for tmp in tJNs])
    
    oupf="./history/PTlog_Theta_lowestT.dat"
    io = open(oupf,"w")
    write(io,nr)
    write(io,Ns)
    write(io,length(iThetas[1]))
    for i=1:Ns
        write(io,iThetas[i])
    end
    close(io)
    oupf="./history/PTlog_llhs.dat"
    io = open(oupf,"w")
    write(io,nr)
    write(io,Ns)
    write(io,Ts)
    for i=1:nr
        write(io,llhs[i])
    end
    close(io)
    
    oupf="./history/PTlog_Ens.dat"
    io = open(oupf,"w")
    write(io,nr)
    write(io,Ns)
    for i=1:num_states
        write(io,Ens[i])
    end
    close(io)
end

function plot_from_history(Hs,target_nuc,tJNs,Erefs,Eexact,intf;path="./history/")
    sntf = Hs[1]
    Anum = parse(Int64, match(reg,target_nuc).match)
    tmp = readsnt(sntf,Anum)
    p_sps = tmp[8]; n_sps = tmp[9]
    label_T1,label_T0 = make_int(p_sps,n_sps)
    SPEs,V1s,V0s = readVint(intf,label_T1,label_T0)
    ln_int = length(SPEs) +length(V1s)+length(V0s)
    Vint = zeros(Float64,ln_int)
    evalVint!(Vint,SPEs,V1s,V0s)
    
    num_states = sum([tmp[2] for tmp in tJNs])   
    inpf= path*"PTlog_Theta_lowestT.dat"
    io = open(inpf,"r")
    nr = read(io,Int64)
    Ns = read(io,Int64)
    Np = read(io,Int64)    
    Thetas = [ Float64[1.0] ] ;deleteat!(Thetas,1)
    for i=1:Ns
        push!(Thetas,[read(io,Float64) for i=1:Np])
    end
    close(io)   
    inpf= path*"PTlog_llhs.dat"
    io=open(inpf,"r")
    nr = read(io,Int64)
    Ns = read(io,Int64)
    Ts = [read(io,Float64) for i=1:nr]   
    llhs = [ Float64[1.0] ];deleteat!(llhs,1)
    for i=1:nr
        push!(llhs,[read(io,Float64) for i=1:Ns])
    end
    close(io)

    AllEns = []
    paths = ["./history_run1/","./history_run2/","./history_run3/"]
    for path in paths
        inpf= path*"PTlog_Ens.dat"
        io = open(inpf,"r")
        nr = read(io,Int64)
        Ns = read(io,Int64)
        Ens = [ Float64[1.0] ];deleteat!(Ens,1)
        for i=1:num_states
            push!(Ens,[read(io,Float64) for i = 1:Ns])
        end
        close(io)
        push!(AllEns,copy(Ens))
    end
    #Thetas,llhs,Ens
    #println("nr $nr Ns $Ns Np $Np Numstates $num_states Ts $Ts")
    plot_MCMC_PT(Thetas,Ts,llhs,tJNs,AllEns,Erefs,Eexact,Vint)
end

function intMCMC(itnum_MCMC,burnin,var_M,iThetas,Theta,c_Theta,Vint,Vopt,
                 nSPEs,nV1s,nV0s,SPEs,V1s,V0s,
                 idx_s_from_i,facs,
                 lnJ,tJNs,HNmats,AllTDs,Allinfos,                 
                 tMat,tildH,Erefs,errors,to)
    llh = -1.e+10
    nllh = -1.e+10
    Accept = true
    Acchit = 0
    llhs = Float64[ ]
    tnum_states = sum([tJNs[i][2] for i =1:lnJ])
    Ens = [ Float64[] for i=1:tnum_states]
    offs = [0]
    for i = 2:lnJ
        push!(offs,sum([tJNs[j][2] for j =1:i-1]))
    end    
    for itM = 1:itnum_MCMC        
        #println("it=$itM")
        @timeit to "int." begin
            proposal_Vint!(nSPEs,nV1s,nV0s,SPEs,V1s,V0s,var_M)
            evalVsnt!(c_Theta,nSPEs,nV1s,nV0s,idx_s_from_i,facs)
        end
        nllh = 0.0
        @timeit to "EC" @inbounds for jidx = 1:lnJ
            @timeit to "Hmats" begin
                num_ev_target = tJNs[jidx][2]
                Hmat = HNmats[jidx][1]
                Linv = HNmats[jidx][2]
                TDs = AllTDs[jidx]
                infos = Allinfos[jidx]
                @inbounds @threads for idx = 1:length(TDs)
                    TD = TDs[idx]
                    i,j,Nij = infos[idx]
                    tmp = dot(TD,c_Theta)
                    Hmat[i,j] = tmp; Hmat[j,i] = tmp
                end
            end
            BLAS.gemm!('N','N',1.0,Linv,Hmat,0.0,tMat)
            BLAS.gemm!('T','N',1.0,Linv,tMat,0.0,tildH)

            @timeit to "Arpack" vals = real.(Arpack.eigs(tildH,nev=num_ev_target,which=:SR,
                                                         tol=1.e-8, maxiter=300))[1]
            nllh += L2_llh(vals,Erefs[jidx],errors[jidx])
            if itM > burnin
                @timeit to "push" for i=1:num_ev_target
                    push!(Ens[offs[jidx]+i],vals[i])
                end
            end
        end
        if nllh > llh
            Accept = true
        elseif nllh-llh > log(rand())
            Accept = true
        else
            Accept = false
        end
        if Accept
            llh = nllh
            Theta .= c_Theta
            SPEs .= nSPEs
            V1s .= nV1s
            V0s .= nV0s
            Acchit += 1
            evalVint!(Vint,SPEs,V1s,V0s)
        end
        @timeit to "int." if itM > burnin 
            push!(iThetas,copy(Vint))
            push!(llhs,llh)
        end
    end
    tx = "Accept rate" * @sprintf("%5.1f ",100*(Acchit-burnin)/(itnum_MCMC-burnin))
    return tx, llhs, Ens
end

function singleH(TD,i,j,Nij,c_Theta, Hmat)
    tmp = dot(TD,c_Theta)
    Hmat[i,j]=tmp;Hmat[j,i]=tmp
    return nothing
end

function intMCMC_PT(itnum_MCMC,burnin,var_M,Ts,
                    iThetas,Thetas,c_Thetas,Vint,
                    RnSPEs,RnV1s,RnV0s,RSPEs,RV1s,RV0s,
                    oSPEs,oV1s,oV0s,
                    idx_s_from_i,facs,
                    lnJ,tJNs,HNmats,AllTDs,Allinfos,                 
                    tMat,tildH,Erefs,errors,to)
    nr = length(Ts) # num_replica
    llhs = [-1.e+30 for i =1:nr]
    nllhs = [-1.e+10 for i =1:nr]
    Accept = false; TAcchit=0
    Acchits = zeros(Int64,nr)
    xi = copy(Thetas[1])
    xj = copy(Thetas[1])
    tnum_states = sum([tJNs[i][2] for i =1:lnJ])
    Ens = [ Float64[] for i=1:tnum_states]
    Evals = deepcopy(Erefs)
    illhs = [ Float64[] for i=1:nr] # llh history for lowest temperature
    offs = [0];for i = 2:lnJ;push!(offs,sum([tJNs[j][2] for j =1:i-1]));end

    # for destructive_eigen my_eigvals
    Dim = size(tMat)[1]
    v0 = [randn() for i=1:Dim]          
          
    vks = [ zeros(Float64,Dim) for i = 1:150 ]
    for i=1:Dim; vks[1][i] = randn(); end
    tmp = dot(vks[1],vks[1])
    vks[1] .*= 1.0/sqrt(tmp)
    num_history = 2
    en_s = [ [ zeros(Float64,tmp[2]) for i =1:num_history] for tmp in tJNs ]
    ####
    
    jobs = [ [0,0] ]; deleteat!(jobs,1)
    for jidx = 1:lnJ
        info = Allinfos[jidx]
        for idx = 1:length(AllTDs[jidx])
            i,j,Nij = info[idx]
            push!(jobs,[jidx,idx,i,j,Nij])
        end
    end
    
    for itM = 1:itnum_MCMC
        if itM % div(itnum_MCMC,10) ==0; println("it ",itM);end
        @inbounds for ridx = 1:nr
            nSPEs = RnSPEs[ridx]; nV1s = RnV1s[ridx]; nV0s = RnV0s[ridx]
            SPEs = RSPEs[ridx]  ;  V1s = RV1s[ridx] ;  V0s = RV0s[ridx]
            Theta = Thetas[ridx]; c_Theta = c_Thetas[ridx]                
            proposal_Vint!(nSPEs,nV1s,nV0s,SPEs,V1s,V0s,var_M*(ridx^(1/3)))
            evalVsnt!(c_Theta,nSPEs,nV1s,nV0s,idx_s_from_i,facs)
            nllhs[ridx] = L2norm(nSPEs,nV1s,nV0s,oSPEs,oV1s,oV0s)
            @timeit to "Hmats" @inbounds @threads for job in jobs
                jidx,idx,mi,mj,Nij = job
                Hmat = @views HNmats[jidx][1]
                Linv = @views HNmats[jidx][2]
                TD = @views AllTDs[jidx][idx]
                info = @views Allinfos[jidx][idx]
                singleH(TD,mi,mj,Nij,c_Theta,Hmat)
            end
            @timeit to "gem&diag" for jidx = 1:lnJ                
                num_ev_target = tJNs[jidx][2]
                Hmat = @views HNmats[jidx][1]
                Linv = @views HNmats[jidx][2]
                Eval = @views Evals[jidx]
                BLAS.gemm!('N','N',1.0,Linv,Hmat,0.0,tMat)
                BLAS.gemm!('T','N',1.0,Linv,tMat,0.0,tildH)                

                #@timeit to "my_eigval" my_eigvals!(tildH,tMat,Dim,vks,en_s[jidx],
                #                                   num_ev_target)                
                #Evals[jidx] .= en_s[jidx][1]               
                #println("En my ", Evals[jidx]) 

                #@timeit to "Arpack" Evals[jidx] .= real.(
                #    Arpack.eigs(tildH,nev=num_ev_target,
                
                @timeit to "Arpack" begin
                    Eval .= Arpack.eigs(tildH,nev=num_ev_target,
                                        which=:SR,tol=1.e-6,maxiter=150)[1]
                end
                #println("En Arpack ", Evals[jidx], "\n")
                nllhs[ridx] += L2_llh(Eval,Erefs[jidx],errors[jidx];T=Ts[ridx])
                if ridx == 1
                    if itM > burnin 
                        for i=1:num_ev_target
                            push!(Ens[offs[jidx]+i],Eval[i])
                        end
                    end
                end
            end
            Accept = false
            if nllhs[ridx] > llhs[ridx]
                Accept =true
            elseif nllhs[ridx]-llhs[ridx] > log(rand())
                Accept = true
            end
            if Accept
                llhs[ridx] = nllhs[ridx]
                Theta .= c_Theta
                SPEs .= nSPEs; V1s .= nV1s;V0s .= nV0s
                evalVint!(Vint,SPEs,V1s,V0s)
                Acchits[ridx] += 1                
            end
            if itM > burnin
                if ridx == 1; push!(iThetas,copy(Vint)); end
                push!(illhs[ridx],llhs[ridx])
                continue
            end
            if itM == burnin;Acchits .= 0;TAcchit=0;end
        end
        ### exchange
        if itM <= burnin || nr == 1; continue ;end
        ri = itM % 2 + 1
        for r1 in ri:2:nr-1
            r2 = r1 + 1
            Ti = Ts[r1]; Tj = Ts[r2]
            xi .= Thetas[r1]; xj .= Thetas[r2]
            llh_i = llhs[r1]; llh_j = llhs[r2]
            logp = (Ti-Tj) * (llh_i/Tj  - llh_j/Ti)
            TAccept = false
            if logp > 0.0
                TAccept=true
            elseif log(rand()) < logp
                TAccept=true
            end
            if TAccept && itM > burnin
                Thetas[r1] .= xj; Thetas[r2] .= xi
                illhs[r1][end] = llh_j
                illhs[r2][end] = llh_i
                TAcchit += 1
                if r1 ==1
                    evalVint!(Vint,RSPEs[r2],RV1s[r2],RV1s[r2])
                    iThetas[end] .= Vint
                end
            end
        end
    end
    rates = (100/(itnum_MCMC-burnin)) .* Acchits
    tx = "Accept rates "
    for i = 1:nr
        tx *= "  "*@sprintf("%5.1f ", rates[i])
    end
    tx *= "\t exchange rate " *@sprintf("%5.1f ", 100*TAcchit/((itnum_MCMC-burnin)*nr))
    return tx, illhs, Ens
end

function my_eigvals!(A,T,Dim,vks,en,num_ev_target;TF=[false],tol=1.e-6)
    T .= 0.0
    for i=1:length(en); en[i] .= 1.e+10;end
    for i=1:Dim; vks[1][i] = randn(); end
    tmp = 1.0 / sqrt(dot(vks[1],vks[1]))
    vks[1] .*= tmp

    for it = 1:length(vks)-1
        v  = vks[it]; Hv = vks[it+1]
        mul!(Hv,A,v)
        talpha = dot(v,Hv)
        T[it,it] = talpha        
        #diagonalize_T!(it,num_ev_target,T,en,2,TF,tol)
        #print_vec("HNLanczos @$it ",en[1])
        if TF[1];break;end
        axpy!(-talpha,v,Hv)
        svks = @views vks[1:it-1]
        ReORTH(it,Hv,svks)
        tbeta = sqrt(dot(Hv,Hv))
        tmp = 1.0/tbeta
        Hv .*= tmp
        T[it+1,it] = tbeta; T[it,it+1] = tbeta
    end
    return nothing
end

function proposal_Theta!(c_Theta,Theta,ln,var_M)
    @inbounds for i=1:ln
        c_Theta[i] = Theta[i] + var_M .* randn()
    end
    return nothing
end
function L2norm(nSPEs,nV1s,nV0s,oSPEs,oV1s,oV0s;lam=1e1)
    s = 0.0
    for i =1:length(nSPEs); s += (nSPEs[i]-oSPEs[i])^2; end    
    for i =1:length(nV1s); s += (nV1s[i]-oV1s[i])^2; end
    for i =1:length(nV0s); s += (nV0s[i]-oV0s[i])^2; end
    return - 0.5 * s * lam
end
function L2_llh(Eval,Eref,err;T=1.0,sigma_sm=0.126) 
    s=0.0
    @inbounds for i = 1:length(Eval)
        #s -=  (Eval[i] - Eref[i])^2 / (sigma_sm^2 + err[i]^2) #i-dependent
        #s -=  (Eval[i] - Eref[i])^2 # i-independent err=1.0
        s -=  (Eval[i] - Eref[i])^2 / (1.0+err[i]^2) # i-independent err=1.0
    end
    return 0.5 * s / (T*length(Eval))
end

function read_exact(target_nuc,targetJ,sntf,lines)
    hit = 0
    for (i,line) in enumerate(lines)
        if occursin("$target_nuc",line)
            if parse(Int,split(lines[i+1],"targetJ")[end]) == targetJ
                hit = 1
                continue
            end
        end
        if hit ==0;continue;end
        csntf = split(sntf,"/")[end]
        if occursin(csntf,line)
            oJs = split(split(split(lines[i+1],"[")[2],"]")[1],",")
            oEs = split(lines[i+2])[2:end]
            Es=[]
            for (idx,J) in enumerate(oJs)
                tJ = Int(2*parse(Float64,J))
                if tJ != targetJ; continue;end
                push!(Es,parse(Float64,oEs[idx]))
            end
            return Es
        end
    end
    return []
end

function plot_EC_scatter(target_nuc,Hs,sumV,tJNs,Dims,exlines)
    js = [ tJNs[i][1] for i=1:length(tJNs)]
    xs = [ Float64[] for i=1:length(js)]
    ys = [ Float64[] for i=1:length(js)]
    minmax=Float64[1.e+10,-1.e+10]
    yerrs = Float64[]
    for tmp in sumV
        sntf,tJ,EnsEC = tmp
        Ens = read_exact(target_nuc,tJ,sntf,exlines)
        tl = Float64[]        
        idx = 0
        for k = 1:length(js)
            if js[k]==tJ
                idx = k;break
            end
        end        
        try
            Ens[1]; EnsEC[1]
        catch
            continue
        end
        for i=1:1 # only the yrast state           
            #if EnsEC[i] < Ens[i]
            #    println("variational err ",
            #            " EC ",EnsEC[i], " Exact ",Ens)
            #end
            push!(xs[idx],Ens[i])
            push!(ys[idx],EnsEC[i])
            push!(tl,Ens[i])
            push!(tl,EnsEC[i])
            push!(yerrs,(EnsEC[i] - Ens[i])/abs(Ens[i]))
        end
        if minimum(tl) < minmax[1]
            minmax[1] = minimum(tl)
        end
        if maximum(tl) > minmax[2]
            minmax[2] = maximum(tl)
        end
    end
    if length(Hs) > 1; println(" typical error ", mean(yerrs));end
    ##scatter
    plt=pyimport("matplotlib.pyplot")
    cm=pyimport("matplotlib.cm")
    patches=pyimport("matplotlib.patches")
    fig = plt.figure(figsize=(6,4))
    ax = fig.add_subplot(111);   #axB = fig.add_subplot(212)
    ax.set_xlabel("Exact (MeV)")
    ax.set_ylabel("EC estimate (MeV)")
    cols = ["red","blue","darkgreen","orange","purple","magenta","cyan"]
    tms = ["o","^",",","d","p","h","8","v","*","1"]
    ax.plot([minmax[1]-5,minmax[2]+5],[minmax[1]-5,minmax[2]+5],
            linestyle="dotted",color="gray",alpha=0.6)
    for jidx = 1:length(js)
        nth = 1 #yrast only
        tl = nothing
        if nth == 1
            tl = "J="*latexstring(string(J2_to_J(js[jidx]))*"_{"*string(nth)*"}")
        end
        ax.scatter(xs[jidx],ys[jidx],marker=tms[jidx],s=20,lw=0.8,
                   zorder=150,color=cols[jidx],facecolor="none",label=tl,alpha=0.5)
    end
    Anum = parse(Int64, match(reg,target_nuc).match)
    nuc_latex = latexstring("{}^{$Anum}")*split(target_nuc,string(Anum))[1]
    
    tl = " $nuc_latex ("*latexstring("N_s")*"="*string(Dims[1])*")"
    ax.text(0.05, 0.9, tl,fontsize=15, transform=ax.transAxes)
    ax.legend(loc="lower right",fontsize=12)
    plt.savefig("./pic/scatter_ECestimates_"*target_nuc*".pdf",bbox_iinches="tight",pad_inches=0)
    plt.close()
end


function ngauss(x,mu,varV)
    return 1.0/(sqrt(2*pi*varV)) * exp( - 0.5* (x-mu)^2 / varV )
end

function plot_MCMC_PT(iThetas,Ts,llhs,tJNs,AllEns,Erefs,Eexact,Vint)
    plt=pyimport("matplotlib.pyplot")
    cm=pyimport("matplotlib.cm")
    patches=pyimport("matplotlib.patches")
    ## loglikelihood
    nr = length(Ts)
    fig = plt.figure(figsize=(10,6))
    ax = fig.add_subplot(111)
    ax.set_xlabel("MC step")
    ax.set_ylabel("loglikelihood")
    #ax.set_ylim(-50,0)
    for i = 1:nr
        ax.plot(llhs[i].*Ts[i],alpha=0.6,lw=0.5,
                label="T="*string(Ts[i]),rasterized=true,zorder=-i)
    end
    #ax.legend()
    plt.savefig("./pic/PT_MCMC_llhs.pdf",bbox_iinches="tight",pad_inches=0.0)
    plt.close()
    
    ## theta plot
    fig = plt.figure(figsize=(10,4))
    ax = fig.add_subplot(111)
    ax.set_xlabel("MC step")
    ax.set_ylabel("SPEs&TBMEs (MeV, .int fmt.)")
    Ns = length(iThetas)
    for i = 1:length(iThetas[1])        
        ax.plot([ iThetas[n][i] for n=1:Ns],alpha=0.5,lw=0.1,rasterized=true)
    end
    plt.savefig("./pic/PT_MCMC_thetas.pdf",bbox_iinches="tight",pad_inches=0.0)
    plt.close()
    fig = plt.figure(figsize=(14,6))
    axs = [ fig.add_subplot(6,11,i) for i =1:66]
    Ns = length(iThetas)
    tbin = 20
    varV= 1.0/10.0
    for i = 1:length(iThetas[1])
        mu = Vint[i];  sigma = sqrt(varV)
        tbin = round(mu-5*sigma):0.1:round(mu+5*sigma)
        axs[i].hist([iThetas[n][i] for n=1:Ns],bins=tbin,density=true,color="darkgreen",alpha=0.8)
        xr = Vint[i]-5*sqrt(varV):0.01:Vint[i]+5*sqrt(varV)
        yr = ngauss.(xr,Vint[i],varV)
        axs[i].plot(xr,yr,linestyle="dashed",color="k")
        axs[i].tick_params(labelleft=false,left=false)
        axs[i].scatter(Vint[i],0,marker="o",color="k",zorder=300)
        axs[i].text(0.03, 0.77,latexstring("n=")*@sprintf("%2i",i),
                    transform=axs[i].transAxes)
    end

    fig.subplots_adjust(wspace=0.0,hspace=0.5)
    plt.savefig("./pic/PT_MCMC_thetas_hist.pdf",bbox_iinches="tight",pad_inches=0.0)
    plt.close()
    
    ### Energy dist. (histgram)
    ln = sum([tmp[2] for tmp in tJNs])
    fig = plt.figure(figsize=(14,6))
    axs = [ fig.add_subplot(3,4,i) for i =1:ln]
    mu = 0.0; sigma=0.0
    cols=["red","blue","green"]
    for i = 1:ln
        for n = 1:length(AllEns)
            if n==1
                mu = mean(AllEns[n][i]);sigma = std(AllEns[n][i])
                axs[i].set_xlim(mu-5*sigma, mu+5*sigma)
                axs[i].tick_params(labelleft=false,left=false)
                #println("axs[i].get_xscale ",axs[i].get_xscale)
                tbin = round(mu-5*sigma):0.5:round(mu+5*sigma)
            end
            #println("para mu: ",mean(AllEns[n][i]),
            #        " std: ", std(AllEns[n][i])," size:" , length(AllEns[n][i]))            
            axs[i].hist(AllEns[n][i],bins=tbin,density=true,edgecolor=cols[n],facecolor="None",
                        alpha=0.6,zorder=300+10*n)
        end
        is, ie = axs[i].get_ylim()
        xstep = 3.0
        axs[i].xaxis.set_ticks(tbin[1]+xstep:xstep:tbin[end]-xstep)
    end
    count = 0
    for (jidx,tmp) in enumerate(tJNs)
        for i=1:tmp[2]
            count += 1
            cJ = tJNs[jidx][1]
            if cJ % 2 == 0
                cJ = latexstring( string(div(cJ,2))*"_{"*string(i)*"}")
            else                
                cJ = latexstring( string(cJ)*"/2_{"*string(i)*"}")       
            end
            tl="";if count==1; tl = "Exp.";end
            tlO="";if count==2; tlO = "USDB";end
            axs[count].scatter(Erefs[jidx][i],0.0,marker="o",color="k",label=tl,zorder=500)
            axs[count].scatter(Eexact[jidx][i],0.0,marker="x",color="r",label=tlO,zorder=600)
            axs[count].text(0.1, 0.8,latexstring("J=")*cJ,
                            transform=axs[count].transAxes)
        end
    end
    axs[1].legend(loc="upper right")
    axs[2].legend(loc="upper right")
    fig.subplots_adjust(wspace=0.0)
    plt.savefig("./pic/PT_MCMC_hist.pdf",bbox_iinches="tight",pad_inches=0.0)
    plt.close()
    # count = 0
    # for (jidx,tmp) in enumerate(tJNs)
    #     for i=1:tmp[2]
    #         count += 1
    #         tx = "E(EC,MCMC) mean: " * @sprintf("%10.4f ", mean(Ens[1][count]))
    #         tx *= " std: "*@sprintf("%10.4f ", std(Ens[count]))
    #         tx *= " Eexp. "*@sprintf("%10.4f ",Erefs[jidx][i])
    #         println(tx)
    #     end
    # end
    return nothing
end

struct Tridx
    i::Int
    j::Int
    Nth::Int
end

function calcOBTD(OBTDs,idxs,p_sps,n_sps,
                  mstates_p,mstates_n,
                  tdims,
                  jocc_p,jocc_n,
                  pbits,nbits,
                  wfs)
    lps=length(p_sps);lns=length(n_sps);lblock = length(nbits)
    #@inbounds @qthreads for k=1:length(idxs)
    @inbounds for k=1:length(idxs)
        #for k=1:length(idxs)       
        tmp = idxs[k]
        i = tmp.i; j=tmp.j; Nij=tmp.Nth
        w_i=wfs[i]; w_j = wfs[j]
        pOBTD = OBTDs[Nij][1]
        nOBTD = OBTDs[Nij][2]
        @inbounds for bi = 1:lblock        
            idim = tdims[bi]
            jocc_p_bi = jocc_p[bi]
            jocc_n_bi = jocc_n[bi]
            pbit = pbits[bi]
            nbit=nbits[bi]
            l_Nn = length(nbit)
            offset = idim - l_Nn
            @inbounds for (pidx,pocc) in enumerate(jocc_p_bi)
                tM = offset + pidx*l_Nn
                @inbounds for (nidx,nocc) in enumerate(jocc_n_bi)
                    Mi = tM + nidx
                    wfprod =  w_i[Mi] .* w_j[Mi]
                    @inbounds for i = 1:lps
                        pOBTD[i] += wfprod .* pocc[i]
                    end                    
                    @inbounds for i = 1:lns
                        nOBTD[i] += wfprod .* nocc[i]
                    end
                    #pOBTD .+= wfprod .* pocc;nOBTD .+= wfprod .* nocc
                    #axpy!(wfprod,pocc,pOBTD);axpy!(wfprod,nocc,nOBTD)     
                end
            end
        end
    end
    if false
        println("OBTDs: ")
        for idx in idxs            
            i = idx.i; j=idx.j; Nij=idx.Nth
            if i == j
                println("<SPE>_[$i] Nij $Nij")
                println(" proton: ",OBTDs[Nij][1])
                println("neutron: ",OBTDs[Nij][2]);println("")
            end                
        end
    end   
    return nothing
end

struct T1ifc
    nth::Int64
    bi::Int64
    i::Int64
    f::Int64
    fac::Float64
end
struct T0ifc
    i::Int64
    f::Int64
    fac::Float64
end

function prepTBTD(tJ,idxs,p_sps,n_sps,
                  mstates_p::Array{Array{Int64,1},1},
                  mstates_n::Array{Array{Int64,1},1},
                  pbits,nbits,
                  labels,TBMEs,
                  wfs::Array{Array{Float64,1}},
                  tdims,Mps)
    lblock=length(pbits)
    lp = length(mstates_p); ln = length(mstates_n)    
    mstates = [mstates_p,mstates_n]
    sps = [p_sps,n_sps]
    loffs = [ 0, length(p_sps)]
    TBTD1 = [ T1ifc[ ] , T1ifc[] ]
    n = nthreads()
    #@threads
    for vrank =1:2 #pp:1, nn:2
        loff = loffs[vrank]
        T1info = TBTD1[vrank]
        vecs= [ [ [ false for i = 1:lp] for j=1:2],
                [ [ false for i = 1:ln] for j=1:2]]
        @inbounds for (i,ME) in enumerate(TBMEs[vrank])
            a,b,c,d,totJ = labels[vrank][i]
            J2  = 2*totJ
            ja,ma_s,ma_idxs = possible_mz(sps[vrank][a-loff],mstates[vrank]) 
            jb,mb_s,mb_idxs = possible_mz(sps[vrank][b-loff],mstates[vrank]) 
            jc,mc_s,mc_idxs = possible_mz(sps[vrank][c-loff],mstates[vrank]) 
            jd,md_s,md_idxs = possible_mz(sps[vrank][d-loff],mstates[vrank])
            @inbounds for (ic,mc) in enumerate(mc_s)
                @inbounds for (id,md) in enumerate(md_s)
                    if c == d && mc >= md; continue;end
                    if abs(mc + md) > J2; continue;end
                    M_ani = mc + md
                    initialize_tvec!(vecs[vrank][1]); vecs[vrank][1][mc_idxs[ic]] = true
                    bit_c = bitarr_to_int(vecs[vrank][1])
                    initialize_tvec!(vecs[vrank][1]); vecs[vrank][1][md_idxs[id]] = true
                    bit_d = bitarr_to_int(vecs[vrank][1])
                    @inbounds for (ia,ma) in enumerate(ma_s)
                        @inbounds for (ib,mb) in enumerate(mb_s)
                            if a == b && ma>=mb; continue;end
                            if ma + mb != M_ani;continue;end
                            initialize_tvec!(vecs[vrank][2]);vecs[vrank][2][ma_idxs[ia]] = true
                            bit_a = bitarr_to_int(vecs[vrank][2])
                            initialize_tvec!(vecs[vrank][2]);vecs[vrank][2][mb_idxs[ib]] = true
                            bit_b = bitarr_to_int(vecs[vrank][2])
                            CG1 = clebschgordan(Float64,ja//2,ma//2,jb//2,mb//2,J2//2,M_ani//2)
                            CG2 = clebschgordan(Float64,jc//2,mc//2,jd//2,md//2,J2//2,M_ani//2)
                            bitlist = bit2b(bit_a,bit_b,bit_c,bit_d)
                            coeff = sqrt( (1.0+deltaf(a,b)) *(1.0+deltaf(c,d)) ) * CG1 * CG2                            
                            if vrank==1 # pp
                                @inbounds for bi = 1:lblock
                                    TF=[true]; ret=[1,-1]; ridx=[-1,-1,-1]
                                    idim = tdims[bi]
                                    l_Nn = length(nbits[bi])
                                    @inbounds for (Npi,Phi) in enumerate(pbits[bi])
                                        TF_connectable(Phi,bitlist,TF)
                                        if TF[1]==false; continue;end
                                        calc_phase!(Phi,bitlist,lp,ret)
                                        bisearch!(pbits[bi],ret[2],ridx)
                                        Npf = ridx[1]
                                        fac = coeff .*ret[1] 
                                        push!(T1info,T1ifc(i,bi,Npi,Npf,fac))
                                    end
                                end                                
                            else ## nn
                                @inbounds for bi = 1:lblock
                                    TF=[true]; ret=[1,-1]; ridx=[-1,-1,-1]
                                    idim = tdims[bi]
                                    l_Nn = length(nbits[bi])
                                    @inbounds for (Nni,Phi) in enumerate(nbits[bi])
                                        TF_connectable(Phi,bitlist,TF)
                                        if TF[1]==false; continue;end
                                        calc_phase!(Phi,bitlist,ln,ret)
                                        bisearch!(nbits[bi],ret[2],ridx)
                                        Nnf = ridx[1]
                                        fac = coeff .*ret[1] 
                                        push!(T1info,T1ifc(i,bi,Nni,Nnf,fac))
                                    end
                                end
                            end                               
                        end
                    end
                end                
            end
        end
    end

    ## pn
    vrank=3
    loff = loffs[2] # for neutron
    delMs = Int64[]
    for j = 1:lp
        for i = 1:lp
            push!(delMs,mstates_p[i][5]-mstates_p[j][5])
        end
    end
    unique!(delMs);sort!(delMs,rev=true)
    maxDeltaM = maximum(delMs)    
    bfs = [ Int64[ ] for i = 1:lblock]
    bis = [ Int64[ ] for i = 1:lblock]
    for i = 1:lblock
        for j = i:lblock
            if Mps[j] - Mps[i] in delMs
                push!(bfs[i],j); push!(bis[j],i)
            end
        end
    end

    vec_ani_p = [false for i = 1:lp];vec_cre_p = [false for i = 1:lp]
    vec_ani_n = [false for i = 1:ln];vec_cre_n = [false for i = 1:ln]
    retM = [0,0,0]
    bif_idxs=[ [0,0,0] ];deleteat!(bif_idxs,1)
    for i = 1:lblock
        for j=i:lblock
            push!(bif_idxs,[i,j,idx_from_ij(i,j,lblock)])
        end
    end
    pjumps = [ [ [ [T0ifc(0,0,0.0)] ] ] for i=1:length(TBMEs[vrank]) ]
    njumps = [ [ [ [T0ifc(0,0,0.0)] ] ] for i=1:length(TBMEs[vrank]) ]
    for i = 1:length(TBMEs[vrank])
        deleteat!(pjumps[i],1)
        deleteat!(njumps[i],1)
    end
    for (nth,ME) in enumerate(TBMEs[vrank])
        a,b,c,d,totJ,vidx = labels[vrank][nth]
        J2  = 2*totJ
        ja,ma_s,ma_idxs = possible_mz(p_sps[a],mstates_p) 
        jc,mc_s,mc_idxs = possible_mz(p_sps[c],mstates_p)
        jb,mb_s,mb_idxs = possible_mz(n_sps[b-loff],mstates_n)         
        jd,md_s,md_idxs = possible_mz(n_sps[d-loff],mstates_n)
        for (ic,mc) in enumerate(mc_s)
            initialize_tvec!(vec_ani_p); vec_ani_p[mc_idxs[ic]] = true
            bit_c = bitarr_to_int(vec_ani_p)
            for (ia,ma) in enumerate(ma_s)
                initialize_tvec!(vec_cre_p); vec_cre_p[ma_idxs[ia]] = true
                bit_a = bitarr_to_int(vec_cre_p)
                Mp = ma - mc
                bisearch!(delMs,Mp,retM);idx = retM[1]
                for (id,md) in enumerate(md_s)
                    if abs(mc + md) > J2; continue;end
                    initialize_tvec!(vec_ani_n); vec_ani_n[md_idxs[id]] = true
                    bit_d = bitarr_to_int(vec_ani_n)
                    for (ib,mb) in enumerate(mb_s) 
                        if mb - md + Mp != 0; continue;end
                        initialize_tvec!(vec_cre_n); vec_cre_n[mb_idxs[ib]] = true
                        bit_b = bitarr_to_int(vec_cre_n)
                        CG1 = clebschgordan(Float64,ja//2,ma//2,jb//2,mb//2,
                                            J2//2,(ma+mb)//2)
                        CG2 = clebschgordan(Float64,jc//2,mc//2,jd//2,md//2,
                                            J2//2,(mc+md)//2)
                        fac = CG1*CG2 #* sqrt(tJ+1) /sqrt(J2+1.0)  
                        pjump = [ T0ifc[] for i=1:length(bif_idxs)]
                        njump = [ T0ifc[] for i=1:length(bif_idxs)]                        
                        TF=[true]; ret_p=[0,0];ret_n=[0,0]
                        ridx=[0,0,0]

                        for bi = 1:lblock
                            idim = tdims[bi]
                            l_Nn_i = length(nbits[bi])
                            for bfidx = 1:length(bfs[bi]) 
                                bf = bfs[bi][bfidx]
                                deltaM = Mps[bf] - Mps[bi]
                                tfac = fac * ifelse(bi==bf,0.5,1.0)
                                if (deltaM in delMs) == false;continue;end
                                fdim = tdims[bf]
                                l_Nn_f = length(nbits[bf])
                                bif_idx = idx_from_ij(bi,bf,lblock)
                                for (Npi,pPhi) in enumerate(pbits[bi])
                                    TF_connectable_1(pPhi,bit_a,bit_c,TF)
                                    if TF[1]==false; continue;end
                                    calc_phase_1!(pPhi,bit_a,bit_c,ret_p)
                                    bisearch!(pbits[bf],ret_p[2],ridx)
                                    if ridx[1] == 0;continue;end
                                    Npf = ridx[1]
                                    push!(pjump[bif_idx],
                                          T0ifc(Npi,Npf,tfac*ret_p[1]))
                                end
                                for (Nni,nPhi) in enumerate(nbits[bi])
                                    TF_connectable_1(nPhi,bit_b,bit_d,TF)
                                    if TF[1]==false; continue;end
                                    calc_phase_1!(nPhi,bit_b,bit_d,ret_n)
                                    bisearch!(nbits[bf],ret_n[2],ridx)
                                    if ridx[1] == 0;continue;end
                                    Nnf = ridx[1]
                                    push!(njump[bif_idx],
                                          T0ifc(Nni,Nnf,1.0*ret_n[1]))
                                end
                                # if ma==mb==mc==md==-1
                                #     if a==1 && b==3 && c==1 && ;
                                #        d==3 && totJ==1 && bi==bf
                                #         println("V($a$b$c$d$totJ) hit ",
                                #                 labels[vrank][nth],
                                #                 " ms $ma $mb $mc $md ",
                                #                 " bi $bi bf $bf bifidx $bif_idx")
                                #         println("pjump ", pjump[bif_idx],
                                #                 " njump ", njump[bif_idx])
                                #         println("")
                                #     end
                                # end
                            end
                        end
                        push!(pjumps[nth],pjump)
                        push!(njumps[nth],njump)
                    end
                end                
            end
        end
    end
    return TBTD1, pjumps,njumps,bif_idxs
end

function calcTBTD(TBTDs,
                  opTBTD1,pjumps,njumps,
                  pbits,nbits,tdims,wfs,idxs,bif_idxs,
                  olabels,oTBMEs,labels,to)
    @inbounds @qthreads for k=1:length(idxs)
        #@inbounds for k=1:length(idxs)
        tmp = idxs[k]
        i = tmp.i; j=tmp.j; Nij=tmp.Nth
        w_i=wfs[i]; w_j = wfs[j]
        coeffs = [ zeros(Float64,length(oTBMEs)) for vrank=1:3]
        #@timeit to "pp/nn" begin
        ## pp/nn
        vrank =1
        for tmp in opTBTD1[vrank] ## pp
            nth = tmp.nth; bi = tmp.bi
            nbit = nbits[bi]
            pbit = pbits[bi]
            lN = length(nbit)
            vidx = labels[vrank][nth][end]
            Npi = tmp.i;Npf = tmp.f;fac = tmp.fac
            tMi = tdims[bi]+ (Npi-1)*length(nbit)
            tMf = tdims[bi]+ (Npf-1)*length(nbit)                
            @inbounds for nidx = 1:length(nbit)
                Mi = tMi+nidx; Mf = tMf+nidx
                coeffs[vrank][vidx] += fac .* (w_i[Mf].*w_j[Mi])
            end
        end
        vrank =2
        for tmp in opTBTD1[vrank] ## nn
            nth = tmp.nth; bi = tmp.bi
            nbit = nbits[bi]
            pbit = pbits[bi]
            lN = length(nbit)
            vidx = labels[vrank][nth][end]
            Nni = tmp.i;Nnf = tmp.f;fac = tmp.fac            
            tMi = tdims[bi]+ Nni -lN 
            tMf = tdims[bi]+ Nnf -lN 
            @inbounds for pidx = 1:length(pbit)
                Mi = tMi + pidx*lN; Mf = tMf + pidx*lN
                coeffs[vrank][vidx] += fac .* (w_i[Mf].*w_j[Mi])
            end 
        end
        #@timeit to "pn" begin        
        ### pn
        vrank=3
        @inbounds for (nth,label) in enumerate(labels[vrank]) 
            vidx = label[end]
            #@inbounds @threads for nth = 1:length(labels[vrank]) 
            #    vidx = labels[vrank][nth][end]
            pjump = pjumps[nth]
            njump = njumps[nth]
            @inbounds for opidx = 1:length(pjump)
                tmp_pj = pjump[opidx]
                tmp_nj = njump[opidx]
                @inbounds for bidx = 1:length(bif_idxs)
                    bi,bf,dummy = bif_idxs[bidx]
                    idim = tdims[bi]; lNi = length(nbits[bi])
                    fdim = tdims[bf]; lNf = length(nbits[bf])
                    pjs = tmp_pj[bidx]
                    njs = tmp_nj[bidx]
                    @inbounds for pj in pjs
                        Npi = pj.i;Npf = pj.f; fac_p = pj.fac
                        tMi = idim + (Npi-1)*lNi
                        tMf = fdim + (Npf-1)*lNf
                        @inbounds for nj in tmp_nj[bidx]                                
                            Nni = nj.i; Nnf=nj.f; fac_n = nj.fac
                            Mi = tMi + Nni;  Mf = tMf + Nnf
                            coeffs[vrank][vidx] +=  fac_p .* fac_n .* (w_i[Mf].*w_j[Mi])
                            coeffs[vrank][vidx] +=  fac_p .* fac_n .* (w_i[Mi].*w_j[Mf]) 
                        end                   
                    end
                end
            end
        end
        for k = 1:3
            TBTDs[Nij] .+= coeffs[k]
        end
    end
    return nothing
end
