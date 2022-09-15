function caliblating_2n3nLECs_byHFMBPT(itnum,optimizer,MPIcomm,chiEFTobj,OPTobj,d9j,HOBs,nucs,HFdata,to,io;Operators=[""])    
    @assert optimizer != "" "optimizer should be BayesOpt/LHS/MCMC"
    if !MPIcomm        
        for it = 1:itnum
            add2n3n(chiEFTobj,to,it)
            @timeit to "Vtrans" dicts_tbme = TMtrans(chiEFTobj,HOBs,to;writesnt=false)
            print_vec("it = "*@sprintf("%8i",it),OPTobj.params,io)
            @timeit to "HF/HFMBPT" hf_main_mem(chiEFTobj,nucs,dicts_tbme,d9j,HOBs,HFdata,to;Operators=Operators)
            if optimizer=="BayesOpt"
                BO_HFMBPT(it,OPTobj,HFdata,to)
            elseif optimizer=="LHS"
                LHS_HFMBPT(it,OPTobj,HFdata,to)
            elseif optimizer=="MCMC"
                MCMC_HFMBPT(it,OPTobj,HFdata,to)
            else
                @error "unknown optimizer: $optimizer " 
                exit()
            end
            updateLECs_in_chiEFTobj!(chiEFTobj,OPTobj.targetLECs,OPTobj.params)
        end
    else
        mpi_hfmbpt(itnum,chiEFTobj,OPTobj,d9j,HOBs,nucs,HFdata,to,io;Operators=Operators)
    end
    return nothing
end

mutable struct MCMCobject
    dim::Int64    
    nstep::Int64
    burnin::Int64
    thining::Int64
    acchit::Int64
    targetLECs::Vector{String}
    params::Vector{Float64}
    params_ref::Vector{Float64}
    sigmas::Vector{Float64}
    cand::Vector{Float64}
    chain::Matrix{Float64}
    history::Vector{Vector{Float64}}
end
mutable struct MPIMCMCobject
    a::Float64
    dim::Int64    
    nstep::Int64
    burnin::Int64
    thining::Int64
    acchit::Int64
    targetLECs::Vector{String}
    params::Vector{Float64}
    params_ref::Vector{Float64}
    cand::Vector{Float64}
    chain::Array{Float64,3}
    history::Vector{Vector{Float64}}
    ens1::Vector{Int64}
    ens2::Vector{Int64}
end

struct LHSobject
    maxDim::Int64
    targetLECs::Vector{String}
    params::Vector{Float64}
    params_ref::Vector{Float64}
    pdomains::Vector{Tuple{Float64, Float64}}
    cand::Vector{Vector{Float64}}
    observed::Vector{Int64}
    unobserved::Vector{Int64}
    history::Vector{Vector{Float64}}
end

struct BOobject
    maxDim::Int64
    targetLECs::Vector{String}
    params::Vector{Float64}
    params_ref::Vector{Float64}
    pdomains::Vector{Tuple{Float64, Float64}}
    pKernel::Vector{Float64}
    Data::Vector{Vector{Float64}}
    cand::Vector{Vector{Float64}}
    observed::Vector{Int64}
    unobserved::Vector{Int64}
    history::Vector{Vector{Float64}}
    Ktt::Matrix{Float64}
    Ktinv::Matrix{Float64}
    Ktp::Matrix{Float64}
    L::Matrix{Float64}
    tMat::Matrix{Float64}
    yt::Vector{Float64}
    yscale::Vector{Float64}
    acquis::Vector{Float64}   
end

function myCholesky!(tmpA,ln,cLL)
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
    return nothing
end

function eval_HFMBPT(it,OPTobj,HFdata,varE,Lam;mcmc=false,io=stdout,debug=false)
    thist = OPTobj.history[it]
    params = OPTobj.params
    params_ref = OPTobj.params_ref
    tvec = params-params_ref
    logprior = -0.5*Lam*dot(tvec,tvec)
    if mcmc
        Lamt = params-params_ref
        tvec[1] *= 4.0
        logprior = -0.5*dot(Lamt,tvec)
    end    
    llh = 0.0
    for (n,tmp) in enumerate(HFdata)
        nuc = tmp.nuc
        data = tmp.data
        dtype = tmp.datatype
        A= nuc.A
        for (i,tdtype) in enumerate(dtype)
            vtho = data[i][1] *ifelse(tdtype=="E",1/A,1.0)
            vexp = data[i][2] *ifelse(tdtype=="E",1/A,1.0)
            tllh = 0.5 * (vtho-vexp)^2 / varE
            llh -= tllh
        end 
    end
    if debug; llh = 0.0;end
    logpost = logprior + llh
    thist[1] = logprior; thist[2] = llh; thist[3] = logpost
    if !mcmc
        println(io,"eval: ","logprior ",@sprintf("%9.2e",logprior),
        "  logllh  ",@sprintf("%9.2e",llh),"  logpost ",@sprintf("%9.2e",logpost))
    else
        println(io,"eval: ","logprior ",@sprintf("%9.2e",logprior),
                "  logllh  ",@sprintf("%9.2e",llh),"  logpost ",@sprintf("%9.2e",logpost),
                " Acc.Rate ",@sprintf("%6.2f",100*OPTobj.acchit/it))
    end
    return nothing
end

function get_LECs_params(op)
    targetLECs = String[]
    params = Float64[]; params_ref=Float64[]; pdomains = Tuple{Float64, Float64}[]
    if op=="2n3nall"
        targetLECs= ["ct1_NNLO","ct3_NNLO","ct4_NNLO","cD","cE"]
        params = zeros(Float64,length(targetLECs))
        params_ref = zeros(Float64,length(targetLECs))
        params_ref[1] = -1.10; params_ref[2] = -5.54; params_ref[3] = 4.17
        pdomains = [ (-1.2,-0.5), (-5.0,-2.0), (2.0,6.0), (-2.0,2.0), (-1.5,1.5) ]
    elseif op=="c34"
        targetLECs= ["ct3_NNLO","ct4_NNLO"]
        params = zeros(Float64,length(targetLECs))
        params_ref = zeros(Float64,length(targetLECs))
        params_ref[1] = -3.2; params_ref[2] = 5.4    
        pdomains = [ (-4.5,-2.0), (3.0,6.0)]
    elseif op == "cDE"
        targetLECs= ["cD","cE"]
        params = zeros(Float64,length(targetLECs))
        params_ref = zeros(Float64,length(targetLECs))
        pdomains = [ (-3.0,3.0), (-1.0,1.0)]
    else
        println("warn: op=$op in get_LECs_paramas is not supported!")
        exit()
    end
    return targetLECs, params,params_ref,pdomains
end

"""
cand:: candidate point given by LatinHypercubeSampling
observed:: list of index of `cand` which has been observed
unobserved:: rest candidate indices
history:: array of [logprior,logllh,logpost] for i=1:num_cand

# for GP
Ktt,Kttinv,Ktp,L:: matrix needed for GP calculation
yt:: mean value of training point, to be mean of initial random point 
yscale:: mean&std of yt that is used to scale/rescale data    
acquis:: vector of acquisition function values
pKernel:: hypara for GP kernel, first one is `tau` and the other ones are correlation lengths
adhoc=> tau =1.0, l=1/domain size
"""
function prepOPT(strLECs::LECs,opt,to,io;num_cand=500,op="2n3nall",optimizer="MCMC",MPIcomm=false)
    LECs = strLECs.vals
    idxLECs = strLECs.idxs
    dLECs = strLECs.dLECs
    if opt == false;return nothing;end
    targetLECs, params,params_ref,pdomains = get_LECs_params(op)
    for (k,target) in enumerate(targetLECs)
        idx = idxLECs[target]
        dLECs[target]= params[k] = LECs[idx]
    end 
    pDim = length(targetLECs)

    if optimizer =="LHS" || optimizer=="BayesOpt"
        @assert num_cand > 2 "itnum must be >2"
        Data = [zeros(Float64,pDim) for i=1:num_cand] 
        gens = 200
        @timeit to "LHCoptim" plan, _ = LHCoptim(num_cand+1,pDim,gens)
        tmp = scaleLHC(plan,pdomains)    
        cand = [ tmp[i,:] for i =1:num_cand+1]
        history = [zeros(Float64,3) for i=1:num_cand]
        observed = Int64[ ]
        unobserved = collect(1:num_cand+1)
        if optimizer=="LHS"
            OPTobj = LHSobject(num_cand,targetLECs,params,params_ref,pdomains,cand,observed,unobserved,history)
            Random.seed!(1234)
            propose_LHS!(1,OPTobj,false)
            for (k,target) in enumerate(targetLECs)
                idx = idxLECs[target]
                LECs[idx] = dLECs[target] = params[k]
            end
            return OPTobj
        else
            Ktt = zeros(Float64,num_cand,num_cand) 
            Ktinv = zeros(Float64,num_cand,num_cand)
            tMat = zeros(Float64,num_cand,num_cand)
            Ktp = zeros(Float64,num_cand,1)
            L = zeros(Float64,num_cand,num_cand)
            yt = zeros(Float64,num_cand)
            yscale = zeros(Float64,2) 
            acquis = zeros(Float64,num_cand); acquis[1] = -1.e+10
            pKernel = ones(Float64,pDim+1)
            for i =1:pDim
                tmp = pdomains[i]
                pKernel[i+1] = abs(tmp[2]-tmp[1])^2
            end            
            OPTobj = BOobject(num_cand,targetLECs,params,params_ref,pdomains,pKernel,Data,cand,observed,unobserved,history,Ktt,Ktinv,Ktp,L,tMat,yt,yscale,acquis)
            Random.seed!(1234)
            propose_LHS!(1,OPTobj,false)
            OPTobj.Data[1] .= OPTobj.params
            for (k,target) in enumerate(targetLECs)
                param = params[k]
                idx = idxLECs[target]        
                LECs[idx] = dLECs[target] = param        
            end
            return OPTobj
        end
    elseif optimizer=="MCMC"
        dim = length(params)
        nstep = num_cand
        thining = 1
        acchit = 0
        burnin = div(nstep,10)
        cand = zeros(Float64,dim)
        chain = zeros(Float64,dim,nstep+1)
        chain[:,1] .= params
        history = [zeros(Float64,3) for i=1:nstep]
        if MPIcomm #AffineInv MCMC
            rank_master = 0
            comm = MPI.COMM_WORLD; npsize = MPI.Comm_size(comm)
            chain = zeros(Float64,nstep,npsize,dim)
            myrank = MPI.Comm_rank(comm)
            for n=1:dim; params[n] += 0.2 * randn(); end
            X0 = @view chain[1,myrank+1,:]; X0 .= params
            if myrank == rank_master
                for src_idx = 1:npsize
                    src_rank = src_idx-1
                    if src_rank == rank_master;continue;end
                    X0_ = @view chain[1,src_idx,:]
                    MPI.Recv!(X0_,src_rank,999,comm)
                end
            else
                MPI.Isend(X0,rank_master,999,comm)
            end
            MPI.Barrier(comm)
            ens1 = collect(2:2:npsize)
            ens2 = collect(3:2:npsize)
            a = 2.0 #GW10 recommendation
            OPTobj = MPIMCMCobject(a,dim,nstep,burnin,thining,acchit,targetLECs,params,params_ref,cand,chain,history,ens1,ens2)
            for (k,target) in enumerate(targetLECs)
                idx = idxLECs[target]
                dLECs[target] = LECs[idx] = params[k]
            end 
            return OPTobj
        else
            sigmas = [0.1, 0.3, 0.3, 0.4, 0.2]
            OPTobj = MCMCobject(dim,nstep,burnin,thining,acchit,targetLECs,params,params_ref,sigmas,cand,chain,history)
            return OPTobj
        end
        #println("params_ref ",params_ref)
    else
        @error "optimizer $optimizer not supported"
    end
end

function propose_LHS!(it,OPTobj,BOproposal)
    params = OPTobj.params
    cand = OPTobj.cand
    obs = OPTobj.observed
    unobs = OPTobj.unobserved
    idx = 0
    if BOproposal == false
        tidx = sample(1:length(unobs))
    else
        tidx = find_max_acquisition(it,OPTobj)
    end 
    idx = unobs[tidx]
    deleteat!(unobs,tidx)
    push!(obs,idx)
    params .= cand[idx]
    return nothing
end 

function propose_MH!(it,OPTobj)
    params = OPTobj.params
    ollh = OPTobj.history[it-1][3]
    nllh = OPTobj.history[it][3]
    lograte = nllh - ollh
    logr = log(rand())
    if logr <= lograte
        OPTobj.chain[:,it] .= params
        OPTobj.acchit += 1
    else
        params .= OPTobj.chain[:,it-1]
        OPTobj.chain[:,it] .= OPTobj.chain[:,it-1]
        OPTobj.history[it] .= OPTobj.history[it-1]
    end
    for n = 1:length(params)
        params[n] += OPTobj.sigmas[n] * randn()
    end
    return nothing
end 
function MCMC_HFMBPT(it,OPTobj,HFdata,to;varE=1.0,Lam=0.1)
    eval_HFMBPT(it,OPTobj,HFdata,varE,Lam;mcmc=true)
    if it >1; propose_MH!(it,OPTobj);end
    return nothing
end

function LHS_HFMBPT(it,LHSobj,HFdata,to;varE=1.0,varR=0.25,Lam=0.1)
    eval_HFMBPT(it,LHSobj,HFdata,varE,Lam)
    propose_LHS!(it,LHSobj,false)
    LHSobj.params .= LHSobj.cand[it]  
    return nothing
end

function BO_HFMBPT(it,BOobj,HFdata,to;var_proposal=0.2,varE=1.0,varR=0.25,Lam=0.1)
    params = BOobj.params
    D = length(params); n_ini_BO = 2*D
    ## Update history[it]
    eval_HFMBPT(it,BOobj,HFdata,varE,Lam)
    if it==n_ini_BO
        @timeit to "Kernel" calcKernel!(it,BOobj;ini=true)        
        tmp = [ BOobj.history[i][3] for i=1:it ]
        ymean = mean(tmp); ystd = std(tmp)
        BOobj.yscale[1] = ymean
        BOobj.yscale[2] = ystd
        yt = @view BOobj.yt[1:it]
        yt .= tmp
        #yt .= (tmp .- ymean) ./ ystd ## normalize
    elseif it > n_ini_BO
        @timeit to "Kernel" calcKernel!(it,BOobj)        
        @timeit to "eval p" evalcand(it,BOobj,to)
    end 
    ## Make proposal
    BOproposal = ifelse(it<=n_ini_BO,false,true)
    propose_LHS!(it,BOobj,BOproposal)
    BOobj.Data[it] .= BOobj.params
    return nothing
end
function calcKernel!(it,BOobj;ini=false,eps=1.e-8)
    Ktt = @view BOobj.Ktt[1:it,1:it]
    obs = BOobj.observed
    cand = BOobj.cand
    pKernel = BOobj.pKernel
    tau = pKernel[1]
    Theta = @view pKernel[2:end]
    pdim = length(Theta)
    tv = @view BOobj.tMat[pdim+1:2*pdim,1:1]
    tv2 = @view BOobj.tMat[2*pdim+1:3*pdim,1:1]
    rTr = @view BOobj.tMat[3*pdim+1:3*pdim+1,1:1]
    if ini 
        for i = 1:it
            c_i = cand[obs[i]]
            for j=i:it
                c_j = cand[obs[j]]
                tv .= c_i
                BLAS.axpy!(-1.0,c_j,tv)
                tv2 .= tv .* Theta
                BLAS.gemm!('T','N',1.0,tv,tv2,0.0,rTr)
                Ktt[i,j] = Ktt[j,i] = exp(-0.5*tau*rTr[1])
            end
            Ktt[i,i] += eps
        end
    else
        i = it
        c_i = cand[obs[i]]
        for j=1:it            
            c_j = cand[obs[j]]
            tv .= c_i
            BLAS.axpy!(-1.0,c_j,tv)
            tv2 .= tv .* Theta
            BLAS.gemm!('T','N',1.0,tv,tv2,0.0,rTr)
            Ktt[i,j] = Ktt[j,i] = exp(-0.5*tau*rTr[1])
        end
        Ktt[i,i] += eps
    end
    ## Calculate Ktt^{-1} 
    Ktinv = @view BOobj.Ktinv[1:it,1:it]
    L = @view BOobj.L[1:it,1:it]
    try
        myCholesky!(Ktt,it,L)
    catch
        println("Theta $Theta")
        for i=1:size(Ktt)[1]
            print_vec("",@view Ktt[i,:])
        end
        exit()
    end
    Linv = inv(L)
    BLAS.gemm!('T','N', 1.0,Linv,Linv,0.0,Ktinv)
    return nothing
end

function calcKtp!(it,xp,BOobj)
    tau = BOobj.pKernel[1]
    Theta = @view BOobj.pKernel[2:end]
    pdim = length(Theta)
    Ktp = @view BOobj.Ktp[1:it,1:1]
    tv = @view BOobj.tMat[pdim+1:2*pdim,1:1]
    tv2 = @view BOobj.tMat[2*pdim+1:3*pdim,1:1]
    rTr = @view BOobj.tMat[3*pdim+1:3*pdim+1,1:1]
    obs = BOobj.observed
    cand = BOobj.cand
    for i = 1:it
        tv .= cand[obs[i]]
        BLAS.axpy!(-1.0,xp,tv)
        tv2 .= tv .* Theta
        BLAS.gemm!('T','N',1.0,tv,tv2,0.0,rTr)
        Ktp[i] = exp(-0.5*tau*rTr[1])
    end
end 

function evalcand(it,BOobj,to;epsilon=1.e-9)
    Ktt = @view BOobj.Ktt[1:it,1:it]
    Ktp = @view BOobj.Ktp[1:it,1:1]
    Ktinv = @view BOobj.Ktinv[1:it,1:it]
    tM1 = @view BOobj.tMat[1:1,1:it]
    tM2 = @view BOobj.tMat[2:2,1:1]    
    yt = @view BOobj.yt[1:it]
    unobs = BOobj.unobserved
    cand = BOobj.cand
    fplus = BOobj.acquis[1]
    tau = BOobj.pKernel[1]
    fAs = @view BOobj.acquis[2:1+length(unobs)]
    fAs .= 0.0
    for (n,idx) in enumerate(unobs)
        xp = cand[idx]
        calcKtp!(it,xp,BOobj)
        BLAS.gemm!('T','N',1.0,Ktp,Ktinv,0.0,tM1)       
        BLAS.gemm!('N','N',1.0,tM1,Ktp,0.0,tM2)
        mup = dot(tM1,yt)
        sigma = tau
        try
            sigma = sqrt(tau + epsilon - tM2[1])
        catch   
            Imat = Matrix{Float64}(I,it,it)
            Mat = Ktt*Ktinv
            tnorm = norm(Imat-Mat,Inf)
            Theta = @view BOobj.pKernel[2:end]
            println("tau-tM2[1] ",tau-tM2[1], "  tnorm $tnorm")
            print_vec("xp",xp)
            print_vec("",Ktp)
            if it < 10
                println("Ktt ",isposdef(Ktt))
                for i=1:size(Ktt)[1]
                    print_vec("",@view Ktt[i,:])
                end
                println("Ktinv")
                for i=1:size(Ktinv)[1]
                    print_vec("",@view Ktinv[i,:])
                end
                println("Mat")
                for i=1:size(Mat)[1]
                    print_vec("",@view Mat[i,:])
                end
                println("Ktp \n$Ktp")
            end
            exit()
        end
        Z = (mup-fplus) / sigma
        fAs[n] = Z*sigma * fPhi(Z) + exp(-0.5*Z^2)/sqrt(2.0*pi)
    end 
    return nothing
end

function find_max_acquisition(it,BOobj)
    unobs = BOobj.unobserved
    fAs = @view BOobj.acquis[2:1+length(unobs)]
    idx = argmax(fAs)
    return idx
end 

function fPhi(Z)
    return  0.5 * erfc(-(Z/sqrt(2.0)))
end


"""
function used for proposals in Affine invariant MCMC
"""
function gz(a=2.0) 
    return (((a-1)*rand() +1)^2) / a
end

function mpi_hfmbpt(itnum,chiEFTobj,OPTobj,d9j,HOBs,nucs,HFdata,to,io;
                    Operators=["Rp2"],rank_master=0,writesnt=false,debug=false)   
    LECs = chiEFTobj.LECs.vals; idxLECs = chiEFTobj.LECs.idxs; dLECs = chiEFTobj.LECs.dLECs
    myrank = MPI.Comm_rank(MPI.COMM_WORLD)       
    comm = MPI.COMM_WORLD;myrank = MPI.Comm_rank(comm);npsize = MPI.Comm_size(comm)
    chain = OPTobj.chain
    for it = 1:itnum
        if it == 1
            if myrank == rank_master
                for dst_idx = 1:npsize                
                    dst_rank = dst_idx-1
                    if dst_rank == rank_master;continue;end
                    X0 = @view chain[1,dst_idx,:]
                    MPI.Isend(X0,dst_rank,999,comm)
                end
            else
                X0 = @view chain[1,myrank+1,:]
                MPI.Recv!(X0,rank_master,999,comm)
                OPTobj.params .= X0
            end
            MPI.Barrier(comm)
            ## update LECs and perform HF-MBPT calculation
            updateLECs_in_chiEFTobj!(chiEFTobj,OPTobj.targetLECs,OPTobj.params)
            calc_vmom_3nf(chiEFTobj,it,to)
            add_V12mom!(chiEFTobj.V12mom,chiEFTobj.V12mom_2n3n)           
            dicts_tbme = TMtrans(chiEFTobj,HOBs,to;writesnt=false)
            if !debug
                hf_main_mem(chiEFTobj,nucs,dicts_tbme,d9j,HOBs,HFdata,to;Operators=Operators,io=io)
            end
            eval_HFMBPT(it,OPTobj,HFdata,0.1,1.0;io=io,debug=debug)
            print_vec("it = "*@sprintf("%8i",it),OPTobj.params,io)
        else
            walker_i = myrank + 1
            candidate = OPTobj.cand
            S_at_t = @view chain[it-1,:,:]
            if myrank == rank_master #master=>worker
                for dst_idx = 1:npsize
                    dst_rank = dst_idx -1
                    if dst_rank == rank_master;continue;end
                    MPI.Isend(S_at_t,dst_rank,99,comm)
                end
            else #worker <= master
                MPI.Recv!(S_at_t,rank_master,99,comm)
            end    
            MPI.Barrier(comm)
            for nbatch = 0:1
                subset = ifelse(nbatch==0,OPTobj.ens1,OPTobj.ens2)
                S_complement = ifelse(nbatch==0,OPTobj.ens2,OPTobj.ens1)        
                if myrank != rank_master 
                    if walker_i % 2 == nbatch
                        Xi = @view S_at_t[walker_i,:]
                        walker_j = sample(S_complement)
                        Xj = @view S_at_t[walker_j,:]               
                        zval = gz(OPTobj.a)
                        candidate .= Xj + zval .*  (Xi - Xj)
                        #print_vec("myrank $myrank walker i/j $walker_i $walker_j Xi $Xi Xj $Xj cand ",candidate)
                        OPTobj.params .= candidate
                        updateLECs_in_chiEFTobj!(chiEFTobj,OPTobj.targetLECs,candidate)
                        calc_vmom_3nf(chiEFTobj,it,to)
                        add_V12mom!(chiEFTobj.V12mom,chiEFTobj.V12mom_2n3n)                      
                        dicts_tbme = TMtrans(chiEFTobj,HOBs,to;writesnt=false)
                        if !debug
                            hf_main_mem(chiEFTobj,nucs,dicts_tbme,d9j,HOBs,HFdata,to;Operators=Operators,io=io)
                        end
                        eval_HFMBPT(it,OPTobj,HFdata,0.1,1.0;io=io,debug=debug)
                        logratio = 1.0
                        oeval = OPTobj.history[it-1]; neval = OPTobj.history[it]
                        logratio = (OPTobj.dim -1) * log(zval) + neval[3] - oeval[3]
                        Accept = ifelse(log(rand())<logratio,true,false)
                        if Accept 
                            Xi .= candidate
                            OPTobj.acchit += 1
                        else
                            neval .= oeval
                            for (k,target) in enumerate(OPTobj.targetLECs)
                                idx = idxLECs[target]
                                LECs[idx] = dLECs[target] = Xi[k] 
                            end   
                        end
                        OPTobj.params .= Xi
                        print_vec("it = "*@sprintf("%8i",it),Xi,io)
                        MPI.Isend(Xi,0,myrank,comm) # worker => master
                    end
                else # master <= worker
                    for walkernum_src in subset
                        src_rank = walkernum_src -1 
                        tmp = @view chain[it,walkernum_src,:]
                        MPI.Recv!(tmp,src_rank,src_rank,comm)
                    end
                end
                MPI.Barrier(comm) 
            end
            MPI.Barrier(comm) 
        end
    end
    MPI.Finalize()
    return nothing
end

"""

Reference:
Goodman & Weare, "Ensemble samplers with affine invariance", Communications in Applied Mathematics and Computational Science, DOI: 10.2140/camcos.2010.5.65, 2010.
"""
function sample_AffineInvMCMC(numwalkers::Int, x0::Matrix{Float64}, 
                              nstep::Integer, thinning::Integer,a::Float64=2.)
	@assert length(size(x0)) == 2
	ndim = size(x0)[1]
    chain = zeros(Float64,ndim,numwalkers,nstep)    
    lLHS = zeros(Float64,1,numwalkers,nstep)
    for n =1:ndim
        tX = @view chain[n,:,1]
        tX .= x0[n,:]
    end
    for walker = 1:numwalkers
        Xi = @view chain[:,walker,1]
        lLHS[1,walker,1] = myllh(Xi)
    end    
    idx_S0 = collect(1:div(numwalkers,2))
    idx_S1 = collect(div(numwalkers,2)+1:numwalkers)
    batch = [(idx_S0,idx_S1),(idx_S1,idx_S0)]
    acchit = 0
    for t = 2:nstep
        for nbatch = 1:2
            idxs_target,idxs_complement = batch[nbatch]
            for walker in idxs_target                
                Xi = @view chain[:,walker,t-1]
                nX = @view chain[:,walker,t]
                llh = lLHS[1,walker,t-1]
                walker2 = sample(idxs_complement)
                Xj = @view chain[:,walker2,t-ifelse(nbatch==1,1,0)]
                z = (((a-1)*rand() + 1)^2) / a
                nX .= Xj + z .* (Xi -Xj)
                nllh = myllh(nX)
                logratio = (ndim-1)*log(z) + nllh - llh
                if log(rand(rng)) < logratio
                    lLHS[1,walker,t] = nllh
                    acchit += 1
                else
                    lLHS[1,walker,t] = llh
                    nX .= Xi    
                end
            end
        end
    end
    println("Acc. rate: ",@sprintf("%6.2f",100*acchit/((nstep-1)*numwalkers)))
    chain,lLHS = flatten_mcmcarray(chain,lLHS)
    return chain, lLHS
end
function flatten_mcmcarray(chain::Array, llhoodvals::Array,order=true)
	numdims, numwalkers, numsteps = size(chain)
	newchain = Array{Float64}(undef, numdims, numwalkers * numsteps)
    for j = 1:numsteps
        for i = 1:numwalkers
            if order
                newchain[:, i + (j - 1) * numwalkers] .= chain[:, i, j]
            else
                newchain[:, i + (j - 1) * numwalkers] .= chain[j,:,i]
            end
        end
    end              
	return newchain, llhoodvals[1:end]
end