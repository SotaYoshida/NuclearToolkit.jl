"""
    make_chiEFTint(;is_plot=false,BO=false,optHFMBPT=false,itnum=20,writesnt=true,nucs=[])

main function in chiEFTint.
This generates NN-potential in momentum space and then transforms it in HO basis.
This function is exported and can be simply called make_chiEFTint() in the run script. 

# Optional arguments: Note that these are mainly for the author's purposes, so you do not specify these.
- `optHFMBPT::Bool` to optimize LECs using BO&HFMBPT
- `itnum::Int=20` number of iteration for BO&HFMBPT
- `is_plot::Bool`, to visualize optimization process of LECs
- `writesnt::Bool`, to write out interaction file in snt (KSHELL) format. ```julia writesnt = false``` case can be usefull when you repeat large number of calculations for different LECs.
"""
function make_chiEFTint(;optHFMBPT=false,itnum=1,is_show=false,writesnt=true,nucs=[],is_plot=false,optimizer="",MPIcomm=false,io=stdout)
    if optimizer!="";optHFMBPT=true;end
    if nucs == []; optHFMBPT=false;end
    if MPIcomm; 
        MPI.Init()
        if !isdir("mpilog");run(`mkdir mpilog`);end
        myrank = MPI.Comm_rank(MPI.COMM_WORLD)
        io = open("./mpilog/log_rank"*string(myrank)*".dat","w")
    else
        io = open("logfile.dat","w")
    end
    chiEFTobj = init_chiEFTparams(;io=io)
    emax = chiEFTobj.emax
    pmax_fm = chiEFTobj.pmax_fm
    n_mesh = chiEFTobj.n_mesh
    Pmax_fm = chiEFTobj.Pmax_fm
    hw = chiEFTobj.hw
    v_chi_order = chiEFTobj.v_chi_order

    to = TimerOutput()
    @timeit to "prep." begin
        dict6j,d6j_nabla,d6j_int = PreCalc6j(emax)
        # prep. momentum mesh
        xr_fm,wr = Gauss_Legendre(0.0,pmax_fm,n_mesh)
        xr = xr_fm * hc
        numst2,dict_numst,arr_numst = bstate(;io=io)
        V12mom = [ zeros(Float64,n_mesh,n_mesh) for i=1:length(numst2)]
        V12mom_2n3n = [ zeros(Float64,n_mesh,n_mesh) for i=1:length(numst2)]
        ## prep. radial functions
        rmass = Mp*Mn/(Mp+Mn)
        br = sqrt(hc^2 /(rmass*hw))
        Rnl = Rnl_all_ab(chiEFTobj,lmax,br,n_mesh,xr_fm)
        ## prep. for valence space oparators
        ntmp = 3; if v_chi_order > 0;ntmp=n_mesh;end
        xrP_fm,wrP = Gauss_Legendre(0.0,Pmax_fm,ntmp)
        xrP = xrP_fm .* hc   
        RNL = Rnl_all_ab(chiEFTobj,lcmax,br,ntmp,xrP_fm)
        ## prep. for partial-wave decompositon
        lsjs = [ [ [J,J,0,J],[J,J,1,J],[J+1,J+1,1,J],[J-1,J-1,1,J],[J+1,J-1,1,J],[J-1,J+1,1,J]] for J = 0:jmax]
        llpSJ_s = [[0,0,1,1],[1,1,1,0],[1,1,0,1],[1,1,1,1],[0,0,0,0],[0,2,1,1],[1,1,1,2]]
        tllsj = zeros(Int,5)
        opfs = [ zeros(Float64,11) for i=1:5]#T,SS,C,LS,SL
        f_ss!(opfs[2]);f_c!(opfs[3])
        ## prep. Gauss point for integrals
        ts, ws = Gauss_Legendre(-1,1,96)
        ## prep. for TBMEs
        jab_max = 4*emax + 2
        Numpn= Dict( [0,0,0,0] => 0 ) ;delete!(Numpn,[0,0,0,0])   
        infos,izs_ab,nTBME = make_sp_state(chiEFTobj,jab_max,Numpn;io=io)
        println(io,"# of channels 2bstate ",length(infos)," #TBME = $nTBME")
        ## prep. integrals for 2n3n
        F0s = zeros(Float64,n_mesh); F1s = zeros(Float64,n_mesh); F2s = zeros(Float64,n_mesh); F3s = zeros(Float64,n_mesh)
        wsyms=[]
        QWs = prep_QWs(chiEFTobj,xr,ts,ws,to)
        if chiEFTobj.calc_3N
            prep_Fis!(chiEFTobj,xr,F0s,F1s,F2s,F3s)
            wsyms = prep_wsyms()
        end
    end
    rdict6j = nothing
    if optHFMBPT
        rdict6j = adhoc_rewrite6jdict(emax,dict6j)
    end
    HFdata = prepHFdata(nucs,"",["E"],"")
    ### LECs
    LECs = Float64[ ]
    idxLECs=Dict{String,Int64}()
    dLECs=Dict{String,Float64}()
    fn_LECs = get_fn_LECs(chiEFTobj.pottype)
    read_LECs!(LECs,idxLECs,dLECs;initialize=true,inpf=fn_LECs)
    
    ## Start Opt stuff    
    OPTobj = prepOPT(LECs,idxLECs,dLECs,optHFMBPT,to,io;num_cand=itnum,optimizer=optimizer,MPIcomm=MPIcomm) 
    d9j = HOBs = nothing
    if optHFMBPT          
        d9j,HOBs = PreCalcHOB(emax,chiEFTobj,d6j_int,to;io=io)
    end
    ## END: BO stuff

    ### Calculation of NN potential and SRG
    X9,U6 = prepareX9U6(2*emax)
    @timeit to "NNcalc" begin 
        if chiEFTobj.calc_NN
            OPEP(chiEFTobj,ts,ws,xr,V12mom,dict_numst,to,lsjs,llpSJ_s,tllsj,opfs,QWs[end])
            LO(chiEFTobj,xr,dLECs,V12mom,dict_numst,to)                
            if chiEFTobj.chi_order >= 1
                NLO(chiEFTobj,xr,dLECs,V12mom,dict_numst,to)
                ### TPE terms (NLO,NNLO,N3LO,N4LO)
                tpe(chiEFTobj,dLECs,ts,ws,xr,V12mom,dict_numst,to,llpSJ_s,lsjs,tllsj,opfs)
            end                
            if chiEFTobj.chi_order >= 3
                N3LO(chiEFTobj,xr,dLECs,V12mom,dict_numst,to)
            end
            if chiEFTobj.chi_order >= 4
                N4LO(chiEFTobj,xr,dLECs,V12mom,dict_numst,to)
            end
        end
    end
    if chiEFTobj.srg
        @timeit to "SRG" SRG(chiEFTobj,xr_fm,wr,V12mom,dict_numst,to)
    end
    if !MPIcomm
        for it = 1:itnum
            if it > 1; for i=1:length(numst2); V12mom_2n3n[i] .= 0.0; end;end
            if chiEFTobj.calc_3N
                calc_vmom_3nf(chiEFTobj,dLECs,ts,ws,xr,V12mom_2n3n,dict_numst,to,F0s,F1s,F2s,F3s,QWs,wsyms,lsjs,llpSJ_s,tllsj)
            end
            add_V12mom!(V12mom,V12mom_2n3n)            
            #transform mom. int. to HO matrix element
            @timeit to "Vtrans" begin
                V12ab = Vrel(chiEFTobj,V12mom_2n3n,numst2,xr_fm,wr,n_mesh,Rnl,to)
                dicts_tbme = TMtrans(chiEFTobj,dLECs,xr,wr,xrP,wrP,Rnl,RNL,nTBME,infos,izs_ab,Numpn,V12ab,arr_numst,dict6j,d6j_nabla,X9,U6,to;writesnt=writesnt)
            end
            ## If you want to optimize (or try samplings) change itnum, insert a function to update/optimize/sample the LECs here     
            if nucs != [ ] && writesnt == false && optHFMBPT 
                print_vec("it = "*@sprintf("%8i",it),OPTobj.params,io)
                @timeit to "HF/HFMBPT" hf_main_mem(chiEFTobj,nucs,dicts_tbme,rdict6j,HFdata,d9j,HOBs,to;Operators=["Rp2"])
                if optimizer=="BayesOpt"
                    BO_HFMBPT(it,OPTobj,HFdata,to)
                elseif optimizer=="LHS"
                    LHS_HFMBPT(it,OPTobj,HFdata,to)
                elseif optimizer=="MCMC"
                    MCMC_HFMBPT(it,OPTobj,HFdata,to)
                end
                for (k,target) in enumerate(OPTobj.targetLECs)
                    idx = idxLECs[target]
                    LECs[idx] = dLECs[target] = OPTobj.params[k] 
                end
            end
            if !optHFMBPT;break;end
        end
    else
        myrank = MPI.Comm_rank(MPI.COMM_WORLD)
        for it = 1:itnum
            mpi_hfmbpt(it,OPTobj,chiEFTobj,LECs,idxLECs,dLECs,nucs,rdict6j,HFdata,d9j,HOBs,to,io,
                       ts,ws,xr,wr,V12mom,V12mom_2n3n,dict_numst,F0s,F1s,F2s,F3s,QWs,wsyms,lsjs,llpSJ_s,tllsj,
                       numst2,xr_fm,n_mesh,xrP,wrP,Rnl,RNL,nTBME,infos,izs_ab,Numpn,arr_numst,dict6j,d6j_nabla,X9,U6;Operators=[])
        end
        MPI.Finalize()
    end
    if io != stdout; close(io);end
    if is_show; show(to, allocations = true,compact = false);println("");end
    return true
end

function gz(a=2.0) 
    return (((a-1)*rand() +1)^2) / a
end

function mpi_hfmbpt(t,OPTobj,chiEFTobj,LECs,idxLECs,dLECs,nucs,rdict6j,HFdata,d9j,HOBs,to,io,
                    ts,ws,xr,wr,V12mom,V12mom_2n3n,dict_numst,F0s,F1s,F2s,F3s,QWs,wsyms,lsjs,llpSJ_s,tllsj,
                    numst2,xr_fm,n_mesh,xrP,wrP,Rnl,RNL,nTBME,infos,izs_ab,Numpn,arr_numst,dict6j,d6j_nabla,X9,U6;
                    Operators=["Rp2"],rank_master=0,writesnt=false,debug=false)
    comm = MPI.COMM_WORLD;myrank = MPI.Comm_rank(comm);npsize = MPI.Comm_size(comm)
    chain = OPTobj.chain
    if t == 1
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
        for (k,target) in enumerate(OPTobj.targetLECs)
            idx = idxLECs[target]
            LECs[idx] = dLECs[target] = OPTobj.params[k]
        end

        for i=1:length(numst2); V12mom_2n3n[i] .= 0.0; end
        calc_vmom_3nf(chiEFTobj,dLECs,ts,ws,xr,V12mom_2n3n,dict_numst,to,F0s,F1s,F2s,F3s,QWs,wsyms,lsjs,llpSJ_s,tllsj)
        add_V12mom!(V12mom,V12mom_2n3n)
        V12ab = Vrel(chiEFTobj,V12mom_2n3n,numst2,xr_fm,wr,n_mesh,Rnl,to)
        dicts_tbme = TMtrans(chiEFTobj,dLECs,xr,wr,xrP,wrP,Rnl,RNL,nTBME,infos,izs_ab,Numpn,V12ab,arr_numst,dict6j,d6j_nabla,X9,U6,to;writesnt=writesnt)        
        if !debug
            hf_main_mem(chiEFTobj,nucs,dicts_tbme,rdict6j,HFdata,d9j,HOBs,to;Operators=Operators,io=io)
        end
        eval_HFMBPT(t,OPTobj,HFdata,0.1,1.0;io=io,debug=debug)
        print_vec("it = "*@sprintf("%8i",1),OPTobj.params,io)
        return nothing
    end
    walker_i = myrank + 1
    candidate = OPTobj.cand
    S_at_t = @view chain[t-1,:,:]
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
                for (k,target) in enumerate(OPTobj.targetLECs)
                    idx = idxLECs[target]
                    LECs[idx] = dLECs[target] = candidate[k] 
                end    
                for i=1:length(numst2); V12mom_2n3n[i] .= 0.0; end
                calc_vmom_3nf(chiEFTobj,dLECs,ts,ws,xr,V12mom_2n3n,dict_numst,to,F0s,F1s,F2s,F3s,QWs,wsyms,lsjs,llpSJ_s,tllsj)
                add_V12mom!(V12mom,V12mom_2n3n)
                V12ab = Vrel(chiEFTobj,V12mom_2n3n,numst2,xr_fm,wr,n_mesh,Rnl,to)
                dicts_tbme = TMtrans(chiEFTobj,dLECs,xr,wr,xrP,wrP,Rnl,RNL,nTBME,infos,izs_ab,Numpn,V12ab,arr_numst,dict6j,d6j_nabla,X9,U6,to;writesnt=writesnt)        
                if !debug
                    hf_main_mem(chiEFTobj,nucs,dicts_tbme,rdict6j,HFdata,d9j,HOBs,to;Operators=Operators,io=io)
                end
                eval_HFMBPT(t,OPTobj,HFdata,0.1,1.0;io=io,debug=debug)
                logratio = 1.0
                if t > 1
                    oeval = OPTobj.history[t-1]; neval = OPTobj.history[t]
                    logratio = (OPTobj.dim -1) * log(zval) + neval[3] - oeval[3]
                end
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
                print_vec("it = "*@sprintf("%8i",t),Xi,io)
                MPI.Isend(Xi,0,myrank,comm) # worker => master
            end
        else # master <= worker
            for walkernum_src in subset
                src_rank = walkernum_src -1 
                tmp = @view chain[t,walkernum_src,:]
                MPI.Recv!(tmp,src_rank,src_rank,comm)
            end
        end
        MPI.Barrier(comm) 
    end
    MPI.Barrier(comm) 
    return nothing
end

function add_V12mom!(V12mom,V12mom_2n3n,a=1.0)
    for i = 1:length(V12mom)
        V12mom_2n3n[i] .+= V12mom[i] * a
    end
    return nothing
end

function get_fn_LECs(pottype)
    fn = ""
    if pottype =="em500n3lo"
        fn = "src/chiEFTint/LECs.jl"
    elseif pottype == "emn500n3lo"
        fn = "src/chiEFTint/LECs_EMN500N3LO.jl"
    elseif pottype == "emn500n4lo"
        fn = "src/chiEFTint/LECs_EMN500N4LO.jl"
    else
        println("warn!! potype=$pottype is not supported now!")
        exit()
    end
    return fn
end

"""
    SRG(xr,wr,V12mom,dict_numst,to)

main function for Similarity Renormalization Group (SRG) transformation of NN interaction in CM-rel momentum space.
"""
function SRG(chiEFTobj,xr,wr,V12mom,dict_numst,to)
    n_mesh = chiEFTobj.n_mesh
    srg_lambda = chiEFTobj.srg_lambda
    l1s = [0,1,2,3,4,5,6,1,1,2,3,4,5,6,0,1,2,3,4,5]
    l2s = [0,1,2,3,4,5,6,1,1,2,3,4,5,6,2,3,4,5,6,7]
    Ss  = [0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1]
    Js  = [0,1,2,3,4,5,6,0,1,2,3,4,5,6,1,2,3,4,5,6]

    ndim = n_mesh*2
    n_ode = div(ndim*(ndim+1),2)
    nthre = nthreads()
    Vs = [ zeros(Float64,ndim,ndim) for i = 1:nthre]
    Ts = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    Hs = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    Hts = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    etas = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    Rs = [zeros(Float64,ndim,ndim) for i = 1:nthre]
    tkeys = [ zeros(Int64,5) for i=1:nthre]
    
    sSRG = (1.0/srg_lambda)^4 # s in literature
    ds = 1.e-4
    srange = 0.0:ds:sSRG
    numit = length(srange)
    nch = length(l1s)
    for nth =1:nch*3
        threid = threadid()
        tkey = tkeys[threid]
        V = Vs[threid]; T = Ts[threid]; H = Hs[threid]
        Ht= Hts[threid];eta = etas[threid]; R = Rs[threid]
        iz = -1
        if nth <= nch
            nothing
        elseif nth <= 2*nch
            iz = 0
        else
            iz = 1
        end
        pnrank = iz + 2
        tdict = dict_numst[pnrank]
        rmass = Mp * Mn / (Mp+Mn) 
        if iz == -1; rmass = 0.5 * Mp;end
        if iz ==  1; rmass = 0.5 * Mn;end
        face =  (2.0*rmass)^2 / hc4

        ich = ifelse(nth%nch!=0,nth%nch,nch)
        l1 = l1s[ich];l2 = l2s[ich];S = Ss[ich];J = Js[ich]
        if abs(iz % 2) == 1 && (l1+S+abs(iz)) % 2 == 0; continue;end
        tkey[1] = 2*iz; tkey[2]=l1; tkey[3]=l2; tkey[4]=S;tkey[5]=J
        if l1 == l2
            V12idx = tdict[tkey]; tv = V12mom[V12idx]
            sV  = @view  V[1:n_mesh,1:n_mesh]          
            sT  = @view  T[1:n_mesh,1:n_mesh]
            sH  = @view  H[1:n_mesh,1:n_mesh]
            sHt = @view Ht[1:n_mesh,1:n_mesh]
            seta = @view eta[1:n_mesh,1:n_mesh]
            tR = @view R[1:n_mesh,1:n_mesh]                
            for (i,x) in enumerate(xr)
                wri = wr[i]
                for (j,y) in enumerate(xr)
                    sH[i,j] = sV[i,j] = tv[i,j] * x * y * sqrt(wri*wr[j])
                end
                sT[i,i] = (x *hc)^2 / (2*rmass)
                sH[i,i] += (x *hc)^2 / (2*rmass)
            end
            if norm(sH-sH',Inf) > 1.e-9; println(" norm(sH-sH') ", norm(sH-sH',Inf));end
            srg_tr(sH,sT,sHt,sV,seta,tR,sSRG,face,ds,numit,to)
            ## overwrite V12
            for (i,x) in enumerate(xr)
                for (j,y) in enumerate(xr) 
                    tv[i,j] = (sV[i,j]-sT[i,j]) / ( x * y * sqrt( wr[i] * wr[j]))
                end
            end
        elseif l1+2==l2
            ## l1-l1, l2-l2, l1-l2, l2-l1 part
            tkey[2]=l1; tkey[3]=l1; V12idx = tdict[tkey]; tv1 = V12mom[V12idx]
            tkey[2]=l2; tkey[3]=l2; V12idx = tdict[tkey]; tv2 = V12mom[V12idx]
            tkey[2]=l1; tkey[3]=l2; V12idx = tdict[tkey]; tv3 = V12mom[V12idx]
            tkey[2]=l2; tkey[3]=l1; V12idx = tdict[tkey]; tv4 = V12mom[V12idx]
            for (i,x) in enumerate(xr)
                for (j,y) in enumerate(xr) 
                    tfac = x * y * sqrt(wr[i] * wr[j])
                    V[i,j] = tv1[i,j] * tfac 
                    V[n_mesh+i,n_mesh+j] = tv2[i,j] * tfac 
                    V[i,n_mesh+j] = tv3[i,j] * tfac 
                    V[n_mesh+i,j] = tv4[i,j] * tfac 
                end
                T[i,i] = T[n_mesh+i,n_mesh+i] = (x *hc)^2 / (2*rmass)
            end
            H .= V;H .+= T
            srg_tr(H,T,Ht,V,eta,R,sSRG,face,ds,numit,to)
            H .= V; H .-= T # Veff = H(s) - T # H is reused as Veff            
            for (i,x) in enumerate(xr)
                for (j,y) in enumerate(xr)
                    deno =  x * y * sqrt(wr[i]*wr[j])
                    tv1[i,j] = H[i,j] / deno # ( x * y * sqrt(wr[i]*wr[j]))
                    tv2[i,j] = H[n_mesh+i,n_mesh+j] /  deno #(x*y *sqrt(wr[i]*wr[j]))
                    tv3[i,j] = H[i,n_mesh+j] / deno #( x*y * sqrt(wr[i]*wr[j]))
                    tv4[i,j] = H[n_mesh+i,j] / deno # ( x*y * sqrt(wr[i]*wr[j]))
                end
            end
        end
    end
    return nothing
end

"""
    commutator(A,B,R,fac)

wrapper function to overwrite ```R``` by ```fac*(AB-BA)```
"""
function commutator(A,B,R,fac)
    BLAS.gemm!('N', 'N',  fac, A, B, 0.0, R)
    BLAS.gemm!('N', 'N', -fac, B, A, 1.0, R)
    return nothing
end

"""
    RKstep(T,Ho,eta,R,faceta,fRK,Ht)

wrapper function to calc. a Runge-Kutta (RK) step
"""
function RKstep(T,Ho,eta,R,faceta,fRK,Ht)
    BLAS.gemm!('N', 'N',  faceta, T, Ho, 0.0, eta)
    BLAS.gemm!('N', 'N', -faceta, Ho, T, 1.0, eta) # =>eta
    BLAS.gemm!('N', 'N',  1.0, eta, Ho, 0.0, R)
    BLAS.gemm!('N', 'N', -1.0, Ho, eta, 1.0, R)    
    BLAS.axpy!(fRK,R,Ht)    
    return nothing
end
function RKstep_mul(T,Ho,eta,R,faceta,fRK,Ht)
    mul!(eta,T,Ho,faceta,0.0)
    mul!(eta,Ho,T,-faceta,1.0)
    mul!(R,eta,Ho,1.0,0.0)
    mul!(R,Ho,eta,-1.0,1.0)
    BLAS.axpy!(fRK,R,Ht)    
    return nothing
end

"""
    srg_tr(Ho,T,Ht,Hs,eta,R,sSRG,face,ds,numit,to; r_err=1.e-8,a_err=1.e-8,tol=1.e-6)

to carry out SRG transformation
"""
function srg_tr(Ho,T,Ht,Hs,eta,R,sSRG,face,ds,numit,to;
                r_err=1.e-8,a_err=1.e-8,tol=1.e-6)
    Hs .= Ho
    Ht .= Hs
    for it = 1:numit
        if it ==numit; ds = sSRG -(numit-1)*ds;end
        Ho .= Hs        
        RKstep(T,Ho,eta,R,face,ds/6.0,Ht)

        Ho .= Hs; BLAS.axpy!(0.5*ds,R,Ho)
        RKstep(T,Ho,eta,R,face,ds/3.0,Ht)

        Ho .= Hs; BLAS.axpy!(0.5*ds,R,Ho)
        RKstep(T,Ho,eta,R,face,ds/3.0,Ht)

        Ho .= Hs; BLAS.axpy!(ds,R,Ho)
        RKstep(T,Ho,eta,R,face,ds/6.0,Ht)

        Hs .= Ht
    end
    return nothing
end

function genLaguerre(n,alpha,x)
    if n==0
        return 1
    elseif n==1
        return -x + (alpha+1)
    else
        s = 0.0
        for i=0:n
            bfac = gamma(n+alpha+1) / ( gamma(n-i+1) * gamma(alpha+i+1))
            s += (-1)^i * x^i / factorial(i) * bfac
        end        
        return s
    end
end

"""
    Gauss_Legendre(xmin,xmax,n;eps=3.e-16) 

to calculate mesh points and weights for Gauss-Legendre quadrature
"""
function Gauss_Legendre(xmin,xmax,n;eps=3.e-16) 
    m = div(n+1,2)
    x = zeros(Float64,n); w = zeros(Float64,n)
    xm = 0.5 * (xmin+xmax); xl = 0.5 * (xmax-xmin)
    p1 = p2 = p3 = pp = z = z1 = 0.0
    for i =1:m
        z1 = 0.0
        z = cos((pi*(i-0.25))/(n+0.5))
        hit = 0        
        while true
            hit += 1
            p1 = 1.0; p2 = 0.0
            for j=1:n
                p3=p2
                p2=p1
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j
            end
            pp = n*(z*p1-p2)/(z^2 -1.0)
            z1 = z
            z -= p1/pp
            if abs(z-z1) < eps; break;end
            if hit > 100; println("warn! in Gausse_Legendre");exit();end
        end
        x[i]=xm-xl*z
        x[n+1-i]=xm+xl*z
        w[i]=2.0*xl/( (1.0-z^2) * pp^2)
        w[n+1-i]=w[i]
    end
    return x,w
end

""" 
    QL(z,J::Int64,ts,ws)

To calculate Legendre functions of second kind by Gauss-Legendre quadrature
"""
function QL(z,J::Int64,ts,ws,QLdict)
    s = 0.0
    if J==0
        s = 0.5 * log( abs((1.0+z)/(1.0-z)) )
    elseif J==1
        s = 0.5 *z * log( abs((1.0+z)/(1.0-z)) ) -1.0
    elseif J==2
        s = 0.25 * (3.0*z^2 -1.0) * log( abs((1.0+z)/(1.0-z)) ) - 1.5*z
    elseif J==3
        s = 0.25 * (5.0*z^3 -3.0*z) * log( abs((1.0+z)/(1.0-z)) ) - 2.5*z^2 + 2.0/3.0
    else
        if J >= 0
            t = get(QLdict[J+1],z,0.0) 
            if t == 0.0
                @inbounds for (i,t) in enumerate(ts)
                    s += ws[i] * (1.0-t*t)^J / (z-t)^(J+1) 
                end 
                s *= 1.0 / 2.0^(J+1)
                QLdict[J+1][z] = s
            else
                return t
            end
        end
    end
    return s
end

function make_sp_state(chiEFTobj,jab_max,Numpn;io=stdout)
    emax = chiEFTobj.emax
    kh = Dict( [0,0] => 0 ) ;delete!(kh,[0,0])
    kn = Dict( [0,0] => 0 ) ;delete!(kn,[0,0])
    kl = Dict( [0,0] => 0 ) ;delete!(kl,[0,0])
    kj = Dict( [0,0] => 0 ) ;delete!(kj,[0,0])    
    maxsps = Int((emax+1) * (emax+2) / 2)
    println(io,"# of sp states $maxsps")    
    n = 0
    for NL=0:emax
        for L=0:NL
            if (NL-L) % 2 == 1;continue;end
            Nn= div(NL-L,2)
            for IS=-1:2:1
                jd=2*L+IS
                if jd < 0; continue;end
                n=n+1
                kh[[-1,n]]=3;  kh[[1,n]]=3
                kn[[-1,n]]=Nn; kn[[1,n]]=Nn
                kl[[-1,n]]=L;  kl[[1,n]]=L
                kj[[-1,n]]=jd; kj[[1,n]]=jd
                Numpn[[-1,Nn,L,jd]]=n
                Numpn[[ 1,Nn,L,jd]]=n
            end
        end
    end 
    infos,izs_ab,nTBME = get_twq_2b(emax,kh,kn,kl,kj,maxsps,jab_max)
    return infos,izs_ab,nTBME
end

function get_twq_2b(emax,kh,kn,kl,kj,maxsps,jab_max)
    infos = [ [0,0,0,0] ];deleteat!(infos,1)
    izs_ab = [ [[0,0,0,0]] ];deleteat!(izs_ab,1)
    ichan = 0
    nTBME = 0
    for izz = -2:2:2
        for ipp = -1:2:1
            for jjx = 0:2*emax+1
                if jjx  > div(jab_max,2);continue;end
                ndim = 0
                tl = [ [0,0,0,0]];deleteat!(tl,1)
                for iza = -1:2:1
                    for izb = iza:2:1
                        if iza + izb != izz;continue;end
                        for ia = 1:maxsps
                            la = kl[[iza, ia]]; ja = kj[[iza, ia]]
                            ibmin = 1
                            if iza == izb;ibmin = ia; end                            
                            for ib = ibmin:maxsps
                                lb = kl[[izb, ib]]; jb = kj[[izb, ib]]                                
                                if (-1)^(la+lb) != ipp;continue;end
                                if tri_check(ja, jb, 2*jjx) == false; continue;end
                                if iza==izb && ia==ib && jjx % 2 == 1;continue;end
                                ndim = ndim + 1
                                push!(tl,[iza,ia,izb,ib])
                            end                            
                        end
                    end
                end
                if ndim != 0
                    ichan += 1
                    #println("ichan $ichan ndim $ndim")
                    push!(infos,[izz,ipp,jjx,ndim])
                    nTBME += Int(ndim*(ndim+1)/2)
                    push!(izs_ab,tl)
                end
            end
        end
    end
    return infos,izs_ab,nTBME
end


function def_sps_snt(emax,target_nlj)
    nljsnt = [ [0,0] ]; deleteat!(nljsnt,1)    
    tzs = [-1,1]
    dict = Dict(0=>0); delete!(dict,0)
    idx = 0; hit = 0
    for pn = 1:2
        tz = tzs[pn]
        for temax = 0:emax
            for l = temax%2:2:temax
                n = div(temax - l,2)
                jmin = 2*l-1
                if jmin < 1; jmin=1;end
                for j=jmin:2:2*l+1
                    push!(nljsnt,[n,l,j,tz])
                    idx += 1
                    for (k,tmp) in enumerate(target_nlj)
                        if tmp[1]==n && tmp[2]==l && tmp[3]==j
                            dict[idx] = k + ifelse(pn==2,length(target_nlj),0)
                        end
                    end
                end
            end
        end
    end
    return nljsnt,dict
end

function freg(p,pp,n)
    if n==0;return 1.0;end
    if n==-1 # for 3D3 pw in older EM interaction
        return 0.5 *( exp( - (p/Lambchi)^4 - (pp/Lambchi)^4 )
                      + exp( - (p/Lambchi)^6 - (pp/Lambchi)^6 ))
    end
    return exp( - (p/Lambchi)^(2*n) - (pp/Lambchi)^(2*n) )
end

"""
    read_LECs!(LECs,idxLECs,dLECs;initialize=false,inpf="src/chiEFTint/LECs.jl")

read LECs from "src/chiEFTint/LECs.jl" (default) and 
store them as ```LECs```(Vector{Float}), ```idxLECs```(Vector{Int}), and ```dLECs```(Dict{str,Float}).
"""
function read_LECs!(LECs,idxLECs,dLECs;initialize=false,inpf="src/chiEFTint/LECs.jl")  
    if initialize
        if !isfile(inpf);inpf="../"*inpf;end # line for test/
        l_E_3F2 = l_E_3F4 = l_E_1F3 = l_e14 = l_e17 = 0.0
        include(inpf)
        if @isdefined(E_3F2); l_E_3F2 = E_3F2;end
        if @isdefined(E_3F4); l_E_3F4 = E_3F4;end
        if @isdefined(E_1F3); l_E_1F3 = E_1F3;end
        if @isdefined(e14); l_e14 = e14;end
        if @isdefined(e17); l_e17 = e17;end
        leclist = [C0_1S0,C0_3S1,C_CSB,C_CIB,
                   C2_3S1,C2_3P0,C2_1P1,C2_3P1,C2_1S0,C2_3SD1,C2_3P2,
                   hD_1S0,D_1S0,D_1P1,D_3P0,D_3P1,D_3P2,hD_3S1,D_3S1,hD_3SD1,D_3SD1,D_3D1,D_1D2,D_3D2,D_3PF2,D_3D3,       
                   l_E_3F2,l_E_1F3,l_E_3F4,
                   c1_NNLO,c2_NNLO,c3_NNLO,c4_NNLO,ct1_NNLO,ct3_NNLO,ct4_NNLO,cD,cE,
                   d12,d3,d5,d145,l_e14,l_e17,
                   c_vs_1,c_vs_2,c_vs_3,c_vs_4,c_vs_5]
        lecname = ["C0_1S0","C0_3S1","C_CSB","C_CIB",
                   "C2_3S1","C2_3P0","C2_1P1","C2_3P1","C2_1S0","C2_3SD1","C2_3P2",
                   "hD_1S0","D_1S0","D_1P1","D_3P0","D_3P1","D_3P2","hD_3S1","D_3S1","hD_3SD1","D_3SD1","D_3D1","D_1D2","D_3D2","D_3PF2","D_3D3",       
                   "E_3F2","E_1F3","E_3F4",
                   "c1_NNLO","c2_NNLO","c3_NNLO","c4_NNLO","ct1_NNLO","ct3_NNLO","ct4_NNLO","cD","cE",
                   "d12","d3","d5","d145","e14","e17",
                   "c_vs_1","c_vs_2","c_vs_3","c_vs_4","c_vs_5"]              
        for (i,lec) in enumerate(leclist)
            tkey = lecname[i]
            dLECs[tkey] = lec
            idxLECs[tkey] = i
            push!(LECs,lec)
        end
    else
        for tkey in keys(idxLECs)
            idx = idxLECs[tkey]
            dLECs[tkey] = LECs[idx]
        end
    end
    return nothing
end

"""
    calc_Vmom!(pnrank,V12mom,tdict,xr,LEC,LEC2,l,lp,S,J,pfunc,n_reg,to;is_3nf=false)    

calc. nn-potential for momentum mesh points
"""
function calc_Vmom!(chiEFTobj,pnrank,V12mom,tdict,xr,LEC,LEC2,l,lp,S,J,pfunc,n_reg,to;is_3nf=false)
    n_mesh = chiEFTobj.n_mesh
    itt = itts[pnrank]; MN = Ms[pnrank]; dwn = 1.0/MN
    V12idx = get(tdict,[itt,lp,l,S,J],-1)
    if V12idx == -1;return nothing;end
    V = V12mom[V12idx]
    @inbounds for i= 1:n_mesh
        x = xr[i]; ex = sqrt(1.0+(dwn.*x)^2)
        @inbounds for j = 1:n_mesh
            y = xr[j]; ey = sqrt(1.0+(dwn.*y)^2)
            ree = 1.0/sqrt(ex*ey)
            if is_3nf;ree=1.0;end
            fac = pfunc(x,y,LEC,LEC2) * freg(x,y,n_reg) *ree
            V[i,j] += fac
        end
    end
    return nothing
end

"""
    Rnl_all_ab(lmax,br,n_mesh,xr_fm)

Returns array for radiul functions (prop to generalized Laguerre polynomials) HO w.f. in momentum space.
Rnlk(l,n,k)=sqrt(br) * R(n,L,Z) *Z with Z=br*k (k=momentum in fm^-1)
"""
function Rnl_all_ab(chiEFTobj,lmax_in,br,n_mesh,xr_fm)
    Nnmax = chiEFTobj.Nnmax
    Rnl = zeros(Float64,Nnmax+1,lmax_in+1,n_mesh)
    for l=0:lmax_in
        for kidx=1:n_mesh
            pb = br * xr_fm[kidx]
            pb2 = pb^2
            fexp = exp(-0.5*pb2)
            fpow = pb^(l+1)
            for n=0:Nnmax
                fac = sqrt(2.0*br*factorial(n) / gamma(n+l+3/2)) * fexp * fpow
                fac *= genLaguerre(n,l+1//2,pb2)
                Rnl[n+1,l+1,kidx]= fac
            end
        end
    end
    return Rnl
end

function bstate(;io=stdout)
    #iz,lz1,lz2,isz,jz
    numst2 = [ [0,0,0,0,0] ];deleteat!(numst2,1)
    dict_numst = [ Dict([0,0]=>1) for i=1:3]
    for i=1:3; delete!(dict_numst[i],[0,0]); end
    arr_numst = [[[[ zeros(Int64,j+iss-abs(j-iss)+1) for ll1=abs(j-iss):j+iss ] for j=0:jmax ] for iss=0:1 ] for pnrank=1:3]
    
    num=0
    ## pp iz = -2
    itt = 1
    for iss=0:1
        for j=0:jmax
            for ll1=abs(j-iss):j+iss
                for ll2=abs(j-iss):j+iss
                    if (-1)^ll1 != (-1)^ll2;continue;end
                    if itt != Int( (1+(-1)^(ll1+iss))/2);continue;end
                    if (-1)^(ll1+iss+itt) != -1;continue;end
                    num=num+1
                    push!(numst2,[-2,ll1,ll2,iss,j])
                    dict_numst[1][[-2,ll1,ll2,iss,j]] = num
                    arr_numst[1][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                end
            end
        end
    end
    ##  nn iz = 2
    pnrank=3
    for iss=0:1
        for j=0:jmax
            for ll1=abs(j-iss):j+iss
                for ll2=abs(j-iss):j+iss
                    if (-1)^ll1 != (-1)^ll2;continue;end
                    if itt != Int( (1+(-1)^(ll1+iss))/2);continue;end
                    if (-1)^(ll1+iss+itt) != -1;continue;end
                    num=num+1
                    push!(numst2,[2,ll1,ll2,iss,j])
                    dict_numst[pnrank][[2,ll1,ll2,iss,j]] = num
                    arr_numst[pnrank][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                end
            end
        end
    end
    ## pn iz = 0
    pnrank=2
    for iss=0:1
        for ktt=1:2
            if (iss==0 && ktt==1) || (iss==1 && ktt==2);itt=1;end
            if (iss==0 && ktt==2) || (iss==1 && ktt==1);itt=0;end
            for j=0:jmax
                for ll1=abs(j-iss):j+iss
                    for ll2=abs(j-iss):j+iss
                        if (-1)^ll1 != (-1)^ll2;continue;end
                        if itt != Int((1+(-1)^(ll1+iss))/2) ;continue;end
                        if (-1)^(ll1+iss+itt) != -1;continue;end
                        num=num+1
                        push!(numst2,[0,ll1,ll2,iss,j])
                        dict_numst[pnrank][[0,ll1,ll2,iss,j]] = num
                        arr_numst[pnrank][iss+1][j+1][ll1-abs(j-iss)+1][ll2-abs(j-iss)+1] = num
                    end
                end
            end
        end
    end
    println(io,"# of two-body states $num")
    return numst2,dict_numst,arr_numst
end

function calc_coulomb(chiEFTobj,xs,ws,Vcoulomb,numst2,nstmax,Rnl;meshp=100)
    hw = chiEFTobj.hw
    Nnmax = chiEFTobj.Nnmax
    rmass = 0.5 * Mp
    brange = sqrt(hc2 / (rmass*hw))
    rmax = 20.0
    r,w = Gauss_Legendre(0.0,rmax,meshp)
    ra1 = zeros(Float64,meshp);rb1 = zeros(Float64,meshp)
    ra2 = zeros(Float64,meshp);rb2 = zeros(Float64,meshp)
    rnl1 = zeros(Float64,meshp)
    rnl2 = zeros(Float64,meshp)
    vcl = zeros(Float64,meshp)
    fmcoul(vcl,r)
    memo1=[0]; memo2=[0]
    for num=1:nstmax
        vcoul = Vcoulomb[num]
        iz12,l1,l2,isz,jj = numst2[num]
        if iz12 !=-2; continue;end
        if l1 !=l2; continue;end
        if jj > jmax;continue;end
        nmax = Nnmax
        memo1[1]=0              
        for n1=0:nmax
            howf(num,brange,n1,l1,r,ra1,rb1,rnl1,meshp,memo1) # => rnl1
            memo2[1]=0
            for n2=0:n1
                howf(num,brange,n2,l2,r,ra2,rb2,rnl2,meshp,memo2) # => rnl2
                scoul = 0.0
                for i=1:meshp
                    scoul += rnl1[i]*rnl2[i]*w[i]*vcl[i]
                end
                vcoul[n1+1,n2+1]= vcoul[n2+1,n1+1]=scoul
            end
        end
    end
    return nothing
end
function fmcoul(vcl,r)
    for (i,x) in enumerate(r)
        vcl[i] = hc *fsalpha /x
    end
    return nothing
end

function howf(num,brange,n,l,r,ra,rb,ret_rnl,meshp,memo)
    #BR,N,L,R,RNL,NMESH
    kmax = l +1
    gam = sqrt(pi)
    for k=1:kmax
        gam *= (0.5+k-1)
    end
    if memo[1] < 1
        for i = 1:meshp
            z = r[i]/brange
            zz = z*z
            zl = 1.5 + l
            ra[i] = sqrt(2.0/(gam*brange)) * z^(l+1) *exp(-0.5*zz)
            rb[i] = ra[i] * (zl-zz)/sqrt(zl)
        end
        memo[1] = 1
    end
    if n < 1; ret_rnl .= ra; return nothing;end
    if n ==1; ret_rnl .= rb; return nothing;end

    if n == memo[1]-1
        ret_rnl .= ra
        return nothing
    end
    if n == memo[1]
        ret_rnl .= rb
        return nothing
    end
    if n > memo[1]
        memo[1] += 1
        a = memo[1]*2+l -0.5
        b = sqrt(memo[1] *(memo[1]+l+0.5))
        c = sqrt( (memo[1]-1)*(memo[1]+l-0.5))
        for i =1:meshp
            zz = (r[i]/brange)^2
            ret_rnl[i] = ((a-zz)*rb[i]-c*ra[i])/b
            ra[i] = rb[i]
            rb[i] = ret_rnl[i]
        end
        while memo[1] < n
            memo[1] += 1
            a = memo[1]*2+l -0.5
            b = sqrt(memo[1] *(memo[1]+0.5))
            c = sqrt( (memo[1]-1)*(memo[1]+l-0.5))
            for i =1:meshp
                zz = (r[i]/brange)^2
                ret_rnl[i] = ((a-zz)*rb[i]-c*ra[i])/b
                ra[i] = rb[i]
                rb[i] = ret_rnl[i]
            end
        end
        return nothing
    end
    if n < memo[1] -1 
        memo[1] = 1
    end
    return nothing
end


function Vrel(chiEFTobj,V12mom,numst2,xr_fm,wr,n_mesh,Rnl,to)    
    Nnmax= chiEFTobj.Nnmax
    nstmax = length(numst2)
    V12ab = [zeros(Float64,Nnmax+1,Nnmax+1) for i=1:nstmax]
    Vcoulomb = [zeros(Float64,Nnmax+1,Nnmax+1) for i=1:nstmax]
    if chiEFTobj.coulomb
        calc_coulomb(chiEFTobj,xr_fm,wr,Vcoulomb,numst2,nstmax,Rnl)
    end
    x = zeros(Float64,Nnmax+1,n_mesh)
    tl = [0,0,0]
    for num = 1:nstmax
        Vtmp = V12mom[num]
        Vab = V12ab[num]
        vcoul = Vcoulomb[num]
        iz,l1,l2,isz,jz = numst2[num]
        @inbounds for n1 = 0:Nnmax
            tR = @views Rnl[n1+1,l1+1,:]
            tx = @views x[n1+1,:]
            @inbounds for k2=1:n_mesh
                sum=0.0
                @inbounds for k1=1:n_mesh
                    vk1k2=Vtmp[k1,k2]
                    if vk1k2 == 0.0; continue;end
                    pkx1=tR[k1]
                    sum += pkx1*vk1k2*xr_fm[k1]*wr[k1]
                end
                #x[n1+1,k2]=sum*wr[k2]*xr_fm[k2]
                tx[k2]=sum*wr[k2]*xr_fm[k2]
            end
        end
        @inbounds for n1=0:Nnmax
            @inbounds for n2=0:Nnmax
                tR = @views Rnl[n2+1,l2+1,:]
                vsum=0.0
                @inbounds for k2=1:n_mesh
                    pkx2=tR[k2]
                    vsum += x[n1+1,k2]*pkx2
                end
                phase=(-1.0)^(n1+n2) 
                t_vcoul = vcoul[n1+1,n2+1]
                Vab[n1+1,n2+1]= phase * vsum + t_vcoul
            end
        end
    end
    return V12ab
end


"""
    hw_formula(A,fnum)

empirical formula for harmonis oscillator parameter hw by mass number A
fnum=2: for sd-shell, Ref. J. Blomqvist and A. Molinari, Nucl. Phys. A106, 545 (1968).
"""
function hw_formula(A,fnum)
    hw = 0.0
    if fnum == 2        
        hw = 45.0 * (A^(-1.0/3.0)) -25.0 * (A^(-2.0/3.0))
    else
        hw = 41.0 * (A^(-1.0/3.0))
    end
    return hw
end


