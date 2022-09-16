"""
Not used
"""

struct dWS3n
  dcgm0::Dict{Vector{Int64},Float64}
  d6j::Dict{Vector{Int64},Float64}
  d9j::Dict{Vector{Int64},Float64}
  dtri::Dict{Vector{Int64},Float64}
  keycg::Vector{Int64}
  key6j::Vector{Int64}
  key9j::Vector{Int64}
end 

struct params3b
    e3max::Int64
    N3max::Int64
    L3max::Int64
    j3max::Int64
    hw::Int64
    TFrenorm::Bool
    meshpoints::Vector{Vector{Float64}}
    LECs::Vector{Float64}
    dWS::dWS3n
end

function test3NF(;target_LECs=["c1_NNLO","c3_NNLO","c4_NNLO","cD","cE"])
    to = TimerOutput()
    LECs = Float64[ ];idxLECs=Dict{String,Int64}();dLECs=Dict{String,Float64}()
    read_LECs!(LECs,idxLECs,dLECs;initialize=true)
    for (k,target) in enumerate(target_LECs)
        idx = idxLECs[target]
        tLEC = LECs[idx]
        dLECs[target] = tLEC 
    end
    ### parameters 
    TFrenorm = false
    Np = Nq = 25; pmax3 = 8.0
    e3max = 8
    L3max = 0 # to be 4
    N3max = 8 # to be 40?
    j3max = 9 # said to be sufficient, but can be e3max*2 + 3
    gl_p,glw_pq = Gauss_Legendre(0.0,pmax3,Np)
    gl_q = copy(gl_p)
    meshpoints = [gl_p,gl_q,glw_pq]
    @timeit to "precalc 6j&9j" dWS = prep_dWS3n(e3max,N3max)
    params = params3b(e3max,N3max,L3max,j3max,hw,TFrenorm,meshpoints,LECs,dWS)

    labframe = true
    if labframe
        #calculation of 3NF in laboratory frame
        @timeit to "ch" calc_channel_T(params,to)
    else
        #solve three-body system
    end
    show(to, allocations = true,compact = false);println("")
    #base3NF = params3b(e3max,hw,TFrenorm)
end

function calc_channel_T(params,to)
    j3max = params.j3max
    JPTs = Vector{Int64}[ ]
    n = 0
    for t = 1:2:3
        for j = 1:2:j3max
            for parity = 1:-2:-1
                n += 1
                push!(JPTs, [j,parity,t])
            end 
        end
    end    
    nch = (j3max +1) * 2
    if n!= nch;println("warn!");end
    cfpdims = Vector{Vector{Int64}}[]
    for ich = 1:nch
        JPT = JPTs[ich]      
        j,p,t = JPT 
        println("ich ",@sprintf("%4i",ich),
        "   J =",@sprintf("%3i", j),
        "   P =",@sprintf("%3i", p),
        "   T =",@sprintf("%3i", t))
        push!(cfpdims, Vector{Int64}[])
        if !check_binexists(ich,JPT,params) 
          @timeit to "get" get_dim_cfp!(cfpdims[ich],JPT,params,to)
          @timeit to "set" set_cfps(cfpdims[ich],JPT,params,to)
        else
           #read_cfps()
        end
    end 
    return nothing
end

function make_cfpbinname(JPT,params;dirname="./cfps/")
  j,p,t=JPT
  Nmax = params.N3max
  fname = dirname*"cfp_j"*string(j)*"p"
  fname *= ifelse(p==1,"+","-")*"t"*string(t)
  fname *= "_Nmax"*string(Nmax)*".bin"
  return fname
end 

function check_binexists(ich,JPT,params;dirname="./cfps/")
  fname = make_cfpbinname(JPT,params;dirname=dirname)
  return isfile(fname)
end


function get_dim_cfp!(cfpdim,JPT,params,to)  
  dWS = params.dWS
  j,p,t = JPT 
  for N = 0:params.N3max
    ndim = 0
    diag = Float64[]
    for n12 = 0:div(N,2)
      for n3 = 0:div(N-2*n12,2)
        for l12 = 0:N-2*n12 -2*n3
          l3 = N- 2*n12 -2*n3 -l12
          if (-1)^(l12+l3) != p;continue;end
          for s12 = 0:1
            for j12 = abs(l12-s12):l12+s12
              for t12 = 0:1
                if !tri_check(2*t12,1,t);continue;end
                if (-1)^(l12+s12+t12) != -1 ;continue;end
                for j3 = abs(2*l3-1):2:2*l3 +1
                    if !tri_check(2*j12, j3, j);continue;end
                    ndim += 1
                    if N==0 && JPT == [1,-1,1]
                      println(" n12 $n12,l12 $l12, s12 $s12 ,",
                      " j12 $j12, t12 $t12, n3 $n3, l3 $l3, j3 $j3 ")
                    end
                    r = anti_op_isospin(dWS,n12,l12,s12,j12,t12,n3,l3,j3,
                                        n12,l12,s12,j12,t12,n3,l3,j3,j,t,
                                        to)
                    push!(diag,r)
                end
              end
            end
          end
        end
      end
    end
    north = Int(round(sum(diag)))
    push!(cfpdim,[ndim,north])
  end  
  return nothing
end

function set_cfps(cfpdim,JPT,params,to)
  fname = make_cfpbinname(JPT,params)
  io = open(fname,"w")
  j,p,t= JPT; Nmax=params.N3max
  write(io, j);write(io,p);write(io,t);write(io,Nmax)
  dWS = params.dWS
  nlsjts = [ zeros(Int64,8,cfpdim[N+1][1]) for N=0:params.N3max]
  j,p,t = JPT 
  for N = 0:params.N3max   
    nphys, north = cfpdim[N+1]
    nlsjt = @views nlsjts[N+1]
    i = 0
    for n12 = 0:div(N,2)
      for n3 = 0:div(N-2*n12,2)
        for l12 = 0:N-2*n12 -2*n3
          l3 = N- 2*n12 -2*n3 -l12
          if (-1)^(l12+l3) != p;continue;end
           for s12 = 0:1
             for j12 = abs(l12-s12):l12+s12
              for t12 = 0:1
                if !tri_check(2*t12,1,t);continue;end
                if (-1)^(l12+s12+t12) != -1 ;continue;end
                for j3 = abs(2*l3-1):2:2*l3 +1
                    if !tri_check(2*j12, j3, j);continue;end
                    i += 1
                    #println("i $i nlsjt $nlsjt l12 $l12 l3 $l3 j3 $j3")
                    nlsjt[1,i] = n12; nlsjt[2,i] = l12
                    nlsjt[3,i] = s12; nlsjt[4,i] = j12
                    nlsjt[5,i] = t12; nlsjt[6,i] = n3
                    nlsjt[7,i] = l3; nlsjt[8,i] = j3                
                end
              end
            end
          end
        end
      end
    end    
    #println("JPT $JPT i $i ndim/north $nphys $north ")
    if nphys == 0;continue;end
    if north == 0;continue;end    
    A = zeros(Float64,nphys,nphys)   
    for ib = 1:nphys # bra
        n12,l12,s12,j12,t12,n3,l3,j3 = @view nlsjt[:,ib]
        for ik = 1:ib # ket
          n45,l45,s45,j45,t45,n6,l6,j6= @view nlsjt[:,ik]
          @timeit to "anti" anti = anti_op_isospin(dWS,
                                  n12,l12,s12,j12,t12,n3,l3,j3,
                                  n45,l45,s45,j45,t45,n6,l6,j6,j,t,to)
          A[ib,ik] = A[ik,ib] = anti
        end
    end     
    @timeit to "eigen" vals,vecs = eigen(A)

    write_cfp_bin(io,north,nphys,vals,vecs,JPT,N,params)
  end
  close(io)
  return nothing
end

function write_cfp_bin(io,north,nphys,vals,vecs,JPT,N,params)
  write(io, nphys);write(io,north);write(io,N)
  cfp = zeros(Float64,north,nphys)
  hit = 0
  for i in eachindex(vals)
    val = vals[i]
    vec = @view vecs[:,i]
    if abs(val - 1.0) < 1.e-5
      hit +=1
      cfp[hit,:] .= vec
      write(io,vec)
    elseif abs(val-0.0) < 1.e-5
      nothing
    else
      println("warn! something is wrong");exit()
    end
  end
  if north != hit; println("warn! something is wrong for north");exit();end
  return nothing
end

"""
  prep_dWS3n(e3max,N3max)

  to prepare Wigner symbols (CG/6j/9j) for three body force
"""
function prep_dWS3n(e3max,N3max)
    dcgm0 = Dict{Vector{Int64},Float64}()
    dtri = Dict{Vector{Int64},Float64}()
    d6j = Dict{Vector{Int64},Float64}()
    d9j = Dict{Vector{Int64},Float64}()
    for l = 0:2*N3max
      for l1 = 0:2*N3max
        for l2 = 0:2*N3max
          key = [l1,l2,l]
          dtri[key] = trinomial(l1,l2,l)
        end
      end
    end
    for l = 0:N3max
      for l1 = 0:N3max
        for l2 = abs(l-l1):l+l1
          if !tri_check(l1,l2,l);continue;end
          key = [l1,l2,l]
          dcgm0[key] = clebschgordan(Float64,l1,0,l2,0,l,0)
        end
      end
    end
    s3 = 1//2
    for s12=0:1
      for s45=0:1
        smin = max(abs(s12-s3),abs(s45-s3))
        smax = min(s12+s3,s45+s3)
        for S = smin:smax
          key = [ 2*s12, 2*S,2*s45] 
          d6j[key] = wigner6j(Float64,1//2,1//2,s12,1//2,S,s45)
        end
      end
    end
    for l12 = 0:N3max
      for s12 = 0:1
        for j12 = abs(l12-s12):l12+s12
          for l3 = 0:N3max-l12
            s3 = 1//2
            for j3 = abs(l3-s3):l3+s3
              for Lam = abs(l12-l3):l12+l3
                if !tri_check(Lam,l12,l3);continue;end
                for S = abs(s12-s3):s12+s3
                  if !tri_check(S,s12,s3);continue;end
                  for J = abs(j12-j3):j12+j3
                    if !tri_check(J,j12,j3);continue;end
                    key = [l12,2*s12,2*j12,
                           l3,       j3*2,
                           Lam,2*S,2*J] # s3=1/2 is trivial, so skipped
                    d9j[key] = wigner9j(l12,s12,j12,l3,s3,j3,Lam,S,J)
                  end
                end
              end
            end
          end
        end
      end
    end
    #wigner9j(l12,s12,j12,l3,1//2,j3//2,lambda,stot//2,jtot//2)
    #wigner6j(Float64,1//2,1//2,s12,1//2,stot//2,s45)
    keycg = zeros(Float64,3); key6j = zeros(Float64,3);key9j = zeros(Float64,8)
    return dWS3n(dcgm0,d6j,d9j,dtri,keycg,key6j,key9j)
end 

"""
  overwrite_key6!(s12,S,s45,key)

overwrite key for 6j symbols
"""
function overwrite_key6!(s12,S,s45,key)
  key[1]=s12; key[2]=S; key[3]=s45
  return nothing
end

"""
  overwrite_key9!(l12,s12,j12,l3,j3,Lam,stot,jtot,key)

overwrite key for 9j symbols
"""
function overwrite_key9!(l12,s12,j12,l3,j3,Lam,stot,jtot,key)
  key[1]=l12; key[2]=s12; key[3]=j12
  key[4]=l3 ; key[5]=j3;  
  key[6]=Lam; key[7]=stot; key[8]=jtot
  return nothing
end

"""
    anti_op_isospin(dWS,n12,l12,s12,j12,t12,n3,l3,j3,n45,l45,s45,j45,t45,n6,l6,j6,jtot,ttot,to)
  
Function to calc. matrix element of antisymmetrizer.  
Detailed derivation can be found in e.g., Eq.(3.119) of Master Thesis by Joachim Langhammer (2010), TUDarmstadt.
"""
function anti_op_isospin(dWS,
                         n12,l12,s12,j12,t12,n3,l3,j3,
                         n45,l45,s45,j45,t45,n6,l6,j6,jtot,ttot,
                         to)
  ex = 0.0
  if (2 * n12 + l12 + 2 * n3 + l3 != 2 * n45 + l45 + 2 * n6 + l6);return 0.0;end
  s_min = max(abs(2*s12 -1), abs(2*s45 -1))
  s_max = min(2*s12+1, 2*s45+1)
  X6 = dWS.d6j; key6 = dWS.key6j
  X9 = dWS.d9j; key9 = dWS.key9j
  for stot = s_min:2:s_max #Note: stot is already doubled
    lambda_min = max(abs(l12-l3),abs(l45-l6),div(abs(jtot-stot),2))
    lambda_max = min(l12+l3, l45+l6, div(jtot+stot,2))
      for lambda = lambda_min:lambda_max      
        tmp  = sqrt((2*j12+1)*(j3+1)*(2*lambda+1)*(stot+1)) 
        tmp *= sqrt((2*j45+1)*(j6+1)*(2*lambda+1)*(stot+1))
        tmp *= sqrt((2*s12+1)*(2*s45+1))     
        try
          overwrite_key9!(l12,2*s12,2*j12,l3,j3,lambda,stot,jtot,key9); tmp *= X9[key9]
        catch 
          println("key9 $key9 l12 $l12 l3 $l3 j3 $j3")
          exit()
        end
        overwrite_key9!(l45,2*s45,2*j45,l6,j6,lambda,stot,jtot,key9); tmp *= X9[key9]      
        overwrite_key6!(2*s12,stot,2*s45,key6);tmp *= X6[key6]
        tmp *= gmosh2(n12, l12, n3, l3, n45, l45, n6, l6, lambda,1.0/3.0,dWS,to)
        ex += tmp
      end
  end
  tmp = (-1.0)^(s12 + t12 + s45 + t45) 
  tmp *= sqrt((2*t12+1)*(2*t45+1)) 
  overwrite_key6!(2*t12,ttot,2*t45,key6)
  tmp *= X6[key6]
  ex *= tmp
  anti = - 2.0 * ex / 3.0
  if (n12 == n45 && l12 == l45 && s12 == s45 && j12 == j45 && 
      t12 == t45 && n3 == n6 && l3 == l6 && j3 == j6) 
      anti = (1.0 - 2.0 * ex ) / 3.0
  end 
  return anti
end