struct ksl_state
    J2::Int    
    prty::Int
    T::String
    Energy::Float64
end

struct kshell_nuc
    nuc::nuclei
    sntf::String
    Egs::Float64
    states::Vector{ksl_state}
end

"""

Example:
>
>Energy levels
>
>N    J prty N_Jp    T     E(MeV)  Ex(MeV)  log-file
>
>1     0 +     1     6   -319.906    0.000  log_Ar48_SDPFSDG_j0p.txt 

"""
function read_kshell_summary(fns::Vector{String};targetJpi="",mode="",nuc="")
    Egs = 1.e+5; sntf=""
    states=ksl_state[]
    for fn in fns
        f = open(fn,"r")
        lines = readlines(f) 
        close(f)
        use2J_insummary = false
        for line in lines
            if occursin("log-file",line) && occursin("2J",line); use2J_insummary = true; end
            if !occursin("log_",line);continue;end
            tl = split(line)
            @assert length(tl) == 8 "unexpected number of elements $(length(tl)) in summary"
            if nuc == ""
                nuc = split(tl[end],"_")[2]            
            else
                @assert occursin(nuc,tl[end]) "nuc = $nuc must be in 'log-file' column"
            end
            if sntf == ""; sntf = split(tl[end],"_")[3];end
            J =tl[2]
            J2 = string(J)
            if !use2J_insummary
                if occursin("/",J2)
                    J2 = parse(Int,split(J2,"/")[1])
                else
                    J2 = 2* parse(Int,J2)
                end
            else
                J2 = parse(Int,J2)
            end
            prty = ifelse(tl[3]=="+",1,-1)
            NJp = parse(Int,tl[4])
            T = tl[5]
            Energy = parse(Float64,tl[6])
            state = ksl_state(J2,prty,T,Energy)
            if !(state in states)
                push!(states,state)
            end 
            tmpJpi = "j"*string(J2)*ifelse(tl[3]=="+","p","n") 
            if targetJpi == "" || targetJpi == tmpJpi
                if Energy < Egs ; Egs = Energy; end
            end
        end
    end
    nuclei = def_nuc(string(nuc),"","")
    return kshell_nuc(nuclei,sntf,Egs,states)
end
