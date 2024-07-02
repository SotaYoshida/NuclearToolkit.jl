struct ksl_state
    J2::Int    
    prty::Int
    T::String
    Energy::Float64
end

struct kshell_nuc
    nuc::nucleus
    sntf::String
    Egs::Float64
    Egs_target::Float64
    states::Vector{ksl_state}
end

"""
    read_kshell_summary(fns::Vector{String};targetJpi="",nuc="")

Reading KSHELL summary file from specified paths, `fns`.
In some beta-decay studies, multiple candidates for parent g.s. will be considered.
In that case, please specify optional argument `targetJpi` e.g., targetJpi="Si36", otherwise the state havbing the lowest energy in summary file(s) will be regarded as the parent ground state.
Example:
```
Energy levels

N    J prty N_Jp    T     E(MeV)  Ex(MeV)  log-file

1     0 +     1     6   -319.906    0.000  log_Ar48_SDPFSDG_j0p.txt 
```
!!! note
    J is doubled in old versions of KSHELL (kshell_ui.py).
"""
function read_kshell_summary(fns::Vector{String};targetJpi="",nuc="")
    Egs = Egs_target = 1.e+5; sntf=""
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
            Egs = min(Energy,Egs)
            if targetJpi == "" 
                Egs_target = Egs
            end
            if targetJpi == tmpJpi
                Egs_target = min(Energy,Egs_target)
            end            
        end
    end
    @assert Egs_target != 1.e+5 "Egs_target error"
    nuclei = def_nuc(string(nuc),"","")
    return kshell_nuc(nuclei,sntf,Egs,Egs_target,states)
end
