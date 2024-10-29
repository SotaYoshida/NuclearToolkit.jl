#using LocalCoverage;  generate_coverage("NuclearToolkit"; run_test=true)
using NuclearToolkit
using Test
using Printf

@testset "NuclearToolkit.jl" begin
    include("chiEFTint_test.jl")
    include("HFMBPT_IMSRG_test.jl")
    include("ShellModel_test.jl")

    rm("tdmat",recursive=true)
    rm("wavsamples",recursive=true)
    rm("appwavs",recursive=true)
    rm("He6_ckpot_j0.wav")
    rm("He6_ckpot_j4.wav")
    rm("interaction_file/ckpot.msnt")

    include("read_3BME_test.jl")

    include("threebody_test.jl")
    rm("A3files",recursive=true)

    include("caliblation_test.jl")

    rm("tbme_emn500n4lo_2n3n_srg10.0hw20emax2.snt.bin")
    rm("tbme_em500n3lo_barehw20emax2.snt")
    rm("tbme_em500n3lo_barehw20emax2.snt.bin")
    rm("vsimsrg_p-shell_coreHe4refHe4_He4_hw20e2_Delta0.0.snt")
    rm("flowOmega",recursive=true)
    rm("pic",recursive=true)

end

