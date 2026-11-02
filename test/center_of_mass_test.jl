# Self-contained oscillator fixtures: no generated nuclear interaction is needed.
function cm_test_space(emax, hw)
    NT = NuclearToolkit
    _, psps, nsps = NT.def_sps(emax)
    sps, dicts1b = NT.make_sps_and_dict_isnt2ims(psps, nsps, emax)
    binfo = NT.basedat(NT.def_nuc("He4", "nucl", ""), "", hw, emax, "nucl")
    dicts = [Dict{Int64,Vector{Vector{Float64}}}() for _ in 1:3]
    _, _, ch1, ch2, _, _ = NT.store_1b2b(sps, dicts1b, dicts, binfo)
    params = NT.init_chiEFTparams(;fn_params="", io=nothing)
    params.Nnmax = 8
    to = NT.TimerOutput()
    dWS = NT.prep_dWS2n(params, to;emax_calc=emax)
    for tbc in ch2.Chan2b
        pnrank = 2 + div(tbc.Tz, 2)
        for (i, bra) in enumerate(tbc.kets), ket in tbc.kets[i:end]
            # kinetic_tb uses proton-first pn pairs; the HF channels are sorted
            # by orbit index. Convert both pairs, including exchange phases.
            abcd = [bra..., ket...]
            phase = 1.0
            for k in (1, 3)
                a, b = abcd[k:k+1]
                if sps[a].tz > sps[b].tz
                    phase *= (-1)^(div(sps[a].j + sps[b].j, 2) + tbc.J + 1)
                    abcd[k], abcd[k+1] = b, a
                end
            end
            orbits = [[sps[a].n, sps[a].l, sps[a].j, sps[a].tz] for a in abcd]
            vpp = hw * phase * NT.kinetic_tb(orbits..., tbc.J, dWS)
            key = NT.get_nkey_from_abcdarr([bra..., ket...])
            push!(get!(dicts[pnrank], key, Vector{Float64}[]),
                  [Float64(tbc.J), vpp / 4, 0.0, 0.0, vpp, 0.0])
        end
    end
    return (;sps, dicts1b, dicts, binfo, ch1, ch2, dWS, to)
end

function cm_test_snt(path, space)
    NT = NuclearToolkit
    (;sps, dicts1b, dicts, binfo) = space
    dim = length(sps) ÷ 2
    open(path, "w") do io
        println(io, "$dim $dim 0 0")
        for i in 1:2dim
            s = sps[dicts1b.snt2ms[i]]
            println(io, "$i $(s.n) $(s.l) $(s.j) $(s.tz)")
        end
        println(io, "0 0") # one-body kinetic energy is constructed by hf_main
        println(io, sum(length(v) for d in dicts for v in values(d)), " 0")
        for d in dicts, (key, entries) in d
            abcd = zeros(Int, 4)
            NT.get_abcdarr_from_intkey!(key, abcd)
            indices = join([dicts1b.ms2snt[a] for a in abcd], " ")
            for v in entries
                # Cancel the intrinsic two-body kinetic term at A=4. This
                # leaves a diagonal oscillator Hamiltonian with a 0s^4 ground
                # state, so adding Lawson must not drive an IMSRG flow.
                println(io, indices, " ", Int(v[1]), " ", -v[5]/4,
                        " 0.0 ", v[5]/binfo.hw)
            end
        end
    end
end

@testset "Center-of-mass Lawson term" begin
    NT = NuclearToolkit
    hw = 20.0
    space = cm_test_space(1, hw)
    (;sps, dicts1b, dicts, binfo, ch1, ch2, dWS, to) = space
    TCM = NT.InitOp(ch1, ch2.Chan2b)
    VCM = NT.InitOp(ch1, ch2.Chan2b)
    NT.CalculateTCM!(TCM, binfo, ch1, ch2.Chan2b, sps)
    NT.Calculate_RCM(binfo, ch1, ch2.Chan2b, sps, VCM, dWS, to;non0_ij=false)
    NT.aOp!(VCM, 0.5 * NT.Mm * hw^2 / NT.hc2)

    @testset "He4 emax=1 reference" begin
        # Tcm = Vcm = 3/4 hw for 0s^4; subtracting 3/2 hw gives zero.
        kinetic = sum(2 * TCM.onebody[pn][1,1] for pn in 1:2)
        potential = sum(2 * binfo.nuc.A * VCM.onebody[pn][1,1] for pn in 1:2)
        @test kinetic ≈ 0.75hw atol=1.e-12
        @test potential ≈ 0.75hw atol=1.e-12
        @test kinetic + potential - 1.5hw ≈ 0.0 atol=1.e-12
    end

    @testset "Complete two-body term and beta scaling" begin
        H0, _, _, _, _, _ = NT.store_1b2b(sps, dicts1b, dicts, binfo)
        for beta in (0.0, 1.0, 2.0)
            HCM = deepcopy(VCM)
            NT.aOp!(HCM, beta)
            target = deepcopy(dicts)
            NT.update_dicts_withHCM!(HCM, ch2, target;BetaCM=beta,A=binfo.nuc.A)
            NT.recalc_v!(binfo.nuc.A, target)
            H, _, _, _, _, _ = NT.store_1b2b(sps, dicts1b, target, binfo)
            for (ch, tbc) in enumerate(ch2.Chan2b)
                cm = H.twobody[ch] - H0.twobody[ch]
                # For beta=1, the CM kinetic term cancels the intrinsic one.
                expected = beta * (binfo.nuc.A * VCM.twobody[ch] - H0.twobody[ch])
                @test cm ≈ expected atol=1.e-12
                for (i, bra) in enumerate(tbc.kets), (j, ket) in enumerate(tbc.kets)
                    nb = sum(2sps[a].n + sps[a].l for a in bra)
                    nk = sum(2sps[a].n + sps[a].l for a in ket)
                    if nb != nk
                        @test cm[i,j] ≈ 0.0 atol=1.e-12
                    end
                end
            end
            # Replacing the CM term must not accumulate it or change Vpp.
            previous = deepcopy(target)
            NT.update_dicts_withHCM!(HCM, ch2, target;BetaCM=beta,A=binfo.nuc.A)
            @test target == previous
        end
    end

    @testset "Radial off-diagonal kinetic term and mass scaling" begin
        radial = cm_test_space(2, hw)
        T = NT.InitOp(radial.ch1, radial.ch2.Chan2b)
        NT.CalculateTCM!(T, radial.binfo, radial.ch1, radial.ch2.Chan2b, radial.sps)
        for pn in 1:2, a in pn:2:length(radial.sps), b in pn:2:length(radial.sps)
            oa, ob = radial.sps[a], radial.sps[b]
            expected = hw / 4 * NT.kinetic_ob([oa.n, oa.l, oa.j], [ob.n, ob.l, ob.j])
            @test T.onebody[pn][(a+1)÷2,(b+1)÷2] ≈ expected atol=1.e-12
        end
        NT.difA_TCM(T, 4, 16)
        oxygen = NT.basedat(NT.def_nuc("O16", "nucl", ""), "", hw, 2, "nucl")
        T16 = NT.InitOp(radial.ch1, radial.ch2.Chan2b)
        NT.CalculateTCM!(T16, oxygen, radial.ch1, radial.ch2.Chan2b, radial.sps)
        @test T.onebody == T16.onebody

        V16 = deepcopy(VCM)
        NT.difA_RCM(V16, 4, 16)
        target = deepcopy(dicts)
        NT.update_dicts_withHCM!(V16, ch2, target;BetaCM=1.0,A=16)
        NT.recalc_v!(16, target)
        for (ch, tbc) in enumerate(ch2.Chan2b)
            for (i, bra) in enumerate(tbc.kets), (j, ket) in enumerate(tbc.kets)
                i > j && continue
                key = NT.get_nkey_from_abcdarr([bra..., ket...])
                v = only(filter(v -> v[1] == tbc.J, target[2+div(tbc.Tz,2)][key]))
                # At beta=1 the two-body kinetic terms cancel at any A.
                @test v[2] ≈ 16 * V16.twobody[ch][i,j] atol=1.e-12
            end
        end

        @testset "r1 dot r2 agrees with direct angular-momentum coupling" begin
            for tbc in radial.ch2.Chan2b
                for bra in tbc.kets, ket in tbc.kets
                    abcd = [bra..., ket...]
                    phase = 1.0
                    for k in (1, 3)
                        a, b = abcd[k:k+1]
                        if radial.sps[a].tz > radial.sps[b].tz
                            phase *= (-1)^(div(radial.sps[a].j + radial.sps[b].j, 2) + tbc.J + 1)
                            abcd[k], abcd[k+1] = b, a
                        end
                    end
                    os = [[radial.sps[a].n, radial.sps[a].l, radial.sps[a].j, radial.sps[a].tz] for a in abcd]
                    normalization = sqrt((1 + (bra[1] == bra[2])) * (1 + (ket[1] == ket[2])))
                    # hcm_2body uses reduced r and nabla matrix elements,
                    # independently of the HOB transformation used by RCM.
                    direct = NT.hcm_2body(vcat(os...)..., tbc.J) / normalization
                    direct += NT.kinetic_tb(os..., tbc.J, radial.dWS)
                    hob = NT.calc_single_r1r2(4, bra, ket, radial.sps, tbc.J, radial.dWS, 1.0, radial.to)
                    @test hob ≈ phase * direct atol=1.e-12
                end
            end
        end
    end

    @testset "hf_main honors fn_params in HF and IMSRG" begin
        mktempdir(@__DIR__) do dir
            snt = joinpath(dir, "oscillator.snt")
            cm_test_snt(snt, space)
            # A conflicting default file detects ignoring fn_params.
            write(joinpath(dir, "optional_parameters.jl"), "BetaCM = 7.0\nemax = 0\n")
            cd(dir) do
                for beta in (0.0, 1.0, 2.0)
                    params = joinpath(dir, "lawson_$beta.jl")
                    write(params, "BetaCM = $beta\nemax = 0\nNnmax = 8\nsmax = 2.0\ndsmax = 0.5\ndenominatorDelta = 0.0\n")
                    hf = hf_main(["He4"], snt, hw, 1;fn_params=params,return_obj=true)
                    @test hf.E0 ≈ 45.0 atol=1.e-10
                    @test hf.EMP2 ≈ 0.0 atol=1.e-12
                    imsrg = hf_main(["He4"], snt, hw, 1;fn_params=params,
                                    return_obj=true,doIMSRG=true,delete_Ops=true)
                    @test imsrg.H.zerobody[1] ≈ 45.0 atol=1.e-10
                end
            end
        end
    end
end
