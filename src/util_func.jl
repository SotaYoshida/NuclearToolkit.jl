const reg = r"[0-9]+"

function latex_nuc(nuc::String)
    A =  match(r"\d+", nuc).match
    el = replace(nuc, r"\d+" => "")
    lnuc = latexstring("{}^{$(A)} \\mathrm{$el}") 
    return lnuc
end

function int2bitstr(x, N)
    return string(bitstring(x)[end-N+1:end] |> s -> lpad(s, N, '0'))
end

const element = Dict{Int, String}()
element[0] = "n"
el_list = ["H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"]
for (Z, e) in enumerate(el_list)
    element[Z] = e
end

function rm_comment(lines)
    nlines = []
    for line in lines
        line = strip(line)
        if length(line) > 0
            if startswith(line,"!")||startswith(line,"#")
                continue
            end
        end
        push!(nlines,line)
    end
    return nlines
end

function rm_nan(array)
    na = []
    for tmp in array
        if tmp != "";push!(na,tmp); end
    end
    return na
end

function print_vec(s,v,io=stdout;ine=false,long=false)
    s *= " "
    for i = 1:length(v)
        if ine
            if long 
                s *= @sprintf "%11.3e" v[i]
            else
                s *= @sprintf "%9.1e" v[i]
            end
        else
            if long
                s *= @sprintf "%15.8f" v[i]
            else  
                s *= @sprintf "%10.4f" v[i]
            end
    	end
    end
    println(io,s)
end

function get_nkey2_arr(tmp;ofst=10^3)
    return tmp[1] + tmp[2] * ofst
end

function get_nkey2(i,j;ofst=10^3)
    return i + j * ofst
end

function get_nkey2_u(i,j)::UInt64
    return (UInt64(i) << 10) +  UInt64(j)
end

function get_nkey3(i,j,k;ofst=10^3)
    return i + ofst * j + ofst^2 * k
end

function get_nkey3_u(i,j,k)::UInt64
    return (UInt64(i) << 20) + (UInt64(j) << 10) +  UInt64(k)
end

function get_nkey3_JPT(arr::Vector{Int64})::UInt64
    return (UInt64(arr[1]+2) << 20) + (UInt64(arr[2]+2) << 10) +  UInt64(arr[3])
end
function get_nkey3_JPT(aTz::Int64,P::Int64,J::Int64)::UInt64
    return (UInt64(aTz+2) << 20) + (UInt64(P+2) << 10) +  UInt64(J)
end

function get_nkey3_ketJ(i,j,J)::UInt64
    return (UInt64(i) << 20) + (UInt64(j) << 10) +  UInt64(J)
end

function get_nkey4(i,j,k,l;ofst=10^3)
    return i + ofst * j + ofst^2 * k + ofst^3 * l
end

function get_nkey6(j1::Int64,j2::Int64,j3::Int64,j4::Int64,j5::Int64,j6::Int64)::UInt64
    return  (UInt64(j1) << 50) + (UInt64(j2) << 40) +(UInt64(j3) << 30) +  (UInt64(j4) << 20) + (UInt64(j5) << 10) +  UInt64(j6)
end

"""
Unlike other wigner symbols, constructing CG coefficients encounter the case with negative indices.
One way to avoid this is to use the following function to get hash key for the CG coefficients.
This function asssume all `j` is >= -3
"""
function get_nkey6_shift(j1::Int64,j2::Int64,j3::Int64,j4::Int64,j5::Int64,j6::Int64; int_shift=3)::UInt64    
    return (UInt64(j1+int_shift) << 50) + (UInt64(j2+int_shift) << 40) +(UInt64(j3+int_shift) << 30) +  (UInt64(j4+int_shift) << 20) + (UInt64(j5+int_shift) << 10) +  UInt64(j6+int_shift)
end

"""
    get_nkey_from_abcdarr(tkey;ofst=1000)

To get integer key from an Int array (with length greater than equal 4)
"""
function get_nkey_from_abcdarr(tkey;ofst=1000)
    return tkey[1] + tkey[2] * ofst + tkey[3] * ofst^2 + tkey[4] * ofst^3 
end

function show_allitems_of_dict(d, title::String="Dictionary Items")
    println("=== $(title) ===")
    println("# Items: $(length(d))")
    for (k,v) in d
        println("$(k) => $(v)")
    end
end

function deltaf(i::Int64,j::Int64)
    ifelse(i==j,1.0,0.0)
end

function show_matrix(text, mat)
    println("$text:")
    for i = 1:size(mat)[1]
        print_vec("",mat[i,:])
    end
    return nothing
end