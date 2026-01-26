using ModelingToolkit, Catalyst, ReactionNetworkImporters
using SIAN, Nemo, StructuralIdentifiability, Groebner
using LinearAlgebra, OffsetArrays, ProgressMeter, Logging, Random, Base.Iterators, Primes, UnicodePlots

include("../make_poly_from_ode.jl")
include("../plug_numbers_in_poly.jl")
include("../quotient_basis.jl")

struct AllMatrices
    n::Int
    m::Int
    range
end

function Base.iterate(S::AllMatrices, state=1)
    if sum(state) > length(S)
        return nothing # Signals the end of iteration
    else
        # Returns a tuple of the current value and the next state
        return (S[state], state + 1)
    end
end

function Base.getindex(S::AllMatrices, i::Int)
    i = i - 1
    A1 = digits(i, base=length(S.range), pad=S.n*S.m)
    A = collect(reshape(A1, S.n, S.m))
end

Base.length(S::AllMatrices) = length(S.range)^(S.n*S.m)
Base.eltype(::Type{AllMatrices}) = Matrix{Int}

struct IteratorProduct{T1,T2}
    I1::T1
    I2::T2
    IteratorProduct(I1::T1, I2::T2) where {T1,T2} = new{T1,T2}(I1,I2)
end

function Base.iterate(S::IteratorProduct, state=1)
    if state > length(S)
        return nothing # Signals the end of iteration
    else
        # Returns a tuple of the current value and the next state
        return (S[state], state + 1)
    end
end

function Base.getindex(S::IteratorProduct, i::Int)
    i1 = div(i - 1, length(S.I2)) + 1
    i2 = i - (i1 - 1)*length(S.I2)
    return (S.I1[i1], S.I2[i2])
end

Base.length(S::IteratorProduct) = length(S.I1)*length(S.I2)
Base.eltype(::Type{IteratorProduct{T1,T2}}) where {T1,T2} = Tuple{eltype(T1),eltype(T2)}

begin
    p = IteratorProduct(1:10, 'a':'z')
    @assert length(p) == 260
    @assert eltype(p) == Tuple{Int,Char}
    @assert collect(p) == [(i, j) for i in 1:10 for j in 'a':'z']
end


function is_correct_stoichiometry(two_mat)
    ok(m) = all(col -> sum(col) > 0, eachcol(m))
    return ok(two_mat[1]) && ok(two_mat[2])
end

function handle_example(X, K_val, K, Y, two_mat, n_observed, index)
    mn = MatrixNetwork(K_val, two_mat[1], two_mat[2]; species = X, params = K)
    prn = loadrxnetwork(mn; name = Symbol(:testnetwork_, index))
    crn = complete(prn.rn)
    
    ode = convert(ODESystem, crn)
    measured_quantities = [
      Y[i] ~ ModelingToolkit.unknowns(ode)[i]
      for i in 1:n_observed
    ]

    res = ode_to_poly(ode, measured_quantities)
    poly = res["polynomial_system"]
    states = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(X)), gens(parent(poly[1])))
    parameters = filter(s -> any(s2 -> startswith(string(s), chopsuffix(string(s2), "(t)")), vcat(K)), gens(parent(poly[1])))
    states = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in states]
    parameters = [StructuralIdentifiability.parent_ring_change(f, parent(poly[1])) for f in parameters]
    parameters = Vector{eltype(states)}(parameters)
    poly_num, subs = plug_numbers_in_poly(poly, states, parameters)

    r_drl, _ = polynomial_ring(base_ring(parent(poly_num[1])), symbols(parent(poly_num[1])), internal_ordering=:degrevlex)
    poly_num = map(f -> StructuralIdentifiability.parent_ring_change(f, r_drl), poly_num)
    poly_num_zp = map(f -> map_coefficients(c -> GF(2^30+3)(c), f), poly_num)
    gb = groebner(poly_num_zp);
    qb = quotient_basis(gb)
    
    crn, measured_quantities, res, poly_num, qb
end

function run_many_times(
                n_species, 
                m_reactions;
                range = 0:1,
                n_observed = 1,
                symbolic_rate = [true for _ in 1:m_reactions],
                how_many = :all # :all or a number
    )
    @eval t = default_t()
    
    str = join(["X"*string(i)*"(t)" for i in 1:n_species], " ")
    expr = "@species $str"
    X = eval(Meta.parse(expr))

    str = join(["k"*string(i)*"" for i in 1:m_reactions], " ")
    expr = "@parameters $str"
    K = eval(Meta.parse(expr))

    str = join(["y"*string(i)*"(t)" for i in 1:n_species], " ")
    expr = "@variables $str"
    Y = eval(Meta.parse(expr))

    Random.seed!(42)
    K_val = copy(K)
    for (i, issym) in enumerate(symbolic_rate)
        if !issym
            K_val[i] = nextprime(2, i+3)
        end
    end
    
    data = []
    two_matrices = IteratorProduct(
        AllMatrices(n_species, m_reactions, range),
        AllMatrices(n_species, m_reactions, range)
    )
    @info "There are $(length(two_matrices)) combinations to consider."

    if how_many == :all
    @showprogress enabled=true showspeed=true for (i, mat) in enumerate(two_matrices)
        if is_correct_stoichiometry(mat)
        result = handle_example(X, K_val, K, Y, mat, n_observed, i)
        push!(data, [mat, result])
        end
    end
    else
    index = rand(1:length(two_matrices), how_many)
    @showprogress enabled=true showspeed=true for idx in index
        mat = two_matrices[idx]
        if is_correct_stoichiometry(mat)
        result = handle_example(X, K_val, K, Y, mat, n_observed, idx)
        push!(data, [mat, result])
        end
    end
    end

    data 
end

#=
data = run_many_times(2,3,range=0:1);
histogram(map(x -> length(x[2][end]), data) |> sort, nbins=4, height=8, width=14)

q4 = data[findall(==(4), map(x -> length(x[2][end]), data))]
all(x -> x[2][3]["identifiability_nemo"]["nonidentifiable"] |> isempty, q4)

q3 = data[findall(==(3), map(x -> length(x[2][end]), data))]
all(x -> x[2][3]["identifiability_nemo"]["nonidentifiable"] |> isempty, q3)

data = run_many_times(2,3,range=0:1, symbolic_rate=[false,false,false]);
=#
