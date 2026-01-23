function quotient_basis(
    poly
)
    R = parent(poly[1])
    @assert internal_ordering(R) == :degrevlex
    present_vars = union(map(vars, poly)...)
    leading_exponents = map(leading_monomial, poly)
    exponents_to_check = Set()
    exponents_checked = Set()
    basis_exponents = Set()
    push!(exponents_to_check, one(R))
    while length(exponents_to_check) > 0
        e = pop!(exponents_to_check)
        push!(exponents_checked, e)
        if !any(le -> divides(e, le)[1], leading_exponents)
            push!(basis_exponents, e)
            for i in 1:length(present_vars)
                next_e = e * present_vars[i]
                if !(next_e in exponents_checked) && !(next_e in exponents_to_check)
                    push!(exponents_to_check, next_e)
                end
            end
        end
    end
    basis_exponents
end
