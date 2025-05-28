function print_counts_per_mag_bin(df::DataFrame, mag::Symbol, k::StepRange{T,T}) where {T<:Real}
    for i ∈ k
        println("N($i<i<$(i+1))=$(sum(i.< df[!,mag] .< i+1))")
    end
end