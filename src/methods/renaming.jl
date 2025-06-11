function rename_mag!(df::DataFrame, mag::Symbol; pattern::String="")
    col = Symbol(names(df)[startswith.(names(df),string(mag)*pattern)][1])
    rename!(df, col => mag)
    return
end

function rename_mag!(df::DataFrame, mags::Vector{Symbol}; pattern::String="")
    for mag ∈ mags
        rename_mag!(df, mag, pattern=pattern)
    end
end
function rename_mag!(df::DataFrame, mag::Symbol, pattern::String)
    col = Symbol(string(mag)*pattern)
    rename!(df, col => mag)
    return
end
function rename_mag!(df::DataFrame, mags::Vector{Symbol}, pattern::String)
    for mag ∈ mags
        rename_mag!(df, mag, pattern)
    end
end