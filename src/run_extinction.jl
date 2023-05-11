"""Run extinction correction for all the streams."""

name_s = ["GD-1", "Pal5", "PS1-A", "Jhelum", "Fjorm-M68"]
file_orig = Vector{Union{Nothing, String}}(nothing, length(name_s))
file_corr = Vector{Union{Nothing, String}}(nothing, length(name_s))
file_phot = Vector{Union{Nothing, String}}(nothing, length(name_s))


for i in eachindex(name_s)
    v = name_files_Gaia_cats(name_s[i]) # v might have more than 2 components.
    file_orig[i] = v[1]
    file_corr[i] = v[2]
end

for i in eachindex(name_s)
    v = name_files_photometry_cats(name_s[i]) # v might have more than 2 components.
    file_phot[i] = v[1]
end


# correct_extinction_Gaia_loop(name_s, file_orig, file_corr)
println("Extinction correction performed.")