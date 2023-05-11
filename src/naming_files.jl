"""Gaia file names for CATS project."""
function name_files_Gaia_cats(name_s::String)
    host_dir = "/home/mmestre/casa/work/data/cats/Gaia"
    file_orig = "$(host_dir)/orig/GaiaDR3-$(name_s)-all.fits"
    file_corr = "$(host_dir)/corr/GaiaDR3-$(name_s)-all_extincorr.fits"
    return [file_orig, file_corr]
end

"""Photometry file names for CATS project."""
function name_files_photometry_cats(name_s::String)
    host_dir = "/home/mmestre/casa/work/data/cats/phot"
    if name_s ∈ Set(["GD-1", "Pal5", "PS1-A"])
        file_phot = "$(host_dir)/PS1DR2-$(name_s)_xm.fits"
    elseif name_s ∈ Set(["Jhelum"])
        file_phot = "$(host_dir)/DES-$(name_s)_xm.fits"
    elseif name_s ∈ Set(["Fjorm-M68"])
        file_phot = nothing
    end
    return [file_phot]
end

"""File names wrapper for all of them."""
function name_files_all(name_s::Vector{String})
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

    return file_orig, file_corr, file_phot
end