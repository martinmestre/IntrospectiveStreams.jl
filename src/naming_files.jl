"""File names for CATS project."""
function name_files_Gaia_cats(name_s::String)
    host_dir = "/home/mmestre/casa/work/data/cats"
    file_orig = "$(host_dir)/GaiaDR3-$(name_s)-all.fits"
    file_corr = "$(host_dir)/GaiaDR3-$(name_s)-all_extincorr.fits"
    return [file_orig, file_corr]
end

function name_files_photometry_cats(name_s::String)
    host_dir = "/home/mmestre/casa/work/data/cats"
    if name_s ∈ Set(["GD-1", "Pal5", "PS1-A"])
        file_phot = "$(host_dir)/PS1DR2-$(name_s)_xm.fits"
    elseif name_s ∈ Set(["Jhelum"])
        file_phot = "$(host_dir)/DES-$(name_s)_xm.fits"
    elseif name_s ∈ Set(["Fjorm-M68"])
        file_phot = nothing
    end
    return [file_phot]
end