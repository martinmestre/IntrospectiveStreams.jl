const base_dir = "/home/mmestre/casa/work/data"
const file_gc  = "$(base_dir)/gc/Baumgardt/gc.txt"

"""Gaia file names for CATS project."""
function name_files_Gaia(name_s::String)
    host_dir = "$(base_dir)/cats/Gaia"
    file_orig = "$(host_dir)/orig/GaiaDR3-$(name_s)-all.fits"
    file_corr = "$(host_dir)/corr/GaiaDR3-$(name_s)-all_extincorr.fits"
    file_filt = "$(host_dir)/filt/GaiaDR3-$(name_s)-all_filt.fits"
    return [file_orig, file_corr, file_filt]
end

"""Photometry file names for CATS project."""
function name_files_photometry(name_s::String)
    host_dir = "$(base_dir)/cats/phot"
    if name_s ∈ Set(["GD-1", "Pal5", "PS1-A"])
        file_phot = "$(host_dir)/PS1DR2-$(name_s)_xm.fits"
    elseif name_s ∈ Set(["Jhelum"])
        file_phot = "$(host_dir)/DES-$(name_s)_xm.fits"
    elseif name_s ∈ Set(["Fjorm-M68"])
        file_phot = nothing
    end
    return [file_phot]
end

"""Gaia file names for CATS project."""
function name_files_isochrone(name_s::String, age::Number, metal::Number)
    host_dir = "$(base_dir)/cats/isochrone"
    file_iso = "$(host_dir)/isochrone_$(name_s)_age$(age)_metal$(metal).csv"
    return [file_iso]
end

"""Plot file names."""
function name_files_plots(name_s::String)
    host_dir = "$(base_dir)/cats/plots"
    file_scatter = "$(host_dir)/scatter_sky_self_frame_$(name_s).pdf"
    return [file_scatter]
end

"""File names wrapper for all of them."""
function name_files_all(name_s::Vector{String}, ages::Vector{<:Number}, metals::Vector{<:Number})
    file_orig = Vector{Union{Nothing, String}}(nothing, length(name_s))
    file_corr = Vector{Union{Nothing, String}}(nothing, length(name_s))
    file_phot = Vector{Union{Nothing, String}}(nothing, length(name_s))
    file_iso  = Vector{Union{Nothing, String}}(nothing, length(name_s))
    file_filt = Vector{Union{Nothing, String}}(nothing, length(name_s))
    file_plot = Vector{Union{Nothing, String}}(nothing, length(name_s))

    for i in eachindex(name_s)
        v = name_files_Gaia(name_s[i]) # v might have more than 3 component.
        file_orig[i] = v[1]
        file_corr[i] = v[2]
        file_filt[i] = v[3]
    end

    for i in eachindex(name_s)
        u = name_files_photometry(name_s[i])
        v = name_files_isochrone(name_s[i], ages[i], metals[i])
        w = name_files_plots(name_s[i])
        file_phot[i] = u[1]
        file_iso[i] = v[1]
        file_plot[i] = w[1]
    end


    return file_orig, file_corr, file_phot, file_iso, file_filt, file_plot
end