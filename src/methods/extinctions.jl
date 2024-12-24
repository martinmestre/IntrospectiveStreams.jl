"""Extinction correction for Gaia magnitudes from Gaia dataset."""
function correct_extinction_Gaia(file_orig::String, file_corr::String)::Nothing
    data = at.Table.read(file_orig)
    g = pyia.GaiaData(data)
    bp0 = g.get_BP0()
    rp0 = g.get_RP0()
    g0 = g.get_G0()
    data["bp"] = bp0
    data["rp"] = rp0
    data["g"] = g0
    data.write(file_corr, format="fits", overwrite=true)
    return nothing
end

"""Performing extinction correction for a given list of streams."""
function correct_extinction_Gaia_loop(name_s::Vector{String}, file_orig::Vector{String}, file_corr::Vector{String})
    for i in 1:length(name_s)
        println("Correcting extinction of stream $(name_s[i])")
        correct_extinction_Gaia(file_orig[i], file_corr[i])
    end
    GC.gc()
end