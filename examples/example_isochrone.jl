"""Example of isocrhone interpolation"""
function example_interpolate_isochrone()
    family = :parsec
    photsys= :hsc
    file_plot = "plots/examples/isochrone_$(family)_$(photsys).pdf"
    file_plot_download = "plots/examples/isochrone_$(family)_$(photsys)_download.pdf"
    log_age = 9.2 # log10(age)
    metal = 0.0  # Metallicity

    # Generate isochrone
    isochrone_df = interpolate_isochrone(family, photsys, log_age, metal; ezpadova_bool=true)
    isochrone_df.color = isochrone_df.gmag - isochrone_df.rmag
    plot_isochrone_cmd(isochrone_df, photsys, file_plot)
    download_df = download_isochrone(family, photsys, 10^log_age, metal)
    download_df.color = download_df.gmag - download_df.rmag
    download_df.label = string.(download_df.label)
    plot_isochrone_cmd(download_df, photsys, file_plot_download)
    println("Isócrona completa")
    return nothing
end