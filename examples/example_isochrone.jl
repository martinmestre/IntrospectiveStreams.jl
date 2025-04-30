"""Example of isocrhone interpolation"""
function example_interpolate_isochrone()
    family = :parsec
    photsys= :hsc
    file_plot = "plots/examples/isochrone_$(family)_$(photsys).pdf"
    file_plot_download = "plots/examples/isochrone_$(family)_$(photsys)_download.pdf"
    log_age = log10(12e9) # log10(age)
    metal =  z2feh_parsec(0.001) # Metallicity
    @show metal
    file_artif = "artifacts/isochrones/parsec/$(photsys)/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
    df_artif = read_parsec_file(file_artif)
    @show df_artif
    # Generate isochrone
    # isochrone_df = interpolate_isochrone(family, photsys, log_age, metal; ezpadova_bool=true)
    # isochrone_df.color = isochrone_df.gmag - isochrone_df.rmag
    # plot_isochrone_cmd(isochrone_df, photsys, file_plot)
    download_df = download_isochrone(family, photsys, 10^log_age, metal)
    download_df.color = download_df.gmag - download_df.rmag
    download_df.label = string.(download_df.label)
    plot_isochrone_cmd(download_df, photsys, file_plot_download)
    println("Isócrona completa")
    return nothing
end