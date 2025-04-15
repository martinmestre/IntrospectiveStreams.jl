"""Example of isocrhone interpolation"""
function example_interpolate_isochrone()
    family = :parsec
    photsys= :hsc
    file_plot = "plots/examples/isochrone_$(family)_$(photsys).pdf"
    log_age = 9.5 # log10(age)
    metal = 0.0  # Metallicity

    # Generate isochrone
    isochrone_df = interpolate_isochrone(family, photsys, log_age, metal; ezpadova_bool=true)
    isochrone_df.color = isochrone_df.gmag - isochrone_df.rmag
    plot_isochrone_cmd(isochrone_df, photsys, file_plot)
    println("Isócrona completa")
    return isochrone_df
end