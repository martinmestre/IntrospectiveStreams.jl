"""Example of isocrhone interpolation"""
function example_interpolate_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R)::DataFrame where {T<:Real,R<:Real}
    file_plot = "plots/examples/isochrone_$(family)_$(photsys).pdf"
    df = interpolate_isochrone(family, photsys, age, metal)
    df_ezp = interpolate_isochrone(family, photsys, age, metal; ezbool=true)
    df_dld = download_isochrone(family, photsys, age, metal)
    # isochrone_df.color = isochrone_df.gmag - isochrone_df.rmag
    # plot_isochrone_cmd(isochrone_df, photsys, file_plot)
    download_df = download_isochrone(family, photsys, 10^log_age, metal)
    download_df.color = download_df.gmag - download_df.rmag
    download_df.label = string.(download_df.label)
    plot_isochrone_cmd(download_df, photsys, file_plot_download)
    println("Isócrona completa")
    return nothing
end