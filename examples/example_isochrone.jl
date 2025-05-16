"""Example of isocrhone interpolation"""
function example_interpolate_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R) where {T<:Real,R<:Real}
    df = interpolate_isochrone(family, photsys, age, metal)
    df_dld = download_isochrone(family, photsys, age, metal)
    df_ezp = interpolate_isochrone(family, photsys, age, metal; ezbool=true)
    df_array = [df, df_dld, df_ezp]
    categorize_phases!.(df_array, :parsec)
    colorear!.(df_array, :gmag, :rmag)
    fig = plot_isochrone_cmd(df_array, :parsec, :hsc, :gmag, :color_gr; only = collect(1:7))
    display(fig)
    println("🔭 Isochrones plotted")
    return nothing
end