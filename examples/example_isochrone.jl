"""Example of isochrone interpolation"""
function example_interpolate_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R) where {T<:Real,R<:Real}
    df = interpolate_isochrone(family, photsys, age, metal)
    # df_dld = download_isochrone(family, photsys, age, metal)
    df_ezp = interpolate_isochrone(family, photsys, age, metal; ezbool=true)
    df_array = [df, df_ezp]
    algos = [:istreams,  :ezpadova]
    categorize_phases!.(df_array, :parsec)
    colorear!.(df_array, :gmag, :rmag)
    parstring = @sprintf("_age%.1f_MH%+.2f", age, metal)
    fig, filename = plot_isochrone_cmd(df_array, algos, :parsec, :hsc, :gmag, :color_gr; only = collect(1:7),paramstring=parstring)
    display(fig)
    dir = joinpath("plots","examples")
    isdir(dir) || mkpath(dir)
    save(joinpath(dir,filename), fig, pt_per_unit=1; pdf_version="1.5")
    println("🔭 Isochrones plotted")
    return nothing
end


