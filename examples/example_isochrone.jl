"""Example of isocrhone interpolation"""
function example_interpolate_isochrone()
    filter="YBC_hsc"
    file_artif = "artifacts/isochrones/parsec/$filter/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
    df_artif = DataFrame(CSV.File(file_artif, delim=" ", ignorerepeated=true, comment="#"))
    target_age = 10.0^10.2 # age (linear)
    target_metallicity = 0.0  # Metalicidad deseada

    # Ver columnas con missing values
    columnas_con_missing(df_artif)

    eliminar_columnas_con_missing!(df_artif)

    # Generar la isócrona completa
    isochrone_df = interpolate_isochrone(df_artif, target_age, target_metallicity)
    println("Isócrona completa:")
    println(isochrone_df)
    return isochrone_df
end