"""Download stellar isochrones"""
function get_isochrone(family::Symbol, filter::String, age::Number, metal::Number;
                        age_scale::String="linear")::DataFrame
    if(family==:mist)
        @assert 5≤log10(age)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
        @assert -4≤metal≤0.5 "Metallicity should satisfy: -4 ≤ FeH ≤ 0.5 (FeH≈MH)."
        println("Note that MIST uses metallicity [FeH] (not Z abundance).")
        df = ezmist.get_one_isochrone(age=age, FeH=metal, v_div_vcrit=0.0,
                    age_scale=age_scale, output_option="photometry",
                    output=filter, Av_value=0.0).to_pandas()|> PyPandasDataFrame |> DataFrame
    elseif(family==:parsec)
        @assert 5≤log10(age)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
        @assert -2.2<metal≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
        println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        df = ezpadova.get_isochrones(age_yr=(age,age,0), MH=(metal,metal,0),
                        model="parsec12s", phot=filter)|> PyPandasDataFrame |> DataFrame
    end
    return df
end
function get_isochrone(family::Symbol, filter::String, age::NTuple{3,Number}, metal::NTuple{3,Number})::DataFrame
    @assert family==:parsec
    @assert  5≤log10(age[1]) && log10(age[2])≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
    @assert -2.2<metal[1] && metal[2]≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
    println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        df = ezpadova.get_isochrones(age_yr=age, MH=metal,model="parsec12s", phot=filter)|> PyPandasDataFrame |> DataFrame
    return df
end


"Interpolate isochrones from the previously downloaded data base (artifacts dir)"
function interpolate_isochrone(family, filter, age, metal)
        @assert family == :parsec "Only Parsec isochrones accepted for the moment"
        if(family==:parsec)
                @assert 5≤log10(age)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
                @assert -2.2<metal≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
                if(filter=="YBC_hsc")
                        file_artif = "artifacts/isochrones/parsec/$filter/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
                        df_artif = DataFrame(CSV.File(file_artif, delim=" ", ignorerepeated=true, comment="#"))
                        return interpolate_isochrone(df_artif, age, metal)
                end
        end
end

# Función para generar la isócrona completa
function interpolate_isochrone(df::DataFrame, target_age::Real, target_metallicity::Real)
    # Encontrar las 4 isócronas más cercanas
    log_age = log10(target_age)
    nearest_isochrones = find_nearest_isochrones(df, log_age, target_metallicity)

    # Obtener los rangos de masa inicial de las isócronas cercanas
    initial_mass_ranges = [unique(df[(df.logAge .== row.logAge) .& (df.MH .== row.MH), :Mini]) for row in eachrow(nearest_isochrones)]

    # Calcular el rango de intersección
    min_masses = [minimum(iso) for iso in initial_mass_ranges]
    max_masses = [maximum(iso) for iso in initial_mass_ranges]
    intersection_range = (maximum(min_masses), minimum(max_masses))  # [max(mínimos), min(máximos)]

    # Seleccionar la isócrona con más puntos dentro del rango de intersección
    iso_with_most_masses = argmax([count(m -> intersection_range[1] <= m <= intersection_range[2], iso) for iso in initial_mass_ranges])

    # Filtrar las masas iniciales de la isócrona seleccionada que están dentro del rango de intersección
    initial_masses = filter(m -> intersection_range[1] <= m <= intersection_range[2], initial_mass_ranges[iso_with_most_masses])


    # Interpolar todos los campos del DataFrame (excepto age, metallicity, initial_mass)
    interpolated_df = DataFrame(initial_mass=initial_masses)
    for col in names(df)
        if col ∉ ["logAge", "MH", "Mini"]
            # Precalcular interpoladores para cada fila de nearest_isochrones
            interpolators = []
            for row in eachrow(nearest_isochrones)
                iso = df[(df.logAge .== row.logAge) .& (df.MH .== row.MH), :]
                sort!(iso, :Mini)  # Ordenar por masa inicial (Mini)
                itp = LinearInterpolation(iso.Mini, iso.magnitude, extrapolation_bc=Line())
                push!(interpolators, itp)
            end

            # Calcular magnitudes para cada initial_mass
            magnitudes = []
            for initial_mass in initial_masses
                mags = []
                for (i, row) in enumerate(eachrow(nearest_isochrones))
                    mag = interpolators[i](initial_mass)  # Usar el interpolador precalculado
                    push!(mags, mag)
                end

                # Interpolación bilineal en edad y metalicidad usando Interpolations.jl
                age_nodes = unique(nearest_isochrones.logAge)
                metallicity_nodes = unique(nearest_isochrones.MH)
                itp = LinearInterpolation((age_nodes, metallicity_nodes), reshape(mags, (2, 2)), extrapolation_bc=Line())
                magnitude = itp(target_age, target_metallicity)
                push!(magnitudes, magnitude)
            end

            interpolated_df[!, col] = magnitudes
        end
    end

    return interpolated_df
end