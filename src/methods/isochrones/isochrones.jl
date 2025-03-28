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



function interpolate_isochrone(df::DataFrame, target_age::Real, target_metallicity::Real)
    # Convert age to logarithmic scale
    log_age = log10(target_age)

    # Get 4 nearest isochrones (already filtered from df)
    nearest_isochrones = find_nearest_isochrones(df, log_age, target_metallicity)

    # Verify we have exactly 4 unique age-metallicity combinations
    unique_pairs = unique(nearest_isochrones[:, [:logAge, :MH]])
    if size(unique_pairs, 1) != 4
        error("Expected exactly 4 unique age-metallicity combinations. Found: $(size(unique_pairs, 1))")
    end

    # Get initial mass ranges from nearest_isochrones (not full df)
    initial_mass_ranges = [unique(nearest_isochrones[(nearest_isochrones.logAge .== row.logAge) .&
                              (nearest_isochrones.MH .== row.MH), :Mini]) for row in eachrow(unique_pairs)]

    # Calculate safe intersection range
    min_masses = [minimum(iso) for iso in initial_mass_ranges]
    max_masses = [maximum(iso) for iso in initial_mass_ranges]
    intersection_range = (maximum(min_masses), minimum(max_masses))

    # Select reference isochrone (with most points in intersection range)
    iso_with_most_masses = argmax([count(m -> intersection_range[1] <= m <= intersection_range[2], iso) for iso in initial_mass_ranges])

    # Filter masses strictly within intersection range
    initial_masses = filter(m -> intersection_range[1] < m < intersection_range[2],
                          initial_mass_ranges[iso_with_most_masses])

    # Prepare output DataFrame
    interpolated_df = DataFrame(initial_mass=initial_masses)

    # Unique nodes for bilinear interpolation (sorted)
    age_nodes = sort(unique(nearest_isochrones.logAge))
    metallicity_nodes = sort(unique(nearest_isochrones.MH))

    # Interpolate each relevant column
    for col in setdiff(names(df), ["logAge", "MH", "Mini"])
        interpolators = Dict()

        # Build interpolators for each unique isochrone using nearest_isochrones
        for row in eachrow(unique_pairs)
            iso = nearest_isochrones[(nearest_isochrones.logAge .== row.logAge) .&
                                  (nearest_isochrones.MH .== row.MH), [:Mini, Symbol(col)]]

            # Clean and sort (with physical tolerance)
            clean_iso = unique(iso, :Mini) |> x -> sort(x, :Mini)

            # Verify sufficient points
            if nrow(clean_iso) < 2
                error("Insufficient points for $col (logAge=$(row.logAge), MH=$(row.MH))")
            end

            # Strict interpolator (no extrapolation)
            interpolators[(row.logAge, row.MH)] = LinearInterpolation(
                clean_iso.Mini,
                clean_iso[!, col],
                extrapolation_bc=Throw()  # Will error if extrapolation attempted
            )
        end

        # Calculate magnitudes for each initial mass
        magnitudes = Float64[]
        for m in initial_masses
            try
                # Build 2x2 magnitude matrix
                mags_matrix = [interpolators[(age, metal)](m) for age in age_nodes, metal in metallicity_nodes]

                # Strict bilinear interpolation
                bilin_itp = LinearInterpolation(
                    (age_nodes, metallicity_nodes),
                    mags_matrix,
                    extrapolation_bc=Throw()
                )
                push!(magnitudes, bilin_itp(log_age, target_metallicity))
            catch e
                error("""Interpolation failed for:
                      Column: $col
                      Mass: $m
                      Target age: $log_age
                      Target metallicity: $target_metallicity
                      Error: $e""")
            end
        end

        interpolated_df[!, col] = magnitudes
    end

    return interpolated_df
end