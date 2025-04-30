"""
    download_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R; age_scale::String="linear")::DataFrame where {T<:Real,R<:Real}
Download a single stellar isochrone.
[age] = Gyr.
"""
function download_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R; age_scale::String="linear")::DataFrame where {T<:Real,R<:Real}
    age_yr = 1e9*age
    if(family==:mist)
        @assert -4≤metal≤0.5 "Metallicity should satisfy: -4 ≤ FeH ≤ 0.5 (FeH≈MH)."
        @assert 5≤log10(age_yr)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
        println("Note that MIST uses metallicity [FeH] (not Z abundance).")
        df = ezmist.get_one_isochrone(age=age_yr, FeH=metal, v_div_vcrit=0.0,
                    age_scale=age_scale, output_option="photometry",
                    output=string(photsys), Av_value=0.0).to_pandas()|> PyPandasDataFrame |> DataFrame
    elseif(family==:parsec)
        @assert -2.2<metal≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
        @assert 0≤age≤13.5 "Age [Gyr] should fulfill: 0 ≤  age [Gyr] ≤ 13.5"
        println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        df = ezpadova.get_isochrones(age_yr=(age_yr, age_yr, 0), MH=(metal, metal, 0),
                        model="parsec12s", photsys_file=string(photsys))|> PyPandasDataFrame |> DataFrame
    end
    return df
end

"""
    download_isochrone(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R})::DataFrame
Download grid of stellar isochrones.
[age] = Gyr.
"""
function download_isochrone(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R})::DataFrame where {T<:Real, R<:Real}
    age_yr = 1e9*age
    @assert family==:parsec
    @assert -2.2<metal[1] && metal[2]≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
    @assert 0≤age[1] && age[2]≤13.5 "Age [Gyr] should fulfill: 0 ≤  age [Gyr] ≤ 13.5"
    println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        df = ezpadova.get_isochrones(age_yr=age_yr, MH=metal,model="parsec12s",
        photsys_file=string(photsys))|> PyPandasDataFrame |> DataFrame
    return df
end


"""
    build_isochrone_grid(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R}) where {T<:Real, R<:Real}
Build and save a large grid of isochrones in (age, FeH) plane to be later used by fuction interpolate_isochrone.
[age] = Gyr.
"""
function build_isochrone_grid(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R}) where {T<:Real, R<:Real}
    dir_path = joinpath("artifacts", "isochrones", string(family), string(photsys))
    age_str = "age_$(round(age[1], digits=1))to$(round(age[2], digits=1))"
    metal_str = "metal_$(round(metal[1], digits=2))to$(round(metal[2], digits=2))"
    filename = "$(age_str)_$(metal_str).jld2"
    file_artif = joinpath(dir_path, filename)

    mkpath(dirname(file_artif))
    jldopen(file_artif, "a+", compress=true) do file
        for a in range(age[1], age[2], step=age[3])
            key_age = "age=$(round(a, digits=1))"
            df = download_isochrone(family, photsys, (a, a, 0), metal)
            metal_groups = groupby(df, :MH)
            for sub_df in metal_groups
                key_MH = "MH=$(round(sub_df.MH[1], digits=2))"
                key = "$key_age/$key_MH"
                file[key] = sub_df
            end
            println("✅ Completed processing for isochrone of $a Gyr")
        end
    end
end


"""
    interpolate_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R; ezpadova_bool::Bool=false)::DataFrame where {T<:Real,R<:Real}Interpolate isochrones from the previously downloaded data base (artifacts dir).
If using 'ezpadova_bool=true' then Python's ezpadova.QuickInterpolator() is used.
If using default value then Julia searches the closest isochrone in (age, FeH) plane.
The corresponding grid of isochrones used in both cases was donwloaded and saved with the
function build_isochrone_grid().
[age] = Gyr.
 """
function interpolate_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R; ezpadova_bool::Bool=false)::DataFrame where {T<:Real,R<:Real}
        @assert family == :parsec "Only Parsec isochrones accepted for the moment"
        @assert 0≤age≤13.5 "Age [Gyr] should fulfill: 0 ≤  age [Gyr] ≤ 13.5"
        @assert -2.2<metal≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."

        if ezpadova_bool
            file_artif = "artifacts/isochrones/$(family)/$(photsys)/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
            quickiso =  ezpadova.QuickInterpolator(file_artif)
            df = quickiso(log(age*1e9), metal) |> PyPandasDataFrame |> DataFrame
            df.label .=  string.(Int.(floor.(df.evol))) # Recompute label so as not to have inerpolated value
            return df
        else
            file_artif = "artifacts/isochrones/$(family)/$(photsys)/family_age=$(age[1]):$(age[2])Gyr_MH=$(metal[1]):$(metal[2]).jld2"
            df_artif = read_parsec_file(file_artif)
            return find_closest_isochrone(df_artif, age, metal)
        end
end




