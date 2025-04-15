"""Download stellar isochrones"""
function download_isochrone(family::Symbol, photsys::String, age::T, metal::R; age_scale::String="linear")::DataFrame where {T<:Real,R<:Real}
    if(family==:mist)
        @assert 5≤log10(age)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
        @assert -4≤metal≤0.5 "Metallicity should satisfy: -4 ≤ FeH ≤ 0.5 (FeH≈MH)."
        println("Note that MIST uses metallicity [FeH] (not Z abundance).")
        df = ezmist.get_one_isochrone(age=age, FeH=metal, v_div_vcrit=0.0,
                    age_scale=age_scale, output_option="photometry",
                    output=photsys, Av_value=0.0).to_pandas()|> PyPandasDataFrame |> DataFrame
    elseif(family==:parsec)
        @assert 5≤log10(age)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
        @assert -2.2<metal≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
        println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        df = ezpadova.get_isochrones(age_yr=(age,age,0), MH=(metal,metal,0),
                        model="parsec12s", photsys_file=photsys)|> PyPandasDataFrame |> DataFrame
    end
    return df
end

function download_isochrone(family::Symbol, photsys::String, age::NTuple{3,Number}, metal::NTuple{3,Number})::DataFrame
    @assert family==:parsec
    @assert  5≤log10(age[1]) && log10(age[2])≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
    @assert -2.2<metal[1] && metal[2]≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
    println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        df = ezpadova.get_isochrones(age_yr=age, MH=metal,model="parsec12s",
        photsys_file=photsys)|> PyPandasDataFrame |> DataFrame
    return df
end

"Interpolate isochrones from the previously downloaded data base (artifacts dir)"
function interpolate_isochrone(family::Symbol, photsys::Symbol, log_age::T, metal::R; ezpadova_bool::Bool=false)::DataFrame where {T<:Real,R<:Real}
        @assert family == :parsec "Only Parsec isochrones accepted for the moment"
        if(family==:parsec)
                @assert 9.2≤log_age≤10.3 "Age should fulfill: 9.2 ≤ log10(age) ≤ 10.3."
                @assert -2.2<metal≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
                if(photsys==:hsc)
                        file_artif = "artifacts/isochrones/parsec/$(photsys)/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
                        if ezpadova_bool
                            quickiso =  ezpadova.QuickInterpolator(file_artif)
                            return quickiso(log_age, metal) |> PyPandasDataFrame |> DataFrame
                        else
                            df_artif = read_parsec_file(file_artif)
                            return interpolate_isochrone(df_artif, log_age, metal)
                        end
                end
        end
end




