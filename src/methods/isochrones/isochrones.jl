"""
    download_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R; age_scale::String="linear")::DataFrame where {T<:Real,R<:Real}
Download a single stellar isochrone.
[age] = Gyr.
"""
function download_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R; age_scale::String="linear")::DataFrame where {T<:Real,R<:Real}
    age_yr = 1e9*age
    if(family==:mist)
        @assert -4≤metal≤0.5 "Metallicity should satisfy: -4 ≤ FeH ≤ 0.5 (FeH≈MH)."
        @assert 5≤log10(age_yr)≤10.3 "Age should fulfill: 5 ≤ log10(age_yr) ≤ 10.3."
        println("Note that MIST uses metallicity [FeH] (not Z abundance).")
        return ezmist.get_one_isochrone(age=age_yr, FeH=metal, v_div_vcrit=0.0,
                    age_scale=age_scale, output_option="photometry",
                    output=string(photsys), Av_value=0.0).to_pandas()|> PyPandasDataFrame |> DataFrame
    elseif(family==:parsec)
        @assert -2.19≤metal≤0.5 "Metallicity should satisfy: -2.19 ≤ FeH ≤ 0.5 (FeH≈MH)."
        @assert 0≤age≤14.0 "Age [Gyr] should fulfill: 0 ≤  age [Gyr] ≤ 14.0"
        println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        return ezpadova.get_isochrones(age_yr=(age_yr, age_yr, 0), MH=(metal, metal, 0),
                        model="parsec12s", photsys_file=string(photsys))|> PyPandasDataFrame |> DataFrame
    end
end

"""
    download_isochrone(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R})::DataFrame
Download grid of stellar isochrones.
[age] = Gyr.
"""
function download_isochrone(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R})::DataFrame where {T<:Real, R<:Real}
    age_yr = 1e9 .*age
    @assert family==:parsec
    @assert -2.19≤metal[1] && metal[2]≤0.5 "Metallicity should satisfy: -2.19 ≤ FeH ≤ 0.5 (FeH≈MH)."
    @assert 1e-9≤age[1] && age[2]≤14.0 "Age [Gyr] should fulfill: 1e-9 ≤  age [Gyr] ≤ 14.0"
    println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
    return ezpadova.get_isochrones(age_yr=age_yr, MH=metal,model="parsec12s", photsys_file=string(photsys))|> PyPandasDataFrame |> DataFrame
    return
end


"""
    build_isochrone_grid(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R}) where {T<:Real, R<:Real}
Build and save a large grid of isochrones in (age, FeH) plane to be later used by fuction interpolate_isochrone. Having all the grid in one single file is not possible because the downloading
may be corrupted in some point.
[age] = Gyr.
Multipledispatch option for building database in chunks of age range.
"""
function build_isochrone_grid(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R}) where {T<:Real, R<:Real}
    dirpath = joinpath(PKGDIR, "artifacts", "isochrones", string(family), string(photsys))
    age_str = @sprintf("ageGyr_%.1fto%.1f", age[1], age[2])
    metal_str = @sprintf("metal_%+.2fto%+.2f", metal[1], metal[2])
    filename = "$(age_str)_$(metal_str).jld2"
    file_artif = joinpath(dirpath, filename)

    mkpath(dirname(file_artif))
    for a in range(age[1], age[2], step=age[3])
        key_age = ""
        jldopen(file_artif, "a+"; compress=true) do file
            key_age = @sprintf("age=%0.1f", a)
            df = download_isochrone(family, photsys, (a, a, 0.0), metal)
            size_bytes = Base.summarysize(df)
            println("DataFrame size: $(size_bytes/1024^2) MB")
            println("Now gouping by metallicity")
            metal_groups = groupby(df, :MH)
            for sub_df_ in metal_groups
                sub_df = copy(sub_df_)
                size_bytes = Base.summarysize(sub_df)
                key_MH = @sprintf("MH=%+.2f", sub_df.MH[1])
                key = "$key_age/$key_MH"
                file[key] = sub_df
            end
        end
        println("Current file size: ", filesize(file_artif) / 1024^2, " MB")
        println("$key_age saved")
    end
    println("✅ Completed process for building isochrone database for $family family
    and $photsys photometry for age = $age.")
    return
end
function build_isochrone_grid(family::Symbol, photsys::Symbol, age::NTuple{3,T}, metal::NTuple{3,R}, chunk_ids::StepRange{I,I}) where {T<:Real, R<:Real, I<:Integer}
    age_nodes = collect(age[1]:age[3]:age[2])
    k = diff(chunk_ids)[1]-1
    for i ∈ collect(chunk_ids)
        age_chunk = (age_nodes[i-k], age_nodes[i], age[3])
        build_isochrone_grid(family, photsys, age_chunk, metal)
    end
    println("🔭 Completed full process for building isochrone database for $family family
    and $photsys photometry in age chunks associated to different files.")
    return
end

"""
    interpolate_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R; ezpadova_bool::Bool=false)::DataFrame where {T<:Real,R<:Real}Interpolate isochrones from the previously downloaded data base (artifacts dir).
If using 'ezpadova_bool=true' then Python's ezpadova.QuickInterpolator() is used.
If using default value then Julia searches the closest isochrone in (age, FeH) plane.
The corresponding grid of isochrones used in both cases was donwloaded and saved with the
function build_isochrone_grid().
[age] = Gyr.
 """
function interpolate_isochrone(family::Symbol, photsys::Symbol, age::T, metal::R; ezbool::Bool=false)::DataFrame where {T<:Real,R<:Real}
    @assert family == :parsec "Only Parsec isochrones accepted for the moment"
    @assert -2.19<metal≤0.5 "Metallicity should satisfy: -2.19 < FeH ≤ 0.5 (FeH≈MH)."
    dirpath = joinpath(PKGDIR,"artifacts", "isochrones", string(family), string(photsys))
    if !ezbool
        @assert 0≤age≤14.0 "Age [Gyr] should fulfill: 0 ≤  age [Gyr] ≤ 14.0"
        df, key = find_nearest_isochrone(dirpath, age, metal)
        println("✅ Isochrone ($family, $photsys) interpolated for age=$age Gyr and MH=$metal using approximation $key.")
    else
        @assert 9.3≤log10(1e9*age)≤10.14 "For ezpadova interpolation age should satisfy 9.3≤log10(age_yr)≤10.14"
        @assert metal<0.41
        file_artif = joinpath(dirpath, "webapi","logAge_9.3to10.14_metal_-2.19to0.41_Niso_390.dat")
        quickiso =  ezpadova.QuickInterpolator(file_artif)
        df = quickiso(log(age*1e9), metal) |> PyPandasDataFrame |> DataFrame
        df.label .=  Int.(floor.(df.evol)) # Recompute label so as not to have inerpolated value
        println("✅ Isochrone ($family, $photsys) interpolated for age=$age Gyr and MH=$metal using
        ezpadova's QuickInterpolator.")
    end
    return df
end




