"""Load dataframes from files."""
function get_dataframes(file_corr, file_phot, file_gc, file_iso, family, age, metal, filter)::Vector{Union{DataFrame,Nothing}}
    f = FITS(file_corr)
    df_astrom = DataFrame(f[2])
    if isnothing(file_phot)
        df_phot = nothing
    else
        f = FITS(file_phot)
        df_phot = DataFrame(f[2])
    end
    df_gc  = DataFrame(CSV.File(file_gc, delim=" ", ignorerepeated=true))
    if isfile(file_iso)
        df_iso = DataFrame(CSV.File(file_iso))
    else
        df_iso = get_isochrone(family, age, metal, filter)
        CSV.write(file_iso, df_iso)
    end
    array_df = [df_astrom, df_phot, df_gc, df_iso]
    return array_df
end


"""Load Galstreams data (one single stream track)."""
function load_stream_track(name_t::String)
    mwsts = galstreams.MWStreams(verbose=false, implement_Off=true)
    # resumen = mwsts.summary |> PyPandasDataFrame |> DataFrame
    # bool_on = resumen.On .== true
    # 𝒯 = resumen.TrackName[bool_on]
    track = mwsts[name_t]
    frame = track.stream_frame
    self_coords = track.track.transform_to(frame)
    ϕ₁ = pyconvert(Vector{Float64}, self_coords.phi1.deg)
    ϕ₂ = pyconvert(Vector{Float64}, self_coords.phi2.deg)
    D = pyconvert(Vector{Float64}, self_coords.distance.value)
    μ₁cosϕ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi1_cosphi2.value)
    μ₁ = μ₁cosϕ₂ ./ cos.(ϕ₂*π/180.)
    μ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi2.value)
    Vᵣ = pyconvert(Vector{Float64}, self_coords.radial_velocity.value)
    df_track = DataFrame(ra=pyconvert(Vector{Float64},track.track.ra.deg),
                        dec=pyconvert(Vector{Float64},track.track.dec.deg),
                        ϕ₁=ϕ₁, ϕ₂=ϕ₂, μ₁cosϕ₂=μ₁cosϕ₂, μ₁=μ₁, μ₂=μ₂, D=D, Vᵣ=Vᵣ)
    return df_track, frame
end