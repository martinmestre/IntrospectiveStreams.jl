"""Load dataframes from files."""
function get_dataframes(file_corr, file_phot, file_gc, file_iso, age, metal, filter)::Vector{Union{DataFrame,Nothing}}
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
        df_iso = get_isochrone(age, metal, filter, "linear")
        CSV.write(file_iso, df_iso)
    end
    array_df = [df_astrom, df_phot, df_gc, df_iso]
    return array_df
end

"""Extinction correction for Gaia magnitudes from Gaia dataset."""
function correct_extinction_Gaia(file_orig::String, file_corr::String)::Nothing
    println(file_orig, file_corr)
    data = at.Table.read(file_orig)
    println("done")
    g = pyia.GaiaData(data)
    bp0 = g.get_BP0()
    rp0 = g.get_RP0()
    g0 = g.get_G0()
    data["bp"] = bp0
    data["rp"] = rp0
    data["g"] = g0
    data.write(file_corr, format="fits", overwrite=true)
    return nothing
end

"""Performing extinction correction for a given list of streams."""
function correct_extinction_Gaia_loop(name_s::Vector{String}, file_orig::Vector{String}, file_corr::Vector{String})
    for i in 1:length(name_s)
        println("Correcting extinction of stream $(name_s[i])")
        correct_extinction_Gaia(file_orig[i], file_corr[i])
    end
    GC.gc()
end

"""Data curation of Gaia dataset."""
function curation!(df::DataFrame, tol::Vector{<:Number})::Nothing
    df.pmra_error_rel = df.pmra_error./df.pmra
    @subset!(df, (:pmra_error_rel .< tol[1]))
    df.pmdec_error_rel = df.pmdec_error./df.pmdec
    @subset!(df, (:pmdec_error_rel .< tol[2]))
    df.parallax_rel_error = df.parallax_error./df.parallax
    @subset!(df, (:parallax_rel_error .< tol[3]))
    return nothing
end

"""CMD filtering (Mutating)."""
function filter_cmd!(df_stream::DataFrame, df_iso::DataFrame, σ_c::Number)::Nothing
    phase_mask = 0 .<= df_iso.phase .< 3
    df_iso = df_iso[phase_mask,:]
    df_iso.color = df_iso.Gaia_BP_EDR3 - df_iso.Gaia_RP_EDR3
    df_iso.left = df_iso.color .- σ_c
    df_iso.right = df_iso.color .+ σ_c
    pol_x = vcat(df_iso.left, reverse(df_iso.right), df_iso.left[1])
    temp = df_iso.Gaia_G_EDR3
    pol_y = vcat(temp, reverse(temp), temp[1])
    polygon = SVector.(pol_x, pol_y)

    df_stream.distmod = pyconvert(Vector{Float64},coord.Distance(Py(df_stream.D)*u.kpc).distmod.value)
    df_stream.g_abs = df_stream.g - df_stream.distmod
    points = [[df_stream.color[i], df_stream.g_abs[i]] for i in 1:nrow(df_stream) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stream, identity(inside))
    return nothing
end



"""Load Galstreams data (one single stream track)."""
function load_stream_track(name_t::String)
    mwsts = galstreams.MWStreams(verbose=false, implement_Off=true)
    resumen = mwsts.summary |> PyPandasDataFrame |> DataFrame
    bool_on = resumen.On .== true
    𝒯 = resumen.TrackName[bool_on]
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

"""Download stellar isochrones"""
function get_isochrone(age::Number, metal::Number,
                           phot::String, age_scale::String)::DataFrame
    println("Note that MIST uses [FeH] (not Z abundance).")
    df = ezmist.get_one_isochrone(age=age, FeH=metal, v_div_vcrit=0.0,
                    age_scale=age_scale, output_option="photometry",
                    output=phot, Av_value=0.0).to_pandas()|> PyPandasDataFrame |> DataFrame
    return df
end

"""Compute stream stars' self coordinates using Galstreams'track and add to dataframe."""
function compute_in_self_coords!(df::DataFrame, frame::Py)::Nothing
    sky_coords = coord.SkyCoord(ra=Py(df.ra)*u.deg, dec=Py(df.dec)*u.deg, pm_ra_cosdec=Py(df.pmra)*u.mas/u.yr, pm_dec=Py(df.pmdec)*u.mas/u.yr, frame="icrs")
    self_coords = sky_coords.transform_to(frame)
    df.ϕ₁ = pyconvert(Vector{Float64}, self_coords.phi1.deg)
    df.ϕ₂ = pyconvert(Vector{Float64}, self_coords.phi2.deg)
    df.μ₁cosϕ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi1_cosphi2.value)
    df.μ₁ = @. df.μ₁cosϕ₂/cos(df.ϕ₂*π/180.0)
    df.μ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi2.value)
    df.D_Π = 1.0./df.parallax
    return nothing
end

"""Compute stream stars' self coordinates using gala's KOPOSOV frame and add to dataframe."""
function compute_in_self_coords!(df::DataFrame)::Py
    sky_coords = coord.SkyCoord(ra=Py(df.ra)*u.deg, dec=Py(df.dec)*u.deg, pm_ra_cosdec=Py(df.pmra)*u.mas/u.yr, pm_dec=Py(df.pmdec)*u.mas/u.yr, frame="icrs")
    kop_frame = galacoord.GD1
    self_coords = sky_coords.transform_to(kop_frame)
    df.ϕ₁ = pyconvert(Vector{Float64}, self_coords.phi1.deg)
    df.ϕ₂ = pyconvert(Vector{Float64}, self_coords.phi2.deg)
    df.μ₁cosϕ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi1_cosphi2.value)
    df.μ₁ = @. df.μ₁cosϕ₂/cos(df.ϕ₂*π/180.0)
    df.μ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi2.value)
    df.D_Π = 1.0./df.parallax
    return kop_frame
end


"Mask out field globular clusters."
function mask_gc!(df_stream::DataFrame, df_gc::DataFrame)
    for i in 1:nrow(df_gc)
        Δra = df_stream.ra.-df_gc.ra[i]
        Δdec = df_stream.dec.-df_gc.dec[i]
        bool_gc = @. sqrt(Δra^2+Δdec^2) > 0.5
        @subset!(df_stream, identity(bool_gc))
    end
end

"""Filter with stream track on the sky."""
function filter_on_sky!(df_stars::DataFrame, df_track::DataFrame, width::Number)
    up = df_track.ϕ₂.+width
    down =  df_track.ϕ₂.-width
    poly_ϕ₁ = vcat(df_track.ϕ₁, reverse(df_track.ϕ₁), df_track.ϕ₁[1])
    poly_ϕ₂ = vcat(down, reverse(up), down[1])
    polygon = SVector.(poly_ϕ₁, poly_ϕ₂)
    points = [[df_stars.ϕ₁[i], df_stars.ϕ₂[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stars, identity(inside))
    return nothing
end

"""Filter with stream on μ-space."""
function filter_on_μ_plane!(df_stars::DataFrame, df_track::DataFrame, Δμ::Number)
    left = df_track.μ₁cosϕ₂.-Δμ
    right =  df_track.μ₁cosϕ₂.+Δμ
    poly_y = vcat(df_track.μ₂, reverse(df_track.μ₂), df_track.μ₂[1])
    poly_x = vcat(left, reverse(right), left[1])
    polygon = SVector.(poly_x, poly_y)
    points = [[df_stars.μ₁cosϕ₂[i], df_stars.μ₂[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stars, identity(inside))
    return nothing
end

"""Non-mutating filter with stream track in any of its dimensions."""
function filter_along_ϕ₁(df_stars::DataFrame, df_track::DataFrame, S::Symbol, σ::Number)::DataFrame
    if S == :ϕ₂
        q🌠 = df_stars.ϕ₂
        q_track = df_track.ϕ₂
    elseif S == :D
        q🌠 = df_stars.D
        q_track = df_track.D
    elseif S == :μ₁cosϕ₂
        q🌠 = df_stars.μ₁cosϕ₂
        q_track = df_track.μ₁cosϕ₂
    elseif S == :μ₁
        q🌠 = df_stars.μ₁
        q_track = df_track.μ₁
    elseif S == :μ₂
        q🌠 = df_stars.μ₂
        q_track = df_track.μ₂
    elseif S == :Vᵣ
        q🌠 = df_stars.radial_velocity
        q_track = df_track.Vᵣ
    end
    up = q_track .+ σ
    down =  q_track .- σ
    poly_ϕ₁ = vcat(df_track.ϕ₁, reverse(df_track.ϕ₁), df_track.ϕ₁[1])
    poly_q = vcat(down, reverse(up), down[1])
    polygon = SVector.(poly_ϕ₁, poly_q)
    points = [[df_stars.ϕ₁[i], q🌠[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    return @subset(df_stars, identity(inside))
end

"""Mutating filter with stream track in any of its dimensions."""
function filter_along_ϕ₁!(df_stars::DataFrame, df_track::DataFrame, S::Symbol, σ::Number)::Nothing
    if S == :ϕ₂
        q🌠 = df_stars.ϕ₂
        q_track = df_track.ϕ₂
    elseif S == :D
        q🌠 = df_stars.D
        q_track = df_track.D
    elseif S == :μ₁cosϕ₂
        q🌠 = df_stars.μ₁cosϕ₂
        q_track = df_track.μ₁cosϕ₂
    elseif S == :μ₁
        q🌠 = df_stars.μ₁
        q_track = df_track.μ₁
    elseif S == :μ₂
        q🌠 = df_stars.μ₂
        q_track = df_track.μ₂
    elseif S == :Vᵣ
        q🌠 = df_stars.radial_velocity
        q_track = df_track.Vᵣ
    end
    println("q_track=$(q_track)")
    up = q_track .+ σ
    down =  q_track .- σ
    poly_ϕ₁ = vcat(df_track.ϕ₁, reverse(df_track.ϕ₁), df_track.ϕ₁[1])
    poly_q = vcat(down, reverse(up), down[1])
    polygon = SVector.(poly_ϕ₁, poly_q)
    points = [[df_stars.ϕ₁[i], q🌠[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stars, identity(inside))
    return nothing
end

"""Filter with constant value."""
function filter_with_ϕ₂!(df::DataFrame, σ::Number)::Nothing
    @subset!(df, abs.(:ϕ₂) .< σ )
    return nothing
end

"""Filter with a box in PM space."""
function filter_box_on_μ_plane!(df::DataFrame, box::Vector{Vector{Float64}})::Nothing
    @subset!(df, box[1][1] .< :μ₁ .< box[1][2])
    @subset!(df, box[2][1] .< :μ₂ .< box[2][2])
    return nothing
end

"""Filter as in PWB18."""
function filter_from_PWB18!(df::DataFrame)::Nothing
    @subset!(df, abs.(:ϕ₂.-0.5) .< 0.75)
    @subset!(df, -50 .< :ϕ₁ .< -10)
    return nothing
end

"""Reflex Correction."""
function reflex_correct!(df::DataFrame, frame::Py)::Nothing
    len = length(df.ϕ₁)
    sky_coords = coord.SkyCoord(phi1=Py(df.ϕ₁)*u.deg, phi2=Py(df.ϕ₂)*u.deg, pm_phi1_cosphi2=Py(df.μ₁cosϕ₂)*u.mas/u.yr, pm_phi2=Py(df.μ₂)*u.mas/u.yr, distance=Py(df.D)*u.kpc, radial_velocity=Py(fill(0.,len))*u.km/u.s, frame=frame)
    vsun = coord.CartesianDifferential(Py([11.1, 220.0+12.24, 7.25])*u.km/u.s)
    rsun = 8.122*u.kpc
    gc_frame = coord.Galactocentric(galcen_distance=rsun, galcen_v_sun=vsun, z_sun=0*u.pc)
    sky_coords_corr = galacoord.reflex_correct(sky_coords, gc_frame)
    df.μ₁cosϕ₂_corr = pyconvert(Vector{Float64}, sky_coords_corr.pm_phi1_cosphi2.value)
    df.μ₁_corr = @. df.μ₁cosϕ₂_corr/cos(df.ϕ₂*π/180.0)
    df.μ₂_corr = pyconvert(Vector{Float64}, sky_coords_corr.pm_phi2.value)
    return nothing
end

"""Clean the cross-match."""
function clean_xmatch!(df::DataFrame)
    # g_0 and g are SPLUS and Gaia mags corrected by extinction.
    unicos = size(unique(df.ID))[1]
    n_multiple = m_multiple = 0
    for i=size(df)[1]:-1:1
        n_repe = count((df.ID.==df.ID[i]))
        if n_repe>1
            n_multiple += 1
            m_multiple += n_repe
            idx = findall(x->x==df.ID[i], df.ID)
            dif_mag = abs.([ df.g_0[idx[1]]-df.g[idx[1]], df.g_0[idx[2]]-df.g[idx[2]] ])
            if dif_mag[1] < dif_mag[2]
                deleteat!(df, [idx[2]])
            else
                deleteat!(df, [idx[1]])
            end
        end
    end
end
