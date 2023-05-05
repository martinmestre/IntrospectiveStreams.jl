"""Extinction correction for Gaia magnitudes from Gaia dataset."""
function extinction_corr(file_orig::String, file_corr::String)::Nothing
    println("file = ", file_corr)
    data = at.Table.read(file_orig)
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

"""Data curation of Gaia dataset."""
function curation!(df::DataFrame, tol::Vector{Number})::Nothing
    df.pmra_error_rel = df.pmra_error./df.pmra
    @subset!(df, (:pmra_error_rel .< tol[1]))
    df.pmdec_error_rel = df.pmdec_error./df.pmdec
    @subset!(df, (:pmdec_error_rel .< tol[2]))
    df.parallax_rel_error = df.parallax_error./df.parallax
    @subset!(df, (:parallax_rel_error .< tol[3]))
    return nothing
end

"""CMD filtering (Mutating)."""
function filter_cmd!(df_stream::DataFrame, df_iso::DataFrame)::Nothing
    phase_mask = 0 .<= df_iso.phase .< 3
    df_iso = df_iso[phase_mask,:]
    df_iso.color = df_iso.Gaia_BP_EDR3 - df_iso.Gaia_RP_EDR3
    df_iso.left = df_iso.color .- 0.1
    df_iso.right = df_iso.color .+ 0.1
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
function get_isochrone(age::Float64, metal::Float64,
                           phot::String, age_scale::String)::DataFrame
    println("Note that MIST uses [FeH] (not Z abundance).")
    df = ezmist.get_one_isochrone(age=age, FeH=metal, v_div_vcrit=0.0,
                    age_scale=age_scale, output_option="photometry",
                    output=phot, Av_value=0.0).to_pandas()|> PyPandasDataFrame |> DataFrame
    return df
end

"""Compute stream stars' self coordinates using Galstreams'track and add to dataframe."""
function compute_in_selfCoords!(df::DataFrame, frame::Py)::Nothing
    sky_coords = coord.SkyCoord(ra=Py(df.ra)*u.deg, dec=Py(df.dec)*u.deg, pm_ra_cosdec=Py(df.pmra)*u.mas/u.yr, pm_dec=Py(df.pmdec)*u.mas/u.yr, frame="icrs")
    self_coords = sky_coords.transform_to(frame)
    df.ϕ₁ = pyconvert(Vector{Float64}, self_coords.phi1.deg)
    df.ϕ₂ = pyconvert(Vector{Float64}, self_coords.phi2.deg)
    df.μ₁cosϕ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi1_cosphi2.value)
    @. df.μ₁ = df.μ₁cosϕ₂/cos(df.ϕ₂*π/180.0)
    df.μ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi2.value)
    @. df.D = 1.0/df.parallax
    return nothing
end

"""Compute stream stars' self coordinates using gala's KOPOSOV frame and add to dataframe."""
function compute_in_selfCoords!(df::DataFrame)::Py
    sky_coords = coord.SkyCoord(ra=Py(df.ra)*u.deg, dec=Py(df.dec)*u.deg, pm_ra_cosdec=Py(df.pmra)*u.mas/u.yr, pm_dec=Py(df.pmdec)*u.mas/u.yr, frame="icrs")
    kop_frame = galacoord.GD1
    self_coords = sky_coords.transform_to(kop_frame)
    df.ϕ₁ = pyconvert(Vector{Float64}, self_coords.phi1.deg)
    df.ϕ₂ = pyconvert(Vector{Float64}, self_coords.phi2.deg)
    df.μ₁cosϕ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi1_cosphi2.value)
    @. df.μ₁ = df.μ₁cosϕ₂/cos(df.ϕ₂*π/180.0)
    df.μ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi2.value)
    @. df.D = 1.0/df.parallax
    return kop_frame
end


"Mask out field globular clusters."
function mask_gc!(df_stream, df_gc)
    for i in 1:nrow(df_gc)
        Δra = df_stream.ra.-df_gc.ra[i]
        Δdec = df_stream.dec.-df_gc.dec[i]
        bool_gc .= sqrt.(Δra.^2 .+ Δdec.^2) .> 0.5
        @subset!(df_stream, identity(bool_gc))
    end
end

"""Filter with stream track on the sky."""
function filter_stream_on_sky!(df_stars::DataFrame, df_track::DataFrame, width::Number)
    up = df_track.ϕ₂.+width
    down =  df_track.ϕ₂.-width
    poly_ϕ₁ = vcat(df_track.ϕ₁, reverse(df_track.ϕ₁), df_track.ϕ₁[1])
    poly_ϕ₂ = vcat(down, reverse(up), down[1])
    polygon = SVector.(poly_ϕ₁, poly_ϕ₂)
    points = [[df_stars.ϕ₁[i], df_stars.ϕ₂[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stars, collect(inside))
    return nothing
end

"""Filter with stream on μ-space."""
function filter_stream_μ_space!(df_stars::DataFrame, df_track::DataFrame, Δμ::Number)
    left = df_track.μ₁cosϕ₂.-Δμ
    right =  df_track.μ₁cosϕ₂.+Δμ
    poly_y = vcat(df_track.μ₂, reverse(df_track.μ₂), df_track.μ₂[1])
    poly_x = vcat(left, reverse(right), left[1])
    polygon = SVector.(poly_x, poly_y)
    points = [[df_stars.μ₁cosϕ₂[i], df_stars.μ₂[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stars, collect(inside))
    return nothing
end

"""Non-mutating filter with stream track in any of its dimensions."""
function filter_with_track(df_stars::DataFrame, df_track::DataFrame, S::Symbol, σ::Number)::DataFrame
    if S == :ϕ₂
        q🌠 = df_stars.ϕ₂
        q_track = df_track.ϕ₂
    elseif S == :D
        q🌠 = 1.0 ./ df_stars.parallax
        q_track = df_track.D
    elseif S == :μ₁cosϕ₂_corr
        q🌠 = df_stars.μ₁cosϕ₂_corr
        q_track = df_track.μ₁cosϕ₂_corr
    elseif S == :μ₁_corr
        q🌠 = df_stars.μ₁_corr
        q_track = df_track.μ₁_corr
    elseif S == :μ₂_corr
        q🌠 = df_stars.μ₂_corr
        q_track = df_track.μ₂_corr
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
    return @subset(df_stars, collect(inside))
end

"""Mutating filter with stream track in any of its dimensions."""
function filter_with_track!(df_stars::DataFrame, df_track::DataFrame, S::Symbol, σ::Number)::Nothing
    if S == :ϕ₂
        q🌠 = df_stars.ϕ₂
        q_track = df_track.ϕ₂
    elseif S == :D
        q🌠 = 1.0 ./ df_stars.parallax
        q_track = df_track.D
    elseif S == :μ₁cosϕ₂_corr
        q🌠 = df_stars.μ₁cosϕ₂_corr
        q_track = df_track.μ₁cosϕ₂_corr
    elseif S == :μ₁_corr
        q🌠 = df_stars.μ₁_corr
        q_track = df_track.μ₁_corr
    elseif S == :μ₂_corr
        q🌠 = df_stars.μ₂_corr
        q_track = df_track.μ₂_corr
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
    @subset!(df_stars, collect(inside))
    return nothing
end

"""Filter with constant value."""
function filter_with_ϕ₂!(df::DataFrame, σ::Number)::Nothing
    @subset!(df, abs.(:ϕ₂) .< σ )
    return nothing
end

"""Filter as in PWB18."""
function filter_PWB18!(df::DataFrame)::Nothing
    @subset!(df, abs.(:ϕ₂.-0.5) .< 0.75)
    @subset!(df, -50 .< :ϕ₁ .< -10)
    return nothing
end

"""Filter with a box in PM space."""
function filter_box_μ(df::DataFrame, box::Vector{Vector{Float64}})::Nothing
    @subset!(df, box[1][1] .< :μ₁_corr .< box[1][2])
    @subset!(df_box, box[2][1] .< :μ₂_corr .< box[2][2])
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
    df.μ₁_corr = df.μ₁cosϕ₂_corr ./ cos.(df.ϕ₂*π/180.)
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


"""Pipeline for a single stream."""
function pipeline(name_s, name_t)

    # """Performing extinction correction. Only once for each stream."""
    # file_orig = "data/$(name_s)/DR3/original/GaiaDR3-$(name_t)-all.fits"
    # file_corr = "data/$(name_s)/DR3/corrected/corr_GaiaDR3-$(name_t)-all.fits"
    # dm.extinction_correction(file_orig, file_corr)
    # GC.gc()

    """Opening the file with extinction corrected magnitudes."""
    file_corr = "data/$(name_s)/DR3/corrected/corr_GaiaDR3-$(name_t)-all.fits"
    f = FITS(file_corr)
    df_stream = DataFrame(f[2])
    println("Fields: ", names(df_stream))
    GC.gc()

    """Remove known globular clusters."""
    df_gc = DataFrame(CSV.File("data/gc_catalog/Baumgardt/orbits_table.txt", delim=" ", ignorerepeated=true))
    dm.rename!(df_gc,[:RA, :DEC] .=> [:ra, :dec])
    dm.mask_gc!(df_stream, df_gc)

    """Load galstreams data."""
    # name_t = "GD-1-I21"
    # name_t = "Fjorm-I21"
    name_t = "M68-P19"
    df_track, self_frame = dm.load_stream_track(name_t)
    D_interp = linear_interpolation(df_track.ϕ₁, df_track.D)  # only activate if needed.
    # %%

    """Compute the stream self-coordinates and reflex correction for both data and track."""

    # kop_frame = dm.compute_in_selfCoords!(df_stream) # A test for Gd-1
    dm.compute_in_selfCoords!(df_stream, self_frame)
    # D_linGD1(x) = 0.05*x+10.0
    @subset!(df_stream, minimum(df_track.ϕ₁) .< :ϕ₁ .< maximum(df_track.ϕ₁))
    # df_stream.D = D_linGD1.(df_stream.ϕ₁)
    df_stream.D = D_interp.(df_stream.ϕ₁)
    dm.reflex_correct!(df_stream, self_frame)
    dm.reflex_correct!(df_track, self_frame)
    # %%

    """CMD filtering."""
    iso_file = "data/products/iso_stream_$(name_s).csv"
    filters = "UBVRIplus"
    df_iso = dm.get_isochrone(11.2e9, -2.2, filters, "linear")
    CSV.write(iso_file, df_iso)
    # df_iso = DataFrame(CSV.File(iso_file))
    df_stream.color = df_stream.bp - df_stream.rp
    df_iso.color = df_iso.bp - df_iso.rp
    dm.filter_cmd!(df_stream, df_iso)
    # %%

    """Curation."""

    # dm.curation!(df_cmd)
    # df_stream.D .= 1.0 ./ df_stream.parallax # not good proxi
    # D_linGD1(x) = 0.05*x+10.0
    @subset!(df_stream, -30 .< :ϕ₁ .< 50)
    df_stream.D = D_linGD1.(df_stream.ϕ₁)
    df_stream.D = D_interp.(df_stream.ϕ₁)
    @subset!(df_stream, -9 .< :μ₁_corr .< -4.5)
    @subset!(df_stream, -1.7 .< :μ₂_corr .< 1)
    @subset!(df_stream, -0.75 .< :color .< 2.)
    @subset!(df_stream, :parallax .< 1.)
    # @subset!(df_stream, :D .> 0 )
    # %%

    """Apply different filters to the stream."""
    file_filt = "data/products/filt_GaiaDR3-$(name_t).fits"
    σ = 0.2

    S = :μ₁_corr
    df_filt = dm.filter_with_track(df_stream, df_track, S, σ)

    S₂ = :μ₂_corr
    dm.filter_with_track!(df_filt, df_track, S₂, σ)

    S₃ = :ϕ₂
    df_filt = dm.filter_with_track(df_stream, df_track, S₃, σ)
    @subset!(df_filt, 0 .< :ϕ₁ .< 70)

    df_filt  = dm.filter_with_ϕ₂(df_stream, σ)

    df_filt  = dm.filter_PWB18(df_stream)


    box = [[3.,12.],[-1,0]]
    box = [[6.5,9.],[-3.8,-2.5]]
    df_box = dm.filter_box_μ(df_stream, box)

    CSV.write(file_filt, df_filt)
    df_filt = CSV.File(file_filt) |> DataFrame
    # %%

    """Do some plots."""

    pm.plot_sky_scatter_selfFrame(df_box, "plots/dr2.pdf", df_track)
    pm.plot_sky_scatter_μ_arrows_selfFrame(df_filt[begin:1:end,:], "plots/sky_scatter_frame_μ_$(name_s)_filt.png", df_track)
    pm.plot_sky_scatter_μ_arrows_corr_selfFrame(df_filt[begin:1:end,:], "plots/sky_scatter_frame_μ_coor_$(name_s)_filt.png", df_track )

    window = [[0.,15.],[-5.,2.]]
    pm.plot_μ_scatter_selfFrame_window(df_filt, df_track, "plots/test.png",  window)
    pm.plot_μ_corr_scatter_selfFrame(df_filt, df_track, "plots/test.png")
    pm.plot_μ_corr_scatter_selfFrame_window(df_filt, df_track, "plots/μ_refCorr_selfFrame_GD-1_DR3.png",  window)
    pm.plot_μ_corr_histo_selfFrame_window(df_filt, df_track, "plots/μ_refCorr_selfFrame_GD-1.png",  window)
    pm.plot_μ_corr_track_selfFrame(df_track, "plots/test.png")


    pm.plot_isochrone_data(df_iso, df_box, "plots/test_cmd.png")
    pm.plot_cmd_histo(df_filt, "plots/test_cmd.png")
end