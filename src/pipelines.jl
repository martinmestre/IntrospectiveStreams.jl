"""Performing extinction correction for a given list of streams."""
function correct_extinction_Gaia_loop(name_s::Vector{String})
    for i in 1:length(name_s)
        println("Correcting extinction of stream $(name_s[i])")
        host_dir = "/home/mmestre/casa/work/data/cats"
        file_orig = "$(host_dir)/GaiaDR3-$(name_s[i])-all.fits"
        file_corr = "$(host_dir)/GaiaDR3-$(name_s[i])-all_extincorr.fits"
        correct_extinction_Gaia(file_orig, file_corr)
    end
    GC.gc()
end

"""Basic filtering with μ and cmd cuts."""
function basic_pipeline(name_s, name_t)


    """Opening the file with extinction corrected magnitudes."""
    host_dir = "/home/mmestre/casa/work/data/cats"
    file_corr = "$(host_dir)/GaiaDR3-$(name_s)-all_extincorr.fits"
    f = FITS(file_corr)
    df_stream = DataFrame(f[2])
    println("Fields: ", names(df_stream))
    GC.gc()

    """Remove known globular clusters."""
    df_gc = DataFrame(CSV.File("data/gc_catalog/Baumgardt/orbits_table.txt", delim=" ", ignorerepeated=true))
    dm.rename!(df_gc,[:RA, :DEC] .=> [:ra, :dec])
    dm.mask_gc!(df_stream, df_gc)

    """Load galstreams data."""
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
