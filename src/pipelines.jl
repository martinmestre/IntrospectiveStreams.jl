"""Basic filtering with μ and cmd cuts."""
function basic_pipeline(array_df, file_filt, name_t, tol_curation, col_bounds, σ_c, σ)
    df_astrom = array_df[1]
    df_phot = array_df[2]
    df_gc   = array_df[3]
    df_iso  = array_df[4]

    # Curation.
    curation!(df_astrom, tol_curation)

    # Remove known globular clusters.
    rename!(df_gc,[:RA, :DEC] .=> [:ra, :dec])
    mask_gc!(df_astrom, df_gc)

    # Load galstreams data.
    df_track, self_frame = load_stream_track(name_t)
    D_interp = linear_interpolation(df_track.ϕ₁, df_track.D)  # only activate if needed.

    # Compute the stream self-coordinates. Do not correct for the solar reflex motion.
    compute_in_self_coords!(df_astrom, self_frame)
    @subset!(df_astrom, minimum(df_track.ϕ₁) .< :ϕ₁ .< maximum(df_track.ϕ₁))
    df_astrom.D_track = D_interp.(df_stream.ϕ₁)

    # CMD filtering.
    df_astrom.Gaia_color = df_stream.bp - df_stream.rp
    @subset!(df_astrom, col_bounds[1] .< :Gaia_color .< col_bounds[2])
    filter_cmd!(df_astrom, df_iso, σ_c)

    # Spatial filtering.
    @subset!(df_astrom, :parallax .< 1.)
    df_filt = filter_with_track(df_stream, df_track, :μ₁, σ)
    filter_with_track!(df_filt, df_track, :μ₂, σ)

     # Saving filtered stream dataset.
    CSV.write(file_filt, df_filt)


    """Do some plots."""
    pm.plot_sky_scatter_selfFrame(df_box, "plots/scatter.pdf", df_track)
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


"""Processing all the pipelines in sequence."""
function basic_pipeline_loop(name_t, file_corr, file_phot, file_gc, file_filt, ages, metals,
    filters, tol_curation, col_bounds, σ_c, σ)
    for i in eachindex(name_t)
        array_df = get_dataframes(file_corr[i], file_phot[i], file_gc[i], file_isochrone[i], ages[i],
                   metals[i], filters[i])
        basic_pipeline(array_df, file_filt[i], name_t[i], tol_curation, col_bounds, σ_c, σ)
    end
    GC.gc()
    return nothing
end