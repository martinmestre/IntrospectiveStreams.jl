"""Basic filtering with μ and cmd cuts."""
function basic_pipeline(array_df, file_filt, file_plot, name_t, tol_curation, col_bounds, σ_c, σ)
    df_astrom = array_df[1]
    df_phot = array_df[2]
    df_gc   = array_df[3]
    df_iso  = array_df[4]

    # Curation.
    # curation!(df_astrom, tol_curation)

    # Remove known globular clusters.
    rename!(df_gc,[:RA, :DEC] .=> [:ra, :dec])
    mask_gc!(df_astrom, df_gc)

    # Load galstreams data.
    df_track, self_frame = load_stream_track(name_t)
    D_interp = linear_interpolation(df_track.ϕ₁, df_track.D)  # only activate if needed.

    # Compute the stream self-coordinates. Do not correct for the solar reflex motion.
    compute_in_self_coords!(df_astrom, self_frame)
    @subset!(df_astrom, minimum(df_track.ϕ₁) .< :ϕ₁ .< maximum(df_track.ϕ₁))
    df_astrom.D = D_interp.(df_astrom.ϕ₁)

    # CMD filtering.
    df_astrom.color = df_astrom.bp - df_astrom.rp
    @subset!(df_astrom, col_bounds[1] .< :color .< col_bounds[2])
    filter_cmd!(df_astrom, df_iso, σ_c)

    # Spatial filtering.
    @subset!(df_astrom, :parallax .< 1.)
    # filter_with_ϕ₂!(df_astrom, 1.0)
    # df_filt = filter_along_ϕ₁(df_astrom, df_track, :μ₁, σ)
    # filter_along_ϕ₁!(df_filt, df_track, :μ₂, σ)

     # Saving filtered stream dataset.
    # CSV.write(file_filt, df_filt)

    """Do some plots."""
    # plot_scatter_on_sky_self_frame(df_astrom, df_track, file_plot)
    window = ((-15,-10),(-5,0))
    box = [[-14,-10.],[-4.,-2.]]
    filter_box_on_μ_plane!(df_astrom, box)
    plot_scatter_on_μ_plane_self_frame(df_astrom, df_track, window, file_plot)
    plot_scatter_on_sky_self_frame(df_astrom, df_track, file_plot)
end

"""Processing all the pipelines in sequence."""
function basic_pipeline_loop(name_t, file_corr, file_phot, file_gc, file_iso, file_filt, file_plot,
    age, metal, filter, tol_curation, col_bounds, σ_c, σ)
    for i in eachindex(name_t)
        array_df = get_dataframes(file_corr[i], file_phot[i], file_gc, file_iso[i], age[i], metal[i], filter[i])
        basic_pipeline(array_df, file_filt[i], file_plot[i], name_t[i], tol_curation, col_bounds, σ_c, σ)
    end
    GC.gc()
    return nothing
end