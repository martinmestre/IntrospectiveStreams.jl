"""Basic filtering with μ and cmd cuts."""
function basic_pipeline(name_s, array_df, file_filt, file_plot, name_t,
    tol_curation, col_bounds, box_μ, σ_c, σ)
    df_astrom = array_df[1]
    df_phot = array_df[2]
    df_gc   = array_df[3]
    df_iso  = array_df[4]

    println("Now in stream $(name_s)")
    # Curation.
    # curation!(df_astrom, tol_curation)

    # Remove known globular clusters.
    rename!(df_gc,[:RA, :DEC] .=> [:ra, :dec])
    mask_gc!(df_astrom, df_gc)

    # Load galstreams data.
    df_track, self_frame = load_stream_track(name_t)
    D_interp = linear_interpolation(df_track.ϕ₁, df_track.D)  # only activate if needed.

    # Compute the stream self-coordinates and correcte for reflex-motion of the ⊙.
    compute_in_self_coords!(df_astrom, self_frame)
    @subset!(df_astrom, minimum(df_track.ϕ₁) .< :ϕ₁ .< maximum(df_track.ϕ₁))
    df_astrom.D = D_interp.(df_astrom.ϕ₁)
    df_astrom2 = deepcopy(df_astrom)
    reflex_correct!(df_astrom, self_frame)
    reflex_correct_steps!(df_astrom2, self_frame)
    Δ₁ = df_astrom.μ₁cosϕ₂_rc - df_astrom2.μ₁cosϕ₂_rc
    Δ₂ = df_astrom.μ₂_rc - df_astrom2.μ₂_rc
    @show maximum.([Δ₁,Δ₂])

    # Raw phase-space filtering
    @subset!(df_astrom, :parallax .< 1.)
    # filter_with_ϕ₂!(df_astrom, 1.0)
    filter_along_ϕ₁!(df_astrom, df_track, :μ₁, σ)
    filter_along_ϕ₁!(df_astrom, df_track, :μ₂, σ)
    # filter_box_on_μ_plane!(df_astrom, box_μ)

    # CMD filtering.
    df_astrom.color = df_astrom.bp - df_astrom.rp
    @subset!(df_astrom, col_bounds[1] .< :color .< col_bounds[2])
    filter_cmd!(df_astrom, df_iso, σ_c)

    # Saving filtered stream dataset.
    CSV.write(file_filt, df_astrom)

    """Do some plots."""
    # window = ((-15,15),(-10,10))
    # plot_scatter_on_μ_plane_self_frame(df_astrom, df_track, file_plot)
    fig = plot_scatter_on_sky_self_frame(name_s, df_astrom, file_plot)
    electrondisplay(fig)
    GC.gc()
end

"""Processing all the pipelines in sequence."""
function basic_pipeline_loop(name_s, name_t, file_corr, file_phot, file_gc, file_iso, file_filt, file_plot,
    family, age, metal, filter, tol_curation, col_bounds, box_μ, σ_c, σ)
    for i in eachindex(name_t)
        array_df = get_dataframes(file_corr[i], file_phot[i], file_gc, file_iso[i], family[i], age[i], metal[i], filter[i])
        basic_pipeline(name_s[i], array_df, file_filt[i], file_plot[i], name_t[i], tol_curation, col_bounds,
        box_μ, σ_c, σ)
    end
    return nothing
end

