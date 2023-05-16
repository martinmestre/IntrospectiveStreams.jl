"""Plot histogram on sky."""
function plot_histog_on_sky(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=200)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_sky(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=200)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot histogram on sky with globular clusters."""
function plot_histog_on_sky_with_gc(df::DataFrame, df_gc::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=200)+data(df_gc)*visual(color="red"))*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")
    plt_M68 = data(df_gc)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")*visual(color="black")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_sky_with_gc(df::DataFrame, df_gc::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=200)+data(df_gc)*visual(color="red"))*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")
    plt_M68 = data(df_gc)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")*visual(color="black")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end


"""Plot histogram on sky using the stream's frame."""
function plot_histog_on_sky_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=100)+data(df_track))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_sky_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=100)+data(df_track))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_sky_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=100)*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_sky_self_frame(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=100)*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end
"""Scatter plot on sky using the stream's frame with and without track."""
function plot_scatter_on_sky_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=3, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_sky_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=3, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_sky_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=3, color=(:black,1))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_sky_self_frame(df::DataFrame, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=3, color=(:black,1))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Scatter plot on sky with PM (μ) arrows using stream's frame."""
function plot_scatter_on_sky_μ_arrows_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, step::Integer, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt,  axis=(;limits=window))
    us = df.μ₁cosϕ₂./cos.(df.ϕ₂*π/180.)
    vs = df.μ₂
    strength = @. sqrt(us^2+vs^2)
    @. us = us/strength
    @. vs = vs/strength
    arrows!(df.ϕ₁, df.ϕ₂, us, vs, arrowsize = strength, lengthscale = 1)
    us = (df_track.μ₁cosϕ₂./cos.(df_track.ϕ₂*π/180.))[begin:step:end]
    vs = df_track.μ₂[begin:step:end]
    strength = @. sqrt(us^2 + vs^2)
    @. us = us/strength
    @. vs = vs/strength
    arrows!(df_track.ϕ₁[begin:step:end], df_track.ϕ₂[begin:step:end], us, vs, arrowsize = strength, lengthscale = 1, color="blue")
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_sky_μ_arrows_self_frame(df::DataFrame, df_track::DataFrame, step::Integer, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)
    us = df.μ₁cosϕ₂./cos.(df.ϕ₂*π/180.)
    vs = df.μ₂
    strength = @. sqrt(us^2+vs^2)
    @. us = us/strength
    @. vs = vs/strength
    arrows!(df.ϕ₁, df.ϕ₂, us, vs, arrowsize = strength, lengthscale = 1)
    us = (df_track.μ₁cosϕ₂./cos.(df_track.ϕ₂*π/180.))[begin:step:end]
    vs = df_track.μ₂[begin:step:end]
    strength = @. sqrt(us^2 + vs^2)
    @. us = us/strength
    @. vs = vs/strength
    arrows!(df_track.ϕ₁[begin:step:end], df_track.ϕ₂[begin:step:end], us, vs, arrowsize = strength, lengthscale = 1, color="blue")
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Scatter plot on sky with reflex-corrected PM (μ) arrows using stream's frame."""
function plot_scatter_on_sky_μ_corr_arrows_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, step::Integer, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt,  axis=(;limits=window))
    us = df.μ₁_corr
    vs = df.μ₂_corr
    strength = @. sqrt(us^2+vs^2)
    @. us = us/strength
    @. vs = vs/strength
    arrows!(df.ϕ₁, df.ϕ₂, us, vs, arrowsize = strength, lengthscale = 1)
    us = df_track.μ₁_corr[begin:step:end]
    vs = df_track.μ₂_corr[begin:step:end]
    strength = @. sqrt(us^2 + vs^2)
    @. us = us/strength
    @. vs = vs/strength
    arrows!(df_track.ϕ₁[begin:step:end], df_track.ϕ₂[begin:step:end], us, vs, arrowsize = strength, lengthscale = 1, color="blue")
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_sky_μ_corr_arrows_self_frame(df::DataFrame, df_track::DataFrame, step::Integer, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)
    us = df.μ₁_corr
    vs = df.μ₂_corr
    strength = @. sqrt(us^2+vs^2)
    @. us = us/strength
    @. vs = vs/strength
    arrows!(df.ϕ₁, df.ϕ₂, us, vs, arrowsize = strength, lengthscale = 1)
    us = df_track.μ₁_corr[begin:step:end]
    vs = df_track.μ₂_corr[begin:step:end]
    strength = @. sqrt(us^2 + vs^2)
    @. us = us/strength
    @. vs = vs/strength
    arrows!(df_track.ϕ₁[begin:step:end], df_track.ϕ₂[begin:step:end], us, vs, arrowsize = strength, lengthscale = 1, color="blue")
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot histogram on μ plane."""
function plot_histog_on_μ_plane(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}},  file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=500)*mapping(:pmra =>L"$μ_{RA}$ [mas/yr]", :pmdec=>L"$μ_{Dec}$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_μ_plane(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=500)*mapping(:pmra =>L"$μ_{RA}$ [mas/yr]", :pmdec=>L"$μ_{Dec}$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot histogram on μ plane using stream's frame with and without track."""
function plot_histog_on_μ_plane_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=500)*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")

    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_μ_plane_self_frame(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=500)*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_μ_plane_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=200)+data(df_track)*visual(markersize=1))*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_on_μ_plane_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=200)+data(df_track)*visual(markersize=1))*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Scatter plot on μ plane in stream's self-frame with and without track."""
function plot_scatter_on_μ_plane_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1.5, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_μ_plane_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_μ_plane_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=1, color=(:black,1))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_μ_plane_self_frame(df::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=1, color=(:black,1))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Scatter plot on reflex-corrected μ plane in stream's self-frame with and without track."""
function plot_scatter_on_μ_corr_plane_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_μ_corr_plane_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_μ_corr_plane_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=1, color=(:black,1))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_scatter_on_μ_corr_plane_self_frame(df::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=1, color=(:black,1))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Track plot on reflex-corrected μ plane in stream's self-frame."""
function plot_track_on_μ_corr_plane_self_frame(df_track::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df_track)*visual(markersize=2,color="red")*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot CMD histogram with or without isochrone."""
function plot_histog_cmd(df::DataFrame, df_iso::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*mapping(:color=>L"$BP-RP$ [Mag]", :g_abs=>L"$G$ [Mag]")*histogram(bins=100)
    plt_iso = data(df_iso)*mapping(:color=>L"$BP-RP$ [Mag]", :Gaia_G_EDR3=>L"$G$ [Mag]")*visual(Lines,color="red")
    plt_iso_bord = data(df_iso)*mapping([:left,:right].=>L"$BP-RP$ [Mag]", :Gaia_G_EDR3=>L"$G$ [Mag]")*visual(Lines,color="black")
    ag = draw!(fig, plt+plt_iso+plt_iso_bord, axis=(;yreversed=true))#, limits=((0,1.5),(14,22))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_histog_cmd(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=300)*mapping(:color=>L"$BP-RP$ [Mag]", :g_abs=>L"$G$ [Mag]")
    ag = draw!(fig, plt, axis=(;yreversed=true))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot single isochrone."""
function plot_isochrone_cmd(df::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*mapping(:color=>L"BP-RP", :Gaia_G_EDR3 =>L"G")*visual(Lines)
    draw!(fig, plt, axis=(;yreversed=true))
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end



