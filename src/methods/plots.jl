
""" function plot_scatter_on_sky(df::DataFrame, coords::Tuple{})
Scatter plot on sky coordinates"""
function plot_scatter_on_sky(df::DataFrame,
                            photsys::Symbol,
                            coord::Tuple{Symbol,Symbol}=(:ra,:dec),
                            label::Tuple{String,String}=("RA~[°]","Dec~[°]"))
    filename = "sky_coords_$(photsys).pdf"
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
    label₁ = L"%$(label[1])"
    label₂ = L"%$(label[2])"
    plt = data(df) *
        visual(markersize=3, color=(:black,0.4)) *
        mapping(coord[1] => label₁, coord[2] => label₂)
    draw!(fig[1,1], plt)
    return fig, filename
end

function plot_scatter_on_sky(df_track::DataFrame, df::DataFrame,
                            photsys::Symbol,
                            coord::Tuple{Symbol,Symbol}=(:ra,:dec),
                            label::Tuple{String,String}=("RA~[°]","Dec~[°]"))
    filename = "sky_coords_$(photsys).pdf"
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
    label₁ = L"%$(label[1])"
    label₂ = L"%$(label[2])"

    plt = data(df) * visual(markersize=3, color=(:black,0.4)) *
            mapping(coord[1] => label₁, coord[2] => label₂)
    plt_track = data(df_track) * visual(Lines, color=wongcolors()[2]) *
            mapping(coord[1] => label₁, coord[2] => label₂)

    range₁ = (minimum(df[!,coord[1]]), maximum(df[!,coord[1]]))
    range₂ = (minimum(df[!,coord[2]]), maximum(df[!,coord[2]]))
    axis = (; limits = (range₁, range₂))
    draw!(fig[1,1], plt+plt_track, axis = axis)
    return fig, filename
end

function plot_scatter_on_sky(df_track::DataFrame, stream_name::Symbol, df::DataFrame,
    photsys::Symbol,
    coord::Tuple{Symbol,Symbol}=(:ra,:dec),
    label::Tuple{String,String}=("\\alpha~[°]","\\delta~[°]");
    ℚ::Symbol=:D, label_ℚ::String="D~[\\textrm{kpc}]")
    filename = "sky_coords_$(photsys)_$(stream_name)_track.pdf"
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30, figure_padding = (10,10,10,10))
    label₁ = L"%$(label[1])"
    label₂ = L"%$(label[2])"

    plt = data(df) * visual(markersize=3, color=(:black,0.4)) *
    mapping(coord[1] => label₁, coord[2] => label₂)
    plt_track = data(df_track) * visual(Scatter) *
    mapping(coord[1] => label₁, coord[2] => label₂, color=ℚ)

    range₁ = (minimum(df[!,coord[1]]), maximum(df[!,coord[1]]))
    range₂ = (minimum(df[!,coord[2]]), maximum(df[!,coord[2]]))
    axis = (; limits = (range₁, range₂), title=L"Orphan-Chenab~stream")
    grid = draw!(fig[1,1], plt+plt_track, scales(Color=(;colormap=:viridis)), axis = axis)
    colorbar!(fig[1,2], grid, label=L"%$(label_ℚ)")
    return fig, filename
end

"""Plot histogram on sky."""
function plot_histog_on_sky(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=200)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_sky(df::DataFrame, file::String="test.pdf")
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=200)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
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
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_sky_with_gc(df::DataFrame, df_gc::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=200)+data(df_gc)*visual(color="red"))*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")
    plt_M68 = data(df_gc)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")*visual(color="black")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end


"""Plot histogram on sky using the stream's frame."""
function plot_histog_on_sky_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=100)+data(df_track))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_sky_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=100)+data(df_track))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_sky_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=100)*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_sky_self_frame(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=100)*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end
"""Scatter plot on sky using the stream's frame with and without track."""
function plot_scatter_on_sky_self_frame(name_s::String, df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=3, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=window,title=name_s))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_scatter_on_sky_self_frame(name_s::String, df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=3, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;title=name_s))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end


"""
    plot_scatter_on_sky_self_frame(name_s::String, df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
"""
function plot_scatter_on_sky_self_frame(name_s::String, df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=3, color=(:black,1))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=window, title=name_s))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_scatter_on_sky_self_frame(name_s::String, df::DataFrame, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=3, color=(:black,1))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;title=name_s))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end


"""Scatter plot on sky with PM (μ) arrows using stream's frame."""
function plot_scatter_on_sky_μ_arrows_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, step::Integer, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
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
    save(file, fig, pt_per_unit=1)
    return fig
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
    save(file, fig, pt_per_unit=1)
    return fig
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
    save(file, fig, pt_per_unit=1)
    return fig
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
    save(file, fig, pt_per_unit=1)
    return fig
end

"""Plot histogram on μ plane."""
function plot_histog_on_μ_plane(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}},  file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=500)*mapping(:pmra =>L"$μ_{RA}$ [mas/yr]", :pmdec=>L"$μ_{Dec}$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_μ_plane(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=500)*mapping(:pmra =>L"$μ_{RA}$ [mas/yr]", :pmdec=>L"$μ_{Dec}$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

"""Plot histogram on μ plane using stream's frame with and without track."""
function plot_histog_on_μ_plane_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=500)*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_μ_plane_self_frame(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=500)*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_μ_plane_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=200)+data(df_track)*visual(markersize=1))*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_on_μ_plane_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=200)+data(df_track)*visual(markersize=1))*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

"""Scatter plot on μ plane in stream's self-frame with and without track."""
function plot_scatter_on_μ_plane_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=3, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_scatter_on_μ_plane_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=3, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_scatter_on_μ_plane_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=3, color=(:black,1))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_scatter_on_μ_plane_self_frame(df::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=3, color=(:black,1))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

"""Scatter plot on reflex-corrected μ plane in stream's self-frame with and without track."""
function plot_scatter_on_μ_corr_plane_self_frame(df::DataFrame, df_track::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_scatter_on_μ_corr_plane_self_frame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_scatter_on_μ_corr_plane_self_frame(df::DataFrame, window::Tuple{Tuple{Number,Number},Tuple{Number,Number}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=1, color=(:black,1))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=window))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_scatter_on_μ_corr_plane_self_frame(df::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*visual(markersize=1, color=(:black,1))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

"""Track plot on reflex-corrected μ plane in stream's self-frame."""
function plot_track_on_μ_corr_plane_self_frame(df_track::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df_track)*visual(markersize=2,color="red")*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
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
    save(file, fig, pt_per_unit=1)
    return fig
end

function plot_histog_cmd(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*histogram(bins=300)*mapping(:color=>L"$BP-RP$ [Mag]", :g_abs=>L"$G$ [Mag]")
    ag = draw!(fig, plt, axis=(;yreversed=true))
    colorbar!(fig[1,2], ag)
    save(file, fig, pt_per_unit=1)
    return fig
end

"""Plot stellar data and single isochrone cmd."""
function plot_cmd(df_stars::DataFrame, df::DataFrame, family::Symbol, photsys::Symbol, mag::Symbol, color::Symbol; only::Vector{T}=Int[]) where {T<:Integer}
    filename = "data_and isochrone_cmd_$(family)_$(photsys)_$(color).pdf"
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
    df_view = isempty(only) ? df : @view df[findall(row -> row.label in only, eachrow(df)), :]
    plt = data(df_view)*mapping(color=>"$(color)", mag =>"uppercase($(mag))", color=:phase=>"Phase")*
    visual(Lines, linewidth=2)
    plt_stars = data(df_stars)*mapping(color=>"$(color)", Symbol(mag,:_abs)=>"uppercase($(mag))")*
    visual(Scatter, markersize=2, color=:black)
    grid = draw!(fig, plt+plt_stars, scales(Color = (; palette = :Set1_9)), axis=(title="$(photsys) CMD", yreversed=true))
    legend!(fig[1,2],grid,; position=:right, titleposition=:top, framevisible=true, padding=5)
    return fig, filename
end

"""Plot single isochrone cmd."""
function plot_cmd(df::DataFrame, family::Symbol, photsys::Symbol, mag::Symbol, color::Symbol; only::Vector{T}=Int[]) where {T<:Integer}
    filename = "isochrone_cmd_$(family)_$(photsys)_$(color).pdf"
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)
    df_view = isempty(only) ? df : @view df[findall(row -> row.label in only, eachrow(df)), :]
    plt = data(df_view)*mapping(color=>"$(color)", mag =>"$(mag)", color=:phase=>"Phase")*
    visual(Lines,linewidth=2)
    grid = draw!(fig, plt, scales(Color = (; palette = :Set1_9)), axis=(title="$(photsys) CMD", yreversed=true))
    legend!(fig[1,2],grid,; position=:right, titleposition=:top, framevisible=true, padding=5)
    return fig, filename
end

function mathL(x::String)
    return L"%$(x)"
end

"""Plot multiple isochrones cmd."""
function plot_cmd(
    df_array::Vector{DataFrame},
    algorithm::Vector{Symbol},
    family::Symbol,
    photsys::Symbol,
    mag::Symbol,
    color::Symbol;
    only::Vector{T}=Int[], paramstring::String=""
) where {T<:Integer}

    filename = "isochrone_cmd_$(family)_$(photsys)_$(color)$(paramstring).pdf"
    size_inches = (12, 9)
    size_pt = 72 .* size_inches
    fig = Figure(size = size_pt, fontsize = 30)

    df_merged = vcat([transform(df[:, [mag, color, :phase, :label]],
    [] => (() -> algo) => :algorithm)
    for (df, algo) in zip(df_array, algorithm)]...)
    df_only = isempty(only) ? df_merged : @view df_merged[findall(row -> row.label in only, eachrow(df_merged)), :]


    label₁ = string(color)[end-1]*"-"*string(color)[end]*" [mag]"
    label₂ = string(mag)[1]*" [mag]"
    spec = data(df_only) *
            mapping(color=>label₁, mag=>label₂, color=:phase, linestyle=:algorithm=>presorted, linewidth=:algorithm=>presorted) *
            visual(Lines)

    # Dibujar sobre el eje existente y pasar paleta directamente
    sc = scales(; Color = (; palette = :Set1_9), LineWidth = (; palette = [1, 4, 1][1:length(algorithm)]), LineStyle=(; palette = [:solid, :dot, :dash][1:length(algorithm)]))
    grid = draw!(fig[1,1], spec, sc; axis=(;yreversed = true))

    # Leyenda en la columna de la derecha
    legend!(fig[1, 2], grid;
        position = :right,
        titleposition = :top,
        framevisible = true,
        padding = 5
    )

    return fig, filename
end


function plot_mags_density(df::DataFrame, mags::Vector{Symbol}, paleta,
                            photsys::Symbol; long=false, kwargs...)
    pattern =join(string.(mags), "")
    filename = "mags_$(photsys)_$(pattern)_density.pdf"
    size_inches = (11, 7)
    size_pt = 72 .* size_inches
    labels = latexstring.(mags)
    if long==true
        df_v = @view df[findall(row -> row.photfilter in mags, eachrow(df)), :]
        plt =   data(df_v) * aog.density(; kwargs...) *
                mapping(:magnitude => L"\mathrm{magnitude}") *
                mapping(color=:photfilter => presorted => L"\mathrm{filter}")
    else
        labels = latexstring.(mags)
        plt =   data(df) * aog.density(; kwargs...) *
                mapping(mags .=> L"\mathrm{magnitude}") *
                mapping(color=dims(1) => renamer(labels) => L"\mathrm{filter}")
    end
    sc = scales(; Color = (; palette = paleta ))
    axis = (; xgridvisible=false, ygridvisible=false, xticks=10:2:30, ylabel=L"\mathrm{PDF}")

    fig = Figure(size = size_pt, fontsize = 30)
    grid = draw!(fig[1,1], plt,sc, axis=axis)
    legend!(fig[1,1], grid; tellwidth=false, halign=:left, valign=:top, margin=(10, 10, 10, 10), patchsize=(20,20))
    return fig, filename
end

function plot_mags_histogram(df::DataFrame, mags::Vector{Symbol}, paleta,
                            photsys::Symbol; long::Bool=false, dodge::Bool=false, ylog::Bool=false, kwargs...)
    pattern =join(string.(mags), "")
    filename = "mags_$(photsys)_$(pattern)_histogram.pdf"
    size_inches = (11, 7)
    size_pt = 72 .* size_inches
    datalims = kwargs[1]

    if long==true
        df_v = @view df[findall(row -> row.photfilter in mags, eachrow(df)), :]
        plt =   data(df_v) * aog.histogram(; kwargs...) *
                mapping(:magnitude => L"\mathrm{magnitude}") *
                mapping(color=:photfilter => presorted => L"\mathrm{filter}") *
                visual(alpha=0.5)
        if dodge
            plt *= mapping(dodge=:photfilter => presorted)
        end
    else
        labels = latexstring.(mags)
        plt =   data(df) * aog.histogram(; kwargs...) *
                mapping(mags .=> L"\mathrm{magnitude}") *
                mapping(color=dims(1) => renamer(labels) => L"\mathrm{filter}") *
                visual(alpha=0.5)
        if dodge
            plt *= mapping(dodge=dims(1))
        end
    end

    sc = scales(; Color = (; palette = paleta ))
    if ylog==true
        exp_rng=range(0,6,step=1)
        ytickpos = (e->10.0.^e).(exp_rng)
        yticknames=replace.(Showoff.showoff(10. .^(exp_rng), :scientific),"1.0×"=> "" )
        axis = (; xgridvisible=false, ygridvisible=false, xticks=10:2:30, ylabel=L"\mathrm{count}", yscale=log, yticks = (ytickpos, yticknames), limits=(datalims,(1,nothing)))
    else
        axis = (; xgridvisible=false, ygridvisible=false, xticks=10:2:30, ylabel=L"\mathrm{count}")
    end

    fig = Figure(size = size_pt, fontsize = 30)
    grid = draw!(fig[1,1], plt,sc, axis=axis)
    legend!(fig[1,1], grid; tellwidth=false, halign=:right, valign=:top, margin=(10, 10, 10, 10), patchsize=(20,20))
    return fig, filename
end