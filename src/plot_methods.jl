"""Plot sky histogram."""
function plot_sky_histo(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")*
                            histogram(bins=200)
    ag=draw!(fig, plt, axis=(;limits=((nothing,nothing),(nothing,nothing))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

function plot_sky_histo_gc(df::DataFrame, df_gc::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=200)+data(df_gc)*visual(color="red"))*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")
    plt_M68 = data(df_gc)*mapping(:ra =>L"RA [$°$]", :dec=>L"Dec [$°$]")*visual(color="black")
    ag = draw!(fig, plt, axis=(;limits=((180,270),(-30,80))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot sky histogram in stream frame."""
function plot_sky_histo_selfFrame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*AlgebraOfGraphics.density()+data(df_track))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
                            # histogram(bins=100)
    ag = draw!(fig, plt, axis=(;limits=((nothing,nothing),(nothing,nothing))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot sky scatter in stream frame."""
function plot_sky_scatter_selfFrame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (6*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black,1))+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    # plt = data(df)*visual(markersize=0.7)*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt, axis=(;limits=((nothing,nothing),(-10,2))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot sky with arrows for μ in stream frame."""
function plot_sky_scatter_μ_arrows_selfFrame(df::DataFrame, df_track::DataFrame, file::String)
    step = 500
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=3)+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt) #, axis=(;limits=((-20,20),(-20,20))))
    us = df.μ₁cosϕ₂./cos.(df.ϕ₂*π/180.)
    vs = df.μ₂
    strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
    us .= us #./strength
    vs .= vs #./strength
    arrows!(df.ϕ₁, df.ϕ₂, us, vs, arrowsize = (3/4)strength, lengthscale = 1)
    us = (df_track.μ₁cosϕ₂./cos.(df_track.ϕ₂*π/180.))[begin:step:end]
    vs = df_track.μ₂[begin:step:end]
    strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
    us .= us #./strength
    vs .= vs #./strength
    arrows!(df_track.ϕ₁[begin:step:end], df_track.ϕ₂[begin:step:end], us, vs, arrowsize = 2strength, lengthscale = 1, color="blue")
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot sky with arrows for reflex-corrected μ in stream frame."""
function plot_sky_scatter_μ_arrows_corr_selfFrame(df::DataFrame, df_track::DataFrame, file::String)
    step = 300
    size_inches = (5*3, 2*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=3)+data(df_track)*visual(markersize=1,color="red"))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    plt = (data(df)*visual(markersize=3))*mapping(:ϕ₁ =>L"ϕ_1 [°]", :ϕ₂=>L"ϕ_2 [°]")
    ag = draw!(fig, plt)#, axis=(;limits=((0,30),(-10,10))))
    us = df.μ₁_corr
    vs = df.μ₂_corr
    strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
    arrows!(df.ϕ₁, df.ϕ₂, us, vs, arrowsize = strength, lengthscale = 0.2)
    us = df_track.μ₁_corr[begin:step:end]
    vs = df_track.μ₂_corr[begin:step:end]
    strength = vec(sqrt.(us .^ 2 .+ vs .^ 2))
    arrows!(df_track.ϕ₁[begin:step:end], df_track.ϕ₂[begin:step:end], us, vs, arrowsize = strength, lengthscale = 0.5, color="cyan")
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot CMD histogram."""
function plot_cmd_histo(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*mapping(:color=>L"$BP-RP$ [Mag]", :g_abs=>L"$G$ [Mag]")*
                                    histogram(bins=300)
    ag = draw!(fig, plt, axis=(;yreversed=true))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot single isochrone."""
function plot_isochrone(df::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*mapping(:color=>L"BP-RP", :Gaia_G_EDR3 =>L"G")*visual(Lines)
    draw!(fig, plt, axis=(;yreversed=true))
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot isochone plus data."""
function plot_isochrone_data(df_iso::DataFrame, df_s::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt_s = data(df_s)*mapping(:color=>L"$BP-RP$ [Mag]", :g_abs=>L"$G$ [Mag]")*histogram(bins=100)
    plt_iso = data(df_iso)*mapping(:color=>L"$BP-RP$ [Mag]", :Gaia_G_EDR3=>L"$G$ [Mag]")*visual(Lines,color="red")
    plt_iso_bord = data(df_iso)*mapping([:left,:right].=>L"$BP-RP$ [Mag]", :Gaia_G_EDR3=>L"$G$ [Mag]")*visual(Lines,color="black")
    ag = draw!(fig, plt_s+plt_iso+plt_iso_bord, axis=(;yreversed=true))#, limits=((0,1.5),(14,22))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot μ space."""
function plot_μ(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*mapping(:pmra =>L"$μ_{RA}$ [mas/yr]", :pmdec=>L"$μ_{Dec}$ [mas/yr]")*
                            histogram(bins=500)
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot μ space inside a window."""
function plot_μ_window(df::DataFrame, window::Vector{Vector{Float64}}, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*mapping(:pmra =>L"$μ_{RA}$ [mas/yr]", :pmdec=>L"$μ_{Dec}$ [mas/yr]")*
                            histogram(bins=2000)
    ag = draw!(fig, plt, axis=(;limits=((window[1][1], window[1][2]),(window[2][1],window[2][2]))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
end

"""Plot μ space in stream's self-frame."""
function plot_μ_selfFrame(df::DataFrame, file::String)
    size_inches = (5*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = data(df)*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")*
                            histogram(bins=500)
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Plot μ space in stream's self-frame inside a window."""
function plot_μ_selfFrame_window(df::DataFrame, window::Vector{Vector{Float64}}, file::String)
    size_inches = (6*3, 5*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*AlgebraOfGraphics.density())*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=((window[1][1], window[1][2]),(window[2][1],window[2][2]))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end
"""Plot μ space in stream's self-frame inside a window plus track."""
function plot_μ_selfFrame_window(df::DataFrame, df_track, window::Vector{Vector{Float64}}, file::String)
    size_inches = (4*3, 4*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*histogram(bins=2000)+data(df_track)*visual(markersize=1))*mapping(:μ₁cosϕ₂ =>L"$μ_1cosϕ_2$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=((window[1][1], window[1][2]),(window[2][1],window[2][2]))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Scatter plot μ space in stream's self-frame inside a window."""
function plot_μ_scatter_selfFrame_window(df::DataFrame, df_track::DataFrame, window::Vector{Vector{Float64}}, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=5)+data(df_track)*visual(markersize=0.5,color="red"))*mapping(:μ₁ =>L"$μ_1$ [mas/yr]", :μ₂=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=((window[1][1], window[1][2]),(window[2][1],window[2][2]))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Scatter plot corrected-μ space in stream's self-frame."""
function plot_μ_corr_scatter_selfFrame(df::DataFrame, df_track::DataFrame, file::String)
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=1, color=(:black, 1))+data(df_track)*visual(markersize=0.5,color="red"))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt)
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Scatter plot corrected-μ space in stream's self-frame inside a window."""
function plot_μ_corr_scatter_selfFrame_window(df::DataFrame, df_track::DataFrame, file::String,  window::Vector{Vector{Float64}})
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*visual(markersize=10, color=(:black, 0.5))+data(df_track)*visual(markersize=0.5,color="red"))*mapping(:μ₁_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=((window[1][1], window[1][2]),(window[2][1],window[2][2]))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Histogram plot corrected-μ space in stream's self-frame inside a window."""
function plot_μ_corr_histo_selfFrame_window(df::DataFrame, df_track::DataFrame, file::String,  window::Vector{Vector{Float64}})
    size_inches = (3*3, 3*3)
    size_pt = 72 .* size_inches
    fig = Figure(resolution = size_pt, fontsize = 30)
    plt = (data(df)*AoG.density()+data(df_track)*visual(markersize=0.5,color="red"))*mapping(:μ₁cosϕ₂_corr =>L"$μ_1$ [mas/yr]", :μ₂_corr=>L"$μ_2$ [mas/yr]")
    ag = draw!(fig, plt, axis=(;limits=((window[1][1], window[1][2]),(window[2][1],window[2][2]))))
    colorbar!(fig[1,2], ag)
    electrondisplay(fig)
    save(file, fig, pt_per_unit=1)
    return nothing
end

"""Scatter plot corrected-μ track in stream's self-frame."""
function plot_μ_corr_track_selfFrame(df_track::DataFrame, file::String)
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


