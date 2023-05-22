### A Pluto.jl notebook ###
# v0.19.25

using Markdown
using InteractiveUtils

# ╔═╡ aa0728bb-dd20-498b-a470-831e238e1bb2
begin
	using Pkg; Pkg.activate("../")
	using Revise
	using IntrospectiveStreams
	using CSV
	using DataFrames
	using CairoMakie, AlgebraOfGraphics
end

# ╔═╡ f6dc3d94-f68d-11ed-1495-c11c11bbb65c
md"# Community Atlas of Stellar Streams"

# ╔═╡ 3cc2a558-622c-4a5a-8085-557c69bac9a8
Makie.inline!(true)

# ╔═╡ 5c7cee94-5226-4376-853d-3f10ba846054
md"## Searching the streams at first order (using PM track from Galstreams)
Here we show the result of making approximate cuts in proper motion (PM) and color-magnitude (CM) for the five streams using Gaia DR3. We also show the case of GD-1 with DR2 for comparison.
The PM cut is guided by the Galstream track. The filtering is done by keeping only stars that satify:
$\forall \phi_1: |\mu_1^{(t)}-\mu_1(\phi_1^{(t)}|<0.7$ and
$\forall \phi_1: |\mu_2^{(t)}-\mu_2(\phi_1^{(t)}|<0.7$
where $\mu_{1,2}$ are the PMs not corrected by the solar reflex-motion.
The process is performed by the function <example> defined below, where the input parameters can be known.
"

# ╔═╡ d0252793-efdb-455b-a23a-02ff8896b34d
begin
    name_s = ["GD-1", "GD-1", "Pal5", "Jhelum", "Fjorm-M68","PS1-A"]
    name_t = ["GD-1-PB18", "GD-1-PB18", "Pal5-PW19", "Jhelum-b-B19", "M68-P19", "PS1-A-B16"]
    age    = [12.0, 12.0, 12.0, 12.0 ,11.2, 12.0]*10^9 #yr
    metal  = [-1.5, -1.5, -1.4, -1.2, -2.2, -1.7]
    filters = fill("UBVRIplus",length(name_s))
    dr     = fill("DR3",length(name_s))
    # dr[1] = "DR2"
    tol_curation = [0.3, 0.3, 1.0]  # tolerances in μ_α*cosδ, μ_δ, Π.
    col_bounds = (-1.0, 4.0)
    box_μ = [[-14,-10.],[-4.,-2.]]
    σ_c = 1
    σ = 0.7

    file_orig, file_corr, file_phot, file_iso, file_filt, file_plot  = name_files_all(dr, name_s, age, metal)
end

# ╔═╡ f5c720f8-d807-4e0b-98f6-252072bdedff
md"Now lets reopen the filtered files and display the plots"

# ╔═╡ cff0f3bf-b79e-40a1-a944-8cd290abcf71
begin
	figs = Vector{Figure}(undef, length(name_s))
	for i ∈ eachindex(name_s)
		df = DataFrame(CSV.File(file_filt[i], delim=",", ignorerepeated=true))
		figs[i]=plot_scatter_on_sky_self_frame(name_s[i], df, file_plot[i])
	end
end

# ╔═╡ 5a70437e-8b9a-41db-91de-91c924832eca
figs[1]

# ╔═╡ 8051526f-566b-48a5-a0d4-43bea8dff698
figs[2]

# ╔═╡ Cell order:
# ╟─f6dc3d94-f68d-11ed-1495-c11c11bbb65c
# ╠═aa0728bb-dd20-498b-a470-831e238e1bb2
# ╠═3cc2a558-622c-4a5a-8085-557c69bac9a8
# ╠═5c7cee94-5226-4376-853d-3f10ba846054
# ╠═d0252793-efdb-455b-a23a-02ff8896b34d
# ╟─f5c720f8-d807-4e0b-98f6-252072bdedff
# ╠═cff0f3bf-b79e-40a1-a944-8cd290abcf71
# ╠═5a70437e-8b9a-41db-91de-91c924832eca
# ╠═8051526f-566b-48a5-a0d4-43bea8dff698
