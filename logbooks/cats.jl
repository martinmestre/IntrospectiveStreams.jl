### A Pluto.jl notebook ###
# v0.19.24

using Markdown
using InteractiveUtils

# ╔═╡ aa0728bb-dd20-498b-a470-831e238e1bb2
begin
	using Pkg; Pkg.activate("../")
	using IntrospectiveStreams
	using CSV
	using DataFrames
	using CairoMakie
end

# ╔═╡ f6dc3d94-f68d-11ed-1495-c11c11bbb65c
md"# Community Atlas of Stellar Streams"

# ╔═╡ 3cc2a558-622c-4a5a-8085-557c69bac9a8
Makie.inline!(true)

# ╔═╡ 5c7cee94-5226-4376-853d-3f10ba846054
md"""## Searching the streams at first order
Here we show the result of making approximate cuts in proper motion (PM) and color-magnitude (CM) for the five streams using Gaia DR3. We also show the case of GD-1 with DR2 for comparison. 
The PM cut is guided by the Galstream track: $\mu^{(t)}$. The filtering is done by keeping only stars that satify:
$\forall \phi_1: |\mu_1^{(⋆)}-\mu_1^{(t)}(\phi_1^{(⋆)})|<0.7$ and
$\forall \phi_1: |\mu_2^{(⋆)}-\mu_2^{(t)}(\phi_1^{(⋆)})|<0.7$
where $\mu_{1,2}$ are the PMs not corrected by the solar reflex-motion, and ⋆ denotes the Gaia data.
The process is performed by the function "example" given below. There you can see the input parameters used for age, metal, $\sigma_c$ (width in color) and $\sigma$ (width in PM). 
"""

# ╔═╡ d0252793-efdb-455b-a23a-02ff8896b34d
begin
    name_s = ["GD-1", "GD-1", "Pal5", "Jhelum", "Fjorm-M68","PS1-A"]
    name_t = ["GD-1-PB18", "GD-1-PB18", "Pal5-PW19", "Jhelum-b-B19", "M68-P19", "PS1-A-B16"]
    age    = [12.0, 12.0, 12.0, 12.0 ,11.2, 12.0]*10^9 #yr
    metal  = [-1.5, -1.5, -1.4, -1.2, -2.2, -1.7]
    filters = fill("UBVRIplus",length(name_s))
    dr     = fill("DR3",length(name_s))
    dr[1] = "DR2"
    tol_curation = [0.3, 0.3, 1.0]  # tolerances in μ_α*cosδ, μ_δ, Π.
    col_bounds = (-1.0, 4.0)
    box_μ = [[-14,-10.],[-4.,-2.]]
    σ_c = 1
    σ = 0.7

    file_orig, file_corr, file_phot, file_iso, file_filt, file_plot  = name_files_all(dr, name_s, age, metal)
	nothing
end

# ╔═╡ f5c720f8-d807-4e0b-98f6-252072bdedff
md"Now lets reopen the filtered files and display the plots without the tracks so as not to cheat the eye."

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

# ╔═╡ 1e581678-aef1-4ae1-a45d-e157266bc0f6
figs[3]

# ╔═╡ 25452558-40b7-4626-a1c8-2d03b493f9fd
figs[4]

# ╔═╡ d779f4e4-0f93-4183-9c01-1d4009f8db8e
figs[5]

# ╔═╡ b6b10200-f743-479d-8aee-668d6938adc3
figs[6]

# ╔═╡ ffbea6cf-7ddd-4449-abe1-52e309c72ca4
begin
	figs_t = Vector{Figure}(undef, length(name_s))
	for i ∈ eachindex(name_s)
		df = DataFrame(CSV.File(file_filt[i], delim=",", ignorerepeated=true))
		df_track, self_frame = load_stream_track(name_t[i]) 
		figs_t[i]=plot_scatter_on_sky_self_frame(name_s[i], df, df_track, file_plot[i])
	end
end

# ╔═╡ e60d7cfd-75fd-411e-8903-6ba689cdf62c
md"""Now lets plot together with the tracks."""

# ╔═╡ b9fc214d-480a-4123-96c6-ac02ffab9d54
figs_t[1]

# ╔═╡ d6f2db68-b4cf-4ef0-a571-f83fb058e4ce
figs_t[2]

# ╔═╡ b472766e-d3db-4231-9653-2f16a9cc7c04
figs_t[3]

# ╔═╡ eafa11bd-09c9-4d24-a7d2-d87e07aefef9
figs_t[4]

# ╔═╡ 1b249ea2-1109-4bdf-9c38-dc8dc8a85054
figs_t[5]

# ╔═╡ 16f8ab26-24ed-40a3-92c7-df8adc5d8464
figs_t[6]

# ╔═╡ Cell order:
# ╟─f6dc3d94-f68d-11ed-1495-c11c11bbb65c
# ╠═aa0728bb-dd20-498b-a470-831e238e1bb2
# ╠═3cc2a558-622c-4a5a-8085-557c69bac9a8
# ╟─5c7cee94-5226-4376-853d-3f10ba846054
# ╠═d0252793-efdb-455b-a23a-02ff8896b34d
# ╟─f5c720f8-d807-4e0b-98f6-252072bdedff
# ╠═cff0f3bf-b79e-40a1-a944-8cd290abcf71
# ╟─5a70437e-8b9a-41db-91de-91c924832eca
# ╟─8051526f-566b-48a5-a0d4-43bea8dff698
# ╟─1e581678-aef1-4ae1-a45d-e157266bc0f6
# ╟─25452558-40b7-4626-a1c8-2d03b493f9fd
# ╟─d779f4e4-0f93-4183-9c01-1d4009f8db8e
# ╟─b6b10200-f743-479d-8aee-668d6938adc3
# ╠═ffbea6cf-7ddd-4449-abe1-52e309c72ca4
# ╟─e60d7cfd-75fd-411e-8903-6ba689cdf62c
# ╠═b9fc214d-480a-4123-96c6-ac02ffab9d54
# ╠═d6f2db68-b4cf-4ef0-a571-f83fb058e4ce
# ╠═b472766e-d3db-4231-9653-2f16a9cc7c04
# ╠═eafa11bd-09c9-4d24-a7d2-d87e07aefef9
# ╠═1b249ea2-1109-4bdf-9c38-dc8dc8a85054
# ╠═16f8ab26-24ed-40a3-92c7-df8adc5d8464
