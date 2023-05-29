# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: -all
#     custom_cell_magics: kql
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.2
#   kernelspec:
#     display_name: Julia 1.9.0
#     language: julia
#     name: julia-1.9
# ---

# %% [markdown]
# # Re-discovering Fjörm stellar stream

# %%
using Pkg; Pkg.activate(".")
using IntrospectiveStreams
using CSV
using FITSIO
using DataFrames, DataFramesMeta
using CairoMakie
using Interpolations
# %%
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

i = 5
name_s = name_s[i]
name_t = name_t[i]
file_corr = file_corr[i]
file_plot = file_plot[i]
age = age[i]
metal = metal[i]
filters = filters[i]
dr = dr[i]
# %%

# Load Gaia data and gc table.
f = FITS(file_corr)
df = DataFrame(f[2])
df_gc  = DataFrame(CSV.File(file_gc, delim=" ", ignorerepeated=true))
# Remove known globular clusters.
rename!(df_gc,[:RA, :DEC] .=> [:ra, :dec])
mask_gc!(df, df_gc)
# Load galstreams data.
df_track, self_frame = load_stream_track(name_t)
D_interp = linear_interpolation(df_track.ϕ₁, df_track.D)
# %%
# Compute the stream self-coordinates and correcte for reflex-motion of the ⊙.
compute_in_self_coords!(df, self_frame)
@subset!(df, minimum(df_track.ϕ₁) .< :ϕ₁ .< maximum(df_track.ϕ₁))
df.D = D_interp.(df.ϕ₁)
reflex_correct!(df, self_frame)
# %%

fig = plot_scatter_on_sky_self_frame(name_s, df, file_plot)


