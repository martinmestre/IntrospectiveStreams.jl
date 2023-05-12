using Pkg
Pkg.activate(".")
using IntrospectiveStreams


name_s = ["GD-1", "Pal5", "PS1-A", "Jhelum", "Fjorm-M68"]
name_t = ["GD-1-PB18", "Pal5-PW19", "Jhelum-b-B19", "PS1-A-B16", "M68-P19"]
age    = [...,11.2e9]
metal  = [...,-2.2]
filters = fill("UBVRIplus",length(name_s))
tol_curation = [0.3, 0.3, 0.3]  # for all the streams the same tolerances in μ_α*cosδ, μ_δ, Π.
col_bounds = (-0.75, 2.0)
σ_c = 0.1
σ = 0.2

file_orig, file_corr, file_phot, file_iso  = name_files_all(name_s)

# correct_extinction_Gaia_loop(name_s, file_orig, file_corr)

basic_pipeline_loop(name_t, file_corr, file_phot, file_gc, file_filt, ages, metals, filters, tol_curation, col_bounds, σ_c, σ)

println("Tasks performed.")