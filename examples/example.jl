
"""Example of package usage."""
function example()

    name_s = ["GD-1", "Pal5", "PS1-A", "Jhelum", "Fjorm-M68"]
    name_t = ["GD-1-PB18", "Pal5-PW19", "Jhelum-b-B19", "PS1-A-B16", "M68-P19"]
    age    = [12.0, 12.0, 12.0, 12.0 ,11.2]*10^9 #yr
    metal  = [-1.5, -1.4, -1.7, -1.2, -2.2]
    name_s = ["GD-1"]
    name_t = ["GD-1-PB18"]
    age    = [12.0]
    metal  = [-1.5]
    filters = fill("UBVRIplus",length(name_s))
    tol_curation = [0.3, 0.3, 0.3]  # for all the streams the same tolerances in μ_α*cosδ, μ_δ, Π.
    col_bounds = (-0.75, 2.0)
    σ_c = 1
    σ = 2

    file_orig, file_corr, file_phot, file_iso, file_filt, file_plot  = name_files_all(name_s, age, metal)
    # correct_extinction_Gaia_loop(name_s, file_orig, file_corr)
    basic_pipeline_loop(name_t, file_corr, file_phot, file_gc, file_iso, file_filt, file_plot, age, metal,  filters, tol_curation, col_bounds, σ_c, σ)

    println("Tasks performed.")
end