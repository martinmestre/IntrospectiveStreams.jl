
"""Example of package usage."""
function example()
    Makie.inline!(true)
    name_s = ["GD-1", "GD-1", "Pal5", "Jhelum", "Fjorm-M68","PS1-A"][2:3]
    name_t = ["GD-1-PB18", "GD-1-PB18", "Pal5-PW19", "Jhelum-b-B19", "M68-P19", "PS1-A-B16"][2:3]
    age    = [12.0, 12.0, 12.0, 12.0 ,11.2, 12.0][2:3]*10^9 #yr
    metal  = [-1.5, -1.5, -1.4, -1.2, -2.2, -1.7][2:3]
    filters = fill("UBVRIplus",length(name_s))
    dr     = fill("DR3",length(name_s))
    # dr[1] = "DR2"
    tol_curation = [0.3, 0.3, 1.0]  # tolerances in μ_α*cosδ, μ_δ, Π.
    col_bounds = (-1.0, 4.0)
    box_μ = [[-14,-10.],[-4.,-2.]]
    σ_c = 1
    σ = 0.7

    file_orig, file_corr, file_phot, file_iso, file_filt, file_plot  = name_files_all(dr, name_s, age, metal)
    # correct_extinction_Gaia_loop(name_s, file_orig, file_corr)
    basic_pipeline_loop(name_s, name_t, file_corr, file_phot, file_gc, file_iso, file_filt, file_plot,
                        age, metal,  filters, tol_curation, col_bounds, box_μ, σ_c, σ)

    println("Tasks performed.")
end