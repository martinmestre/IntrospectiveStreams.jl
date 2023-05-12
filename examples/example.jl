using Pkg
Pkg.activate(".")
using IntrospectiveStreams


name_s = ["GD-1", "Pal5", "PS1-A", "Jhelum", "Fjorm-M68"]
name_t = ["GD-1-I21", "Pal5-PW19", "Jhelum-I21", "PS1-A-B16", "M68-P19"]

file_orig, file_corr, file_phot = name_files_all_cats(name_s)

# correct_extinction_Gaia_loop(name_s, file_orig, file_corr)

df = multi_pipeline(file_corr, file_phot, name_t)

println("Tasks performed.")