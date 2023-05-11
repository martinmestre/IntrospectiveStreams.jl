"""Example about using this package."""
using Pkg
Pkg.activate(".")
using IntrospectiveStreams

name_s = ["GD-1", "Pal5", "PS1-A", "Jhelum", "Fjorm-M68"]

file_orig, file_corr, file_phot = name_files_all(name_s)

# correct_extinction_Gaia_loop(name_s, file_orig, file_corr)

println("Tasks performed.")