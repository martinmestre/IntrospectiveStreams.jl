"""Run basic_pipeline for each stream."""

name_s = ["GD-1", "Pal5", "Jhelum", "PS1-A", "Fjorm-M68"]
name_t = ["GD-1-I21", "Pal5-PW19", "Jhelum-I21", "PS1-A-B16", "M68-P19"]

for i in eachindex(name_s)
    basic_pipeline(name_s[i],name_t[i])
end