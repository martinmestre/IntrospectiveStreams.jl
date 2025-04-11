using DataFrames
using CSV
using ScatteredInterpolation
using CairoMakie
using AlgebraOfGraphics
using Statistics  # Para usar var()


"""
    load_isochrone_data(filepath::String)

Load isochrone data from a file.
"""
function load_isochrone_data(filepath::String)
    df = DataFrame(CSV.File(filepath, delim=" ", ignorerepeated=true, comment="#", silencewarnings=true))
    return df
end

"""Remove columns with any missing value."""
function eliminar_columnas_con_missing!(df::DataFrame)
    n_antes = ncol(df)
    cols_con_missing = [col for col in names(df) if any(ismissing.(df[!, col]))]

    if !isempty(cols_con_missing)
        println("Eliminando columnas con missing: ", join(cols_con_missing, ", "))
        select!(df, Not(cols_con_missing))
        println("Se eliminaron ", n_antes - ncol(df), " columnas")
    else
        println("No se encontraron columnas con missing")
    end

    return df  # Devuelve el mismo DataFrame modificado
end

"""
    _bracket(value::Real, sorted_seq::AbstractVector{<:Real})

Find the interval in sorted_seq that brackets the given value.
Returns a vector with 1 or 2 elements.
"""
function _bracket(value::Real, sorted_seq::AbstractVector{<:Real})
    if value ≤ first(sorted_seq)
        return [first(sorted_seq)]
    elseif value ≥ last(sorted_seq)
        return [last(sorted_seq)]
    end

    idx = searchsortedfirst(sorted_seq, value)
    return unique([sorted_seq[idx-1], sorted_seq[idx]])
end

"""
    _add_evolution_phase(iso::DataFrame)

Add continuous evolution phase values to an isochrone DataFrame.
"""
function _add_evolution_phase(iso::DataFrame)
    iso = copy(iso)
    iso.evol = zeros(nrow(iso))

    for label in unique(iso.label)
        idx = iso.label .== label
        n_points = sum(idx)
        delta = 1.0 / n_points
        evol_values = label .+ (0:delta:(1-0.5*delta))[1:n_points]
        iso.evol[idx] .= evol_values
    end

    return iso
end

"""
    _get_interp_points(iso::DataFrame)

Get interpolation points (logAge, MH, evol) from an isochrone DataFrame.
"""
function _get_interp_points(iso::DataFrame)
    hcat(iso.logAge, iso.MH, iso.evol)
end

"""
    interpolate_isochrone(family::Symbol, filter::String, age::Real, metal::Real)

High-level interface to interpolate isochrones from the PARSEC database.
"""
function interpolate_isochrone(family::Symbol, filter::String, age::Real, metal::Real)
    @assert family == :parsec "Only PARSEC isochrones are currently supported"
    @assert 9.2 ≤ log10(age) ≤ 10.3 "Age should fulfill: 9.2 ≤ log10(age) ≤ 10.3"
    @assert -2.2 < metal ≤ 0.5 "Metallicity should satisfy: -2.2 < [Fe/H] ≤ 0.5"

    if filter == "hsc"
        file_artif = "artifacts/isochrones/parsec/$filter/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
        df_artif = load_isochrone_data(file_artif)
        return interpolate_isochrone(df_artif, age, metal)
    else
        error("Filter $filter not implemented")
    end
end
"""
    interpolate_isochrone(df::DataFrame, target_age::Real, target_metallicity::Real)

Interpolate isochrones at given (age, metallicity) from a DataFrame of PARSEC isochrones.

# Arguments
- `df`: DataFrame containing isochrone data with columns "logAge", "MH", and "label"
- `target_age`: Age in years (linear scale) to interpolate to
- `target_metallicity`: Metallicity [Fe/H] to interpolate to

# Returns
- DataFrame with interpolated isochrone at requested parameters
"""

function interpolate_isochrone(df::DataFrame, target_age::Real, target_metallicity::Real)

    int_columns = [col for col in names(df) if eltype(df[!, col]) <: Integer]
    println("Columnas enteras: ", int_columns)
    # --- 1. Preprocesamiento seguro ---
    target_logAge = log10(target_age)

    # Eliminar solo columnas no esenciales con missing (excepto las de interpolación)
    essential_cols = ["logAge", "MH", "label"]
    df = _add_evolution_phase(df)  # Añade columna 'evol' continua
    non_essential = setdiff(names(df), [essential_cols; "evol"])
    cols_to_drop = [col for col in non_essential if any(ismissing, df[!, col])]
    if !isempty(cols_to_drop)
        select!(df, Not(cols_to_drop))
    end

    # --- 2. Preparación de coordenadas ---
    logAges = sort(unique(df.logAge))
    MHs = sort(unique(df.MH))
    columns = setdiff(names(df), ["logAge", "MH", "label", "evol", "pmode"])  # Propiedades a interpolar

    # --- 3. Recolección de puntos CON 'evol' ---
    interp_points = []
    interp_values = Dict(col => [] for col in columns)

    for logAge in _bracket(target_logAge, logAges), MH in _bracket(target_metallicity, MHs)
        iso = filter(row -> isapprox(row.logAge, logAge, rtol=1e-3) &&
                          isapprox(row.MH, MH, rtol=1e-3), df)
        isempty(iso) && continue

        # Usamos evol (no label) como tercera coordenada
        points = hcat(iso.logAge, iso.MH, iso.evol)  # ¡Clave! 3 columnas: logAge, MH, evol
        push!(interp_points, points)

        for col in columns
            push!(interp_values[col], iso[!, col])
        end
    end

    isempty(interp_points) && error("No interpolation points found for target age and metallicity")

    # --- 4. Interpolación por propiedad ---
    all_points = vcat(interp_points...)
    phases = range(0, 9, length=1000)
    query_points = hcat(
        fill(target_logAge, length(phases)),
        fill(target_metallicity, length(phases)),
        phases
    )'  # Mantenemos la transpuesta para la interpolación

    # Crear el DataFrame con las columnas fijas primero
    evol_phases = collect(phases)
    result_df = DataFrame(
        logAge = fill(target_logAge, length(phases)),
        MH = fill(target_metallicity, length(phases)),
        evol = evol_phases,  # Convertido a Vector explícitamente
        label = floor.(Int, evol_phases)
    )

    for col in columns
        # Concatenar todos los valores para esta columna
        all_values = vcat(interp_values[col]...)

        # Verificar dimensiones
        if length(all_values) != size(all_points, 1)
            @warn "Skipping column $col due to length mismatch (expected $(size(all_points,1)), got $(length(all_values)))"
            continue
        end

        # Crear interpolación
        itp = interpolate(Shepard(), all_points', all_values)

        # Único cambio esencial: usar vec() para convertir el resultado matricial a vector
        result_df[!, col] = vec(evaluate(itp, query_points))
    end

    return result_df
end




"""
    plot_isochrone_cmd(df::DataFrame, telescope::Symbol, file::String)

Plot the isochrone in CMD space for the specified telescope.
"""
function plot_isochrone_cmd(df::DataFrame, telescope::Symbol, file::String)
    # Configuración básica del gráfico
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1],
              xlabel="G - R",
              ylabel="G",
              title="Isochrone CMD")

    # Mapeo simple de colores por fase evolutiva
    colors = cgrad(:viridis, length(unique(df.label)), categorical=true)

    # Graficar cada fase evolutiva por separado
    for (i, label) in enumerate(sort(unique(df.label)))
        subset = df[df.label .== label, :]
        lines!(ax, subset.gmag - subset.rmag, subset.gmag,
               color=colors[i], label="Phase $label")
    end

    # Leyenda y guardado
    Legend(fig[1, 2], ax)
    save(file, fig)
    return fig
end

"""
    example_interpolate_isochrone()

Example usage with Subaru HSC filters.
"""
function example_interpolate_isochrone()
    camera = "hsc"
    file_artif = "artifacts/isochrones/parsec/$camera/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
    df_artif = load_isochrone_data(file_artif)

    file_plot = "plots/isochrone_$camera.pdf"
    target_age = 10.0^9.5
    target_metallicity = 0.0  # Solar metallicity

    # Interpolate and plot
    isochrone_df = interpolate_isochrone(df_artif, target_age, target_metallicity)
    fig = plot_isochrone_cmd(isochrone_df, :Subaru, file_plot)

    return fig
end

