using DataFrames
using CSV
using Interpolations
using LinearAlgebra

# Cargar el DataFrame desde un archivo
filter = "G"  # Filtro deseado
file_artif = "artifacts/isochrones/parsec/$filter/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
df_artif = DataFrame(CSV.File(file_artif, delim=" ", ignorerepeated=true, comment="#"))

# Función para encontrar las 4 isócronas más cercanas usando mapeo de la grilla
function find_nearest_isochrones(df, target_age, target_metallicity)
    # Verificar que la edad y metalicidad objetivo estén dentro del rango
    age_min, age_max = extrema(df.age)
    metallicity_min, metallicity_max = extrema(df.metallicity)

    if target_age < age_min || target_age > age_max ||
       target_metallicity < metallicity_min || target_metallicity > metallicity_max
        error("La edad o metalicidad objetivo están fuera del rango para interpolar.")
    end

    # Mapear la grilla de edad y metalicidad a índices enteros
    unique_ages = sort(unique(df.age))
    unique_metallicities = sort(unique(df.metallicity))

    # Encontrar los índices de los nodos más cercanos
    age_idx = searchsortedfirst(unique_ages, target_age)
    metallicity_idx = searchsortedfirst(unique_metallicities, target_metallicity)

    # Obtener las cuatro combinaciones de edad y metalicidad más cercanas
    age_nodes = unique_ages[age_idx-1:age_idx]
    metallicity_nodes = unique_metallicities[metallicity_idx-1:metallicity_idx]

    # Crear un DataFrame con las cuatro isócronas más cercanas
    nearest_isochrones = df[
        (df.age .∈ [age_nodes]) .& (df.metallicity .∈ [metallicity_nodes]), :
    ]

    return nearest_isochrones
end

# Función para interpolar en masa inicial para una isócrona específica
function interpolate_in_initial_mass(iso, initial_mass)
    # Ordenar por masa inicial
    sort!(iso, :initial_mass)

    # Crear el interpolador en masa inicial
    itp = LinearInterpolation(iso.initial_mass, iso.magnitude, extrapolation_bc=Line())

    return itp(initial_mass)
end

# Función para generar la isócrona completa
function generate_complete_isochrone(df, target_age, target_metallicity)
    # Encontrar las 4 isócronas más cercanas
    nearest_isochrones = find_nearest_isochrones(df, target_age, target_metallicity)

    # Obtener los rangos de masa inicial de las isócronas cercanas
    initial_mass_ranges = [unique(df[(df.age .== row.age) .& (df.metallicity .== row.metallicity), :initial_mass]) for row in eachrow(nearest_isochrones)]

    # Calcular la intersección de los rangos de masas iniciales
    common_initial_masses = intersect(initial_mass_ranges...)

    # Seleccionar la isócrona con más valores de masas iniciales en el rango común
    iso_with_most_masses = df[(df.age .== nearest_isochrones.age[1]) .& (df.metallicity .== nearest_isochrones.metallicity[1]), :]
    for row in eachrow(nearest_isochrones)
        iso = df[(df.age .== row.age) .& (df.metallicity .== row.metallicity), :]
        if length(intersect(iso.initial_mass, common_initial_masses)) > length(intersect(iso_with_most_masses.initial_mass, common_initial_masses))
            iso_with_most_masses = iso
        end
    end

    # Usar el espaciamiento de la isócrona seleccionada
    selected_masses = sort(iso_with_most_masses.initial_mass)
    selected_spacing = diff(selected_masses)

    # Generar el rango de masas iniciales con el espaciamiento seleccionado
    min_mass = minimum(common_initial_masses)
    max_mass = maximum(common_initial_masses)
    initial_masses = collect(min_mass:selected_spacing[1]:max_mass)

    # Interpolar todos los campos del DataFrame (excepto age, metallicity, initial_mass)
    interpolated_df = DataFrame(initial_mass=initial_masses)
    for col in names(df)
        if col ∉ ["age", "metallicity", "initial_mass"]
            magnitudes = []
            for initial_mass in initial_masses
                # Interpolar en masa inicial para cada isócrona cercana
                mags = []
                for row in eachrow(nearest_isochrones)
                    iso = df[(df.age .== row.age) .& (df.metallicity .== row.metallicity), :]
                    mag = interpolate_in_initial_mass(iso, initial_mass)
                    push!(mags, mag)
                end

                # Interpolación bilineal en edad y metalicidad usando Interpolations.jl
                age_nodes = unique(nearest_isochrones.age)
                metallicity_nodes = unique(nearest_isochrones.metallicity)
                itp = LinearInterpolation((age_nodes, metallicity_nodes), reshape(mags, (2, 2)), extrapolation_bc=Line())
                magnitude = itp(target_age, target_metallicity)
                push!(magnitudes, magnitude)
            end
            interpolated_df[!, col] = magnitudes
        end
    end

    return interpolated_df
end

# Ejemplo de uso
target_age = 1.5  # Edad deseada
target_metallicity = 0.015  # Metalicidad deseada

# Generar la isócrona completa
isochrone_df = generate_complete_isochrone(df_artif, target_age, target_metallicity)
println("Isócrona completa:")
println(isochrone_df)