# The four metallicity-abundance transformations below were taken from the
# ugali.py code (https://github.com/DarkEnergySurvey/ugali).
"""Z to FeH transformation for MESA isochrones (MIST)"""
function z2feh_mist(z)
        # Section 3.1 of Choi et al. 2016 (https://arxiv.org/abs/1604.08592)
        z₀ = z                # Initial metal abundance
        yₚ = 0.249            # Primordial He abundance (Planck 2015)
        c  = 1.5              # He enrichment ratio

        y₀ = yₚ + c * z₀
        x₀ = 1 - y₀ - z₀

        z☼ = 0.0142           # Solar metal abundance
        y☼ = 0.2703           # Solar He abundance (Asplund 2009)
        x☼ = 1 - y☼ - z☼

        return log10( z₀/z☼ * x☼/x₀)
end

"""FeH to Z transformation for MESA isochrones (MIST)"""
function feh2z_mist(feh)
        # Section 3.1 of Choi et al. 2016 (https://arxiv.org/abs/1604.08592)
        yₚ = 0.249            # Primordial He abundance (Planck 2015)
        c  = 1.5              # He enrichment ratio

        z☼ = 0.0142           # Solar metal abundance
        y☼ = 0.2703           # Solar He abundance (Asplund 2009)
        x☼ = 1 - y☼ - z☼

        return (1 - yₚ)/( (1 + c) + (x☼/z☼) * 10^(-feh))
end

"""Z to FeH transformation for PADOVA isochrones (PARSEC)"""
function z2feh_parsec(z)
        # Taken from Table 3 and Section 3 of Bressan et al. 2012
        # Confirmed in Section 2.1 of Marigo et al. 2017
        z₀ = z                # Initial metal abundance
        yₚ = 0.2485            # Primordial He abundance (komatsu 2011)
        c  = 1.78             # He enrichment ratio

        y₀ = yₚ + c * z₀
        x₀ = 1 - y₀ - z₀

        z☼ = 0.01524           # Solar metal abundance
        y☼ = 0.2485           # Solar He abundance (Caffau 2011)
        x☼ = 1 - y☼ - z☼

        return log10( z₀/z☼ * x☼/x₀)
end

"""FeH to Z transformation for PADOVA isochrones (PARSEC)"""
function feh2z_parsec(feh)
        # Taken from Table 3 and Section 3 of Bressan et al. 2012
        # Confirmed in Section 2.1 of Marigo et al. 2017
        yₚ = 0.2485           # Primordial He abundance
        c  = 1.78             # He enrichment ratio

        z☼ = 0.01524          # Solar metal abundance
        y☼ = 0.2485           # Solar He abundance (Caffau 2011)
        x☼ = 1 - y☼ - z☼

        return (1 - yₚ)/( (1 + c) + (x☼/z☼) * 10^(-feh))
end


"""
    read_parsec_file(filename)
Deprecated use. Now files are saved in JLD2 format instead of CSV.
It can be use to see the files that ezpadova QuickInterpolator() opens.
"""
function read_parsec_file(filename)
    # Leer todas las líneas para procesar los comentarios
    lines = readlines(filename)

    # Buscar la última línea de comentario que parece contener nombres de columnas
    header_line = ""
    data_start = 0

    for (i, line) in enumerate(lines)
        if startswith(line, "#")
            # Eliminamos el # y espacios iniciales
            content = strip(replace(line, r"^#\s*" => ""))
            # Si la línea contiene al menos 3 palabras, asumimos que es el encabezado
            if length(split(content)) >= 3
                header_line = content
            end
        elseif data_start == 0
            data_start = i
            break  # Salimos al encontrar la primera línea de datos
        end
    end

    if isempty(header_line)
        error("No se encontró línea de encabezado válida en el archivo")
    end

    # Procesar los nombres de columnas
    column_names = split(header_line)

    # Leer los datos
    df = CSV.read(filename, DataFrame;
                delim=' ',
                ignorerepeated=true,
                comment="#",  # Ahora como String
                header=false,
                skipto=data_start)

    # Asignar nombres de columnas (solo hasta el número de columnas disponibles)
    rename!(df, Symbol.(column_names[1:min(end, ncol(df))]))

    return df
end


function list_age_metal_keys(dirpath::String)
    @show dirpath
    allfiles = filter(!isdir, readdir(dirpath, join=true))
    ages = Float64[]
    metals = Float64[]
    file_ids = Int64[]
    first_file_processed = false
    for i ∈ eachindex(allfiles)
        jldopen(allfiles[i], "r") do file
            # Extraer y convertir las edades
            grupos = keys(file)
            append!(ages, parse.(Float64, replace.(grupos, "age=" => "")))
            append!(file_ids, fill(i, length(grupos)))
            if !first_file_processed
                # Extraer y convertir las metalicidades del primer grupo
                primer_grupo = file[first(grupos)]
                subgrupos = keys(primer_grupo)
                metals = parse.(Float64, replace.(subgrupos, "MH=" => ""))
                first_file_processed = true
            end
        end
    end
    return ages, metals, file_ids, allfiles
end


# Función para encontrar las 4 isócronas más cercanas usando find_interval_uniform
function find_nearest_isochrone(dirpath::String, age::T, metal::R) where {T<:Real,R<:Real}
    ages, metals, file_ids, allfiles = list_age_metal_keys(dirpath)
    age_idx = argmin(abs.(ages.-age))
    metal_idx = argmin(abs.(metals.-metal))
    key_age = @sprintf("age=%0.1f", ages[age_idx])
    key_metal = @sprintf("MH=%+.2f", metals[metal_idx])
    filepath = allfiles[file_ids[age_idx]]
    key = "$key_age/$key_metal"
    df = load(filepath, key)
    return df, key
end


function colorear!(df::DataFrame, mag₁::Symbol, mag₂::Symbol)
    colname = Symbol("color_"*string(mag₁)[1]*string(mag₂)[1])
    df[!, colname] = df[!, mag₁] - df[!,mag₂]
    return nothing
end

function get_evolutionary_phase(label::Int)
    return get(PHASE_MAPPING, label, "Unknown")  # "Unknown" si el label no existe
end

function categorize_phases!(df::DataFrame, family::Symbol)
    if family == :parsec
          phase_map = Dict(
                0 => "0-PMS",
                1 => "1-MS",
                2 => "2-SGB",
                3 => "3-RGB",
                4 => "4-CHeB (i)",
                5 => "5-CHeB (b)",
                6 => "6-CHeB (r)",
                7 => "7-EAGB",
                8 => "8-TPAGB",
                9 => "9-post-AGB"
                )
    end
    get_evolutionary_phase(label::Int) = get(phase_map, label, "Unknown")
    df.phase = get_evolutionary_phase.(df.label)
    return nothing
end

# creo que no la uso mas.
function enumerate_algos(algo::Symbol)::Int64
    if algo == :istreams
        n = 1
    elseif algo == :download
        n= 2
    elseif algo == :ezpadova
        n = 3
    else
        n = 0
    end
    return Int64(n)
end
