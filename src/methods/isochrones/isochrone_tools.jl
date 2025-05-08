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



function list_age_metal_keys(filepath)
    jldopen(filepath, "r") do file
        # Extraer y convertir las edades
        grupos = keys(file)
        edades = parse.(Float64, replace.(grupos, "age=" => ""))

        # Extraer y convertir las metalicidades del primer grupo
        primer_grupo = file[first(grupos)]
        subgrupos = keys(primer_grupo)
        metalicidades = parse.(Float64, replace.(subgrupos, "MH=" => ""))

        return edades, metalicidades
    end
end



# Función para encontrar las 4 isócronas más cercanas usando find_interval_uniform
function find_nearest_isochrone(file_artif::String, age::T, metal::R) where {T<:Real,R<:Real}
    age_keys, metal_keys = list_age_metal_keys(file_artif)
    age_idx = argmin(abs.(age_keys.-age))
    metal_idx = argmin(abs.(metal_keys.-metal))
    key_age = @sprintf("age=%0.1f", age_keys[age_idx])
    key_metal = @sprintf("MH=%+.2f", metal_keys[metal_idx])
    key = "$key_age/$key_metal"
    return load(file_artif, key), key
end



