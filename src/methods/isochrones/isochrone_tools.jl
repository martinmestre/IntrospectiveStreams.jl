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

function find_interval_uniform(A, x)
    # Verificar que x esté dentro del rango de A
    if x < first(A) || x > last(A)
        error("x está fuera del rango de A. x debe estar en [$(first(A)), $(last(A))].")
    end

    # Caso especial: x es igual al último elemento de A
    if x == last(A)
        return length(A) - 1  # Devuelve el índice del penúltimo elemento
    end

    # Calcular el paso (espaciado) de la grilla
    step = (last(A) - first(A)) / (length(A) - 1)

    # Calcular el índice i usando floor
    i = floor(Int, (x - first(A)) / step) + 1

    return i  # Devuelve el índice calculado
end

# Función para encontrar las 4 isócronas más cercanas usando find_interval_uniform
function find_nearest_isochrones(df, target_age, target_metallicity)
    # Mapear la grilla de edad y metalicidad a índices enteros
    unique_ages = sort(unique(df.logAge))
    unique_metallicities = sort(unique(df.MH))

    # Encontrar los índices de los nodos más cercanos usando find_interval_uniform
    age_idx = find_interval_uniform(unique_ages, target_age)
    metallicity_idx = find_interval_uniform(unique_metallicities, target_metallicity)

    # Obtener las cuatro combinaciones de edad y metalicidad más cercanas
    age_nodes = unique_ages[age_idx:age_idx+1]
    metallicity_nodes = unique_metallicities[metallicity_idx:metallicity_idx+1]

    # Crear un DataFrame con las cuatro isócronas más cercanas
    nearest_isochrones = df[
        (df.logAge .∈ [age_nodes]) .& (df.MH .∈ [metallicity_nodes]), :
    ]

    return nearest_isochrones
end






