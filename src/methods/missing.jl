function columnas_con_missing(df::DataFrame)
    # Encuentra columnas con al menos un missing
    cols_con_missing = String[]

    for col in names(df)
        if any(ismissing.(df[!, col]))
            push!(cols_con_missing, col)
        end
    end

    # Mostrar resultados
    if isempty(cols_con_missing)
        println("✅ El DataFrame no contiene columnas con valores missing")
    else
        println("⚠️ Columnas con valores missing (", length(cols_con_missing), "):")
        for col in cols_con_missing
            n_missing = sum(ismissing.(df[!, col]))
            pct_missing = round(n_missing / nrow(df) * 100, digits=2)
            println("- ", col, ": ", n_missing, " missing (", pct_missing, "%)")
        end
    end

    return cols_con_missing
end

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

function dropinfinite!(df::DataFrame, cols::Vector{Symbol}=[:magnitude])
    # Combine column checks with OR using proper broadcasting
    inf_mask = reduce((a,b) -> a .| b, [isinf.(df[!, col]) for col in cols])
    delete!(df, inf_mask)
    return nothing
end

function dropinfinite(df::DataFrame, cols::Vector{Symbol}=[:magnitude])
    # Combine column checks with OR using proper broadcasting
    inf_mask = reduce((a,b) -> a .| b, [isinf.(df[!, col]) for col in cols])
    df_c = df[Not(inf_mask),:]
    return df_c
end
