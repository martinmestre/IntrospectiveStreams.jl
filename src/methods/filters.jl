"""Data curation of Gaia dataset."""
function curation!(df::DataFrame, tol::Vector{<:Number})::Nothing
    df.pmra_error_rel = df.pmra_error./df.pmra
    @subset!(df, (:pmra_error_rel .< tol[1]))
    df.pmdec_error_rel = df.pmdec_error./df.pmdec
    @subset!(df, (:pmdec_error_rel .< tol[2]))
    df.parallax_rel_error = df.parallax_error./df.parallax
    @subset!(df, (:parallax_rel_error .< tol[3]))
    return nothing
end

"""CMD filtering (Mutating)."""
function filter_cmd!(df_stream::DataFrame, df_iso::DataFrame, df_track::DataFrame, photfilter::Symbol, color::Symbol, σ_c::Number; label_lim::Tuple{I,I}=(0,3)) where {I<:Integer}
    phase_mask = label_lim[1] .<= df_iso.label .<= label_lim[2]
    df_isom = df_iso[phase_mask,:]
    df_isom.left = df_isom[!,color] .- σ_c
    df_isom.right = df_isom[!,color] .+ σ_c
    pol_x = vcat(df_isom.left, reverse(df_isom.right), df_isom.left[1])
    mag = df_isom[!, photfilter]
    pol_y = vcat(mag, reverse(mag), mag[1])
    polygon = SVector.(pol_x, pol_y)

    # The line below remains from Gaia df usage.
    # df_stream.distmod = pyconvert(Vector{Float64},coord.Distance(Py(df_stream.D)*u.kpc).distmod.value)
    interpolate_distance!(df_stream, df_track) # Generates :D and :distmod columns at df_stream
    photfilter_abs = Symbol(photfilter,:_abs)
    df_stream[!, photfilter_abs] = df_stream[!, photfilter] - df_stream.distmod
    points = [[df_stream[!, color][i], df_stream[!, photfilter_abs][i]] for i in 1:nrow(df_stream) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stream, identity(inside))
    return nothing
end


"Mask out field globular clusters."
function mask_gc!(df_stream::DataFrame, df_gc::DataFrame)
    for i in 1:nrow(df_gc)
        Δra = df_stream.ra.-df_gc.ra[i]
        Δdec = df_stream.dec.-df_gc.dec[i]
        bool_gc = @. sqrt(Δra^2+Δdec^2) > 0.5
        @subset!(df_stream, identity(bool_gc))
    end
end

"""Filter with stream track on the sky."""
function filter_on_sky!(df_stars::DataFrame, df_track::DataFrame, width::Number)
    up = df_track.ϕ₂.+width
    down =  df_track.ϕ₂.-width
    poly_ϕ₁ = vcat(df_track.ϕ₁, reverse(df_track.ϕ₁), df_track.ϕ₁[1])
    poly_ϕ₂ = vcat(down, reverse(up), down[1])
    polygon = SVector.(poly_ϕ₁, poly_ϕ₂)
    points = [[df_stars.ϕ₁[i], df_stars.ϕ₂[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stars, identity(inside))
    return nothing
end

"""Filter with stream on μ-space."""
function filter_on_μ_plane!(df_stars::DataFrame, df_track::DataFrame, Δμ::Number)
    left = df_track.μ₁cosϕ₂.-Δμ
    right =  df_track.μ₁cosϕ₂.+Δμ
    poly_y = vcat(df_track.μ₂, reverse(df_track.μ₂), df_track.μ₂[1])
    poly_x = vcat(left, reverse(right), left[1])
    polygon = SVector.(poly_x, poly_y)
    points = [[df_stars.μ₁cosϕ₂[i], df_stars.μ₂[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stars, identity(inside))
    return nothing
end

"""Non-mutating filter with stream track in any of its dimensions."""
function filter_along_ϕ₁(df_stars::DataFrame, df_track::DataFrame, S::Symbol, σ::Number)::DataFrame
    if S == :ϕ₂
        q🌠 = df_stars.ϕ₂
        q_track = df_track.ϕ₂
    elseif S == :D
        q🌠 = df_stars.D
        q_track = df_track.D
    elseif S == :μ₁cosϕ₂
        q🌠 = df_stars.μ₁cosϕ₂
        q_track = df_track.μ₁cosϕ₂
    elseif S == :μ₁
        q🌠 = df_stars.μ₁
        q_track = df_track.μ₁
    elseif S == :μ₂
        q🌠 = df_stars.μ₂
        q_track = df_track.μ₂
    elseif S == :Vᵣ
        q🌠 = df_stars.radial_velocity
        q_track = df_track.Vᵣ
    end
    up = q_track .+ σ
    down =  q_track .- σ
    poly_ϕ₁ = vcat(df_track.ϕ₁, reverse(df_track.ϕ₁), df_track.ϕ₁[1])
    poly_q = vcat(down, reverse(up), down[1])
    polygon = SVector.(poly_ϕ₁, poly_q)
    points = [[df_stars.ϕ₁[i], q🌠[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    return @subset(df_stars, identity(inside))
end

"""Mutating filter with stream track in any of its dimensions."""
function filter_along_ϕ₁!(df_stars::DataFrame, df_track::DataFrame, S::Symbol, σ::Number)::Nothing
    if S == :ϕ₂
        q🌠 = df_stars.ϕ₂
        q_track = df_track.ϕ₂
    elseif S == :D
        q🌠 = df_stars.D
        q_track = df_track.D
    elseif S == :μ₁cosϕ₂
        q🌠 = df_stars.μ₁cosϕ₂
        q_track = df_track.μ₁cosϕ₂
    elseif S == :μ₁
        q🌠 = df_stars.μ₁
        q_track = df_track.μ₁
    elseif S == :μ₂
        q🌠 = df_stars.μ₂
        q_track = df_track.μ₂
    elseif S == :Vᵣ
        q🌠 = df_stars.radial_velocity
        q_track = df_track.Vᵣ
    end
    up = q_track .+ σ
    down =  q_track .- σ
    poly_ϕ₁ = vcat(df_track.ϕ₁, reverse(df_track.ϕ₁), df_track.ϕ₁[1])
    poly_q = vcat(down, reverse(up), down[1])
    polygon = SVector.(poly_ϕ₁, poly_q)
    points = [[df_stars.ϕ₁[i], q🌠[i]] for i in 1:nrow(df_stars) ]
    inside = [inpolygon(p, polygon; in=true, on=false, out=false) for p in points]
    @subset!(df_stars, identity(inside))
    return nothing
end

"""Filter with constant value."""
function filter_with_ϕ₂!(df::DataFrame, σ::Number)::Nothing
    @subset!(df, abs.(:ϕ₂) .< σ )
    return nothing
end

"""Filter with a box in PM space."""
function filter_box_on_μ_plane!(df::DataFrame, box::Vector{Vector{Float64}})::Nothing
    @subset!(df, box[1][1] .< :μ₁ .< box[1][2])
    @subset!(df, box[2][1] .< :μ₂ .< box[2][2])
    return nothing
end

"""Filter as in PWB18."""
function filter_from_PWB18!(df::DataFrame)::Nothing
    @subset!(df, abs.(:ϕ₂.-0.5) .< 0.75)
    @subset!(df, -50 .< :ϕ₁ .< -10)
    return nothing
end



"""Clean the cross-match."""
function clean_xmatch!(df::DataFrame)
    # g_0 and g are SPLUS and Gaia mags corrected by extinction.
    unicos = size(unique(df.ID))[1]
    n_multiple = m_multiple = 0
    for i=size(df)[1]:-1:1
        n_repe = count((df.ID.==df.ID[i]))
        if n_repe>1
            n_multiple += 1
            m_multiple += n_repe
            idx = findall(x->x==df.ID[i], df.ID)
            dif_mag = abs.([ df.g_0[idx[1]]-df.g[idx[1]], df.g_0[idx[2]]-df.g[idx[2]] ])
            if dif_mag[1] < dif_mag[2]
                deleteat!(df, [idx[2]])
            else
                deleteat!(df, [idx[1]])
            end
        end
    end
end

