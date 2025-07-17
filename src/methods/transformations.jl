"""Compute stream stars' self coordinates using Galstreams'track and add to dataframe."""
function compute_in_self_coords!(df::DataFrame, frame::Py)::Nothing
    sky_coords = coord.SkyCoord(ra=Py(df.ra)*u.deg, dec=Py(df.dec)*u.deg, pm_ra_cosdec=Py(df.pmra)*u.mas/u.yr, pm_dec=Py(df.pmdec)*u.mas/u.yr, frame="icrs")
    self_coords = sky_coords.transform_to(frame)
    df.ϕ₁ = pyconvert(Vector{Float64}, self_coords.phi1.deg)
    df.ϕ₂ = pyconvert(Vector{Float64}, self_coords.phi2.deg)
    df.μ₁cosϕ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi1_cosphi2.value)
    df.μ₁ = @. df.μ₁cosϕ₂/cos(df.ϕ₂*π/180.0)
    df.μ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi2.value)
    df.D_Π = 1.0./df.parallax
    return nothing
end

"""Compute stream stars' self coordinates using gala's KOPOSOV frame and add to dataframe."""
function compute_in_self_coords!(df::DataFrame)::Py
    sky_coords = coord.SkyCoord(ra=Py(df.ra)*u.deg, dec=Py(df.dec)*u.deg, pm_ra_cosdec=Py(df.pmra)*u.mas/u.yr, pm_dec=Py(df.pmdec)*u.mas/u.yr, frame="icrs")
    kop_frame = galacoord.GD1
    self_coords = sky_coords.transform_to(kop_frame)
    df.ϕ₁ = pyconvert(Vector{Float64}, self_coords.phi1.deg)
    df.ϕ₂ = pyconvert(Vector{Float64}, self_coords.phi2.deg)
    df.μ₁cosϕ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi1_cosphi2.value)
    df.μ₁ = @. df.μ₁cosϕ₂/cos(df.ϕ₂*π/180.0)
    df.μ₂ = pyconvert(Vector{Float64}, self_coords.pm_phi2.value)
    df.D_Π = 1.0./df.parallax
    return kop_frame
end

"""Reflex Correction step by step, using the self_frame at SkyCoord object,
without calling Gala."""
function reflex_correct_steps!(df::DataFrame, frame::Py)::Nothing
    len = length(df.ϕ₁)
    sky_coords = coord.SkyCoord(phi1=Py(df.ϕ₁)*u.deg, phi2=Py(df.ϕ₂)*u.deg, pm_phi1_cosphi2=Py(df.μ₁cosϕ₂)*u.mas/u.yr, pm_phi2=Py(df.μ₂)*u.mas/u.yr, distance=Py(df.D)*u.kpc, radial_velocity=Py(fill(0.,len))*u.km/u.s, frame=frame)
    v☼ = coord.CartesianDifferential(Py([11.1, 220.0+12.24, 7.25])*u.km/u.s)
    r☼ = 8.122*u.kpc
    gc_frame = coord.Galactocentric(galcen_distance=r☼, galcen_v_sun=v☼, z_sun=0*u.pc)
    observed = sky_coords.transform_to(gc_frame)
    rep = observed.cartesian.without_differentials()
    rep = rep.with_differentials(observed.cartesian.differentials["s"] + v☼)
    sky_coords_rc = gc_frame.realize_frame(rep).transform_to(frame)
    # sky_coords_rc = galacoord.reflex_correct(sky_coords, gc_frame)
    df.μ₁cosϕ₂_rc = pyconvert(Vector{Float64}, sky_coords_rc.pm_phi1_cosphi2.value)
    df.μ₁_rc = @. df.μ₁cosϕ₂_rc/cos(df.ϕ₂*π/180.0)
    df.μ₂_rc = pyconvert(Vector{Float64}, sky_coords_rc.pm_phi2.value)
    return nothing
end

"""Reflex Correction calling Gala, using the ICRS frame for SkyCoord."""
function reflex_correct!(df::DataFrame, frame::Py)::Nothing
    len = length(df.ϕ₁)
    sky_coords_icrs = coord.SkyCoord(ra=Py(df.ra)*u.deg, dec=Py(df.dec)*u.deg, pm_ra_cosdec=Py(df.pmra)*u.mas/u.yr, pm_dec=Py(df.pmdec)*u.mas/u.yr, distance=Py(df.D)*u.kpc, radial_velocity=Py(fill(0.,len))*u.km/u.s, frame="icrs")
    v☼ = coord.CartesianDifferential(Py([11.1, 220.0+12.24, 7.25])*u.km/u.s)
    r☼ = 8.122*u.kpc
    gc_frame = coord.Galactocentric(galcen_distance=r☼, galcen_v_sun=v☼, z_sun=0*u.pc)
    sky_coords_icrs_rc = galacoord.reflex_correct(sky_coords_icrs, gc_frame)
    sky_coords_rc = sky_coords_icrs_rc.transform_to(frame)
    df.μ₁cosϕ₂_rc = pyconvert(Vector{Float64}, sky_coords_rc.pm_phi1_cosphi2.value)
    df.μ₁_rc = @. df.μ₁cosϕ₂_rc/cos(df.ϕ₂*π/180.0)
    df.μ₂_rc = pyconvert(Vector{Float64}, sky_coords_rc.pm_phi2.value)
    return nothing
end

"""
    interpolate_distance_from_track!(df_stream::T, df_track::T)::Nothing where {T::DataFrame}
This function receives a stream df that has been already reduced to be within the "ra" domain of the track df.
"""
function interpolate_distance!(df_stream::DataFrame, df_track::DataFrame, var::Symbol=:ra)
    # itp = cubic_spline_interpolation(df_track[!,var], df_track.D)
    itp = linear_interpolation((df_track[!,var],), df_track.D)
    df_stream.D = itp.(df_stream[!,var])
    df_stream.distmod = 5log10.(df_stream.D) .+ 10
    return nothing
end