
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


"""Download stellar isochrones"""
function get_isochrone(family::Symbol, age::Number, metal::Number,
                           filter::String; age_scale::String="linear")::DataFrame
    if(family==:mist)
        @assert 5≤log10(age)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
        @assert -4≤metal≤0.5 "Metallicity should satisfy: -4 ≤ FeH ≤ 0.5 (FeH≈MH)."
        println("Note that MIST uses metallicity [FeH] (not Z abundance).")
        df = ezmist.get_one_isochrone(age=age, FeH=metal, v_div_vcrit=0.0,
                    age_scale=age_scale, output_option="photometry",
                    output=filter, Av_value=0.0).to_pandas()|> PyPandasDataFrame |> DataFrame
    elseif(family==:parsec)
        @assert 5≤log10(age)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
        @assert -2.2<metal≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
        println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        df = ezpadova.get_isochrones(age_yr=(age,age,0), MH=(metal,metal,0),
                        model="parsec12s", phot=filter)|> PyPandasDataFrame |> DataFrame
    end
    return df
end
function get_isochrone(family::Symbol, age::NTuple{3,Number}, metal::NTuple{3,Number}, filter::String)::DataFrame
    @assert family==:parsec
    @assert  5 ≤ log10(age[1]) && log10(age[2]) ≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
    @assert -2.2 < metal[1] && metal[2] ≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
    println("Note that Parsec uses metallicity [M/H]=[FeH] (using Z needs to modify get_isochrone function).")
        df = ezpadova.get_isochrones(age_yr=age, MH=metal,model="parsec12s", phot=filter)|> PyPandasDataFrame |> DataFrame
    return df
end


"Interpolate isochrones from the previously downloaded data base (artifacts dir)"
function interpolate_isochrone(family, age, metal, filter)
        @assert family == :parsec "Only Parsec isochrones accepted for the moment"
        if(family==:parsec)
                @assert 9.2≤log10(age)≤10.3 "Age should fulfill: 5 ≤ log10(age) ≤ 10.3."
                @assert -2.2<metal≤0.5 "Metallicity should satisfy: -2.2 < FeH ≤ 0.5 (FeH≈MH)."
                if(filter=="YBC_hsc")
                        file_artif = "artifacts/isochrones/parsec/SubaruHSC/family_MH_-2.2_0.5_logAge_9.2_10.3.dat"
                        df_artif = DataFrame(CSV.File(file_artif, delim=" ", ignorerepeated=true, comment="#"))
                        points = [df_artif.MH df_artif.logAge df_artif.Mini]'
                        samples = df_artif.gmag
                        itp = ScatteredInterpolation.interpolate(Multiquadratic(), points, samples)
                        interpolated = evaluate(itp, [0.0, 10.0, 1.0])
                        @show itp interpolated
                        return
                end
        end
end