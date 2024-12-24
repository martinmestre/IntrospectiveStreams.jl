@testset "MetallicityTransformationsUsedInMist" begin
    n = 50
    for i ∈ 1:n
        feh = -4 + rand()*(1-(-4))
        z = feh2z_mist(feh)
        feh₂ = z2feh_mist(z)
        @test feh₂ ≈ feh rtol=5.0e-7
        z₂ = feh2z_mist(feh)
        @test z₂ ≈ z rtol=5.0e-7
    end
end

@testset "MetallicityTransformationsUsedInParsec" begin
    n = 50
    for i ∈ 1:n
        feh = -4 + rand()*(1 - (-4))
        z = feh2z_parsec(feh)
        feh₂ = z2feh_parsec(z)
        @test feh₂ ≈ feh rtol=5.0e-5
        z₂ = feh2z_parsec(feh)
        @test z₂ ≈ z rtol=5.0e-7
    end
end

@testset "DownloadIsochroneEzMist" begin
    n=2
    exp_l = 9
    exp_h = 10.3
    metal_l = -4
    metal_h = 0.5
    f_exp(x) = (exp_h-exp_l)*x+exp_l
    f_metal(x) = (metal_h-metal_l)*x+metal_l
    for i ∈ 1:n
        age = 10^f_exp(rand())
        metal = f_metal(rand())
        @show age, metal
        family, age, metal, filter = :mist, age, metal, "UBVRIplus"
        df_iso = get_isochrone(family, age, metal, filter)
        @test typeof(df_iso) == DataFrame
    end
end

@testset "DownloadIsochroneEzParsec" begin
    n=2
    exp_l = 9
    exp_h = 10.3
    metal_l = -2.19999
    metal_h = 0.5
    f_exp(x) = (exp_h-exp_l)*x+exp_l
    f_metal(x) = (metal_h-metal_l)*x+metal_l
    for i ∈ 1:n
        age = 10^f_exp(rand())
        metal = f_metal(rand())
        @show age, metal
        family, age, metal, filter = :parsec, age, metal, "YBC_hsc"
        df_iso = get_isochrone(family, age, metal, filter)
        @test typeof(df_iso) == DataFrame
    end
end

@testset "InterpolateParsecIsochone" begin
    n=2
    exp_l = 9
    exp_h = 10.3
    metal_l = -2.19999
    metal_h = 0.5
    f_exp(x) = (exp_h-exp_l)*x+exp_l
    f_metal(x) = (metal_h-metal_l)*x+metal_l
    for i ∈ 1:n
        age = 10^f_exp(rand())
        metal = f_metal(rand())
        @show age, metal
        family, age, metal, filter = :parsec, age, metal, "YBC_hsc"
        df_download = get_isochrone(family, age, metal, filter)
        df_interpol = interpolate_isochrone(family, age, metal, filter)
        @test df_download ≈ df_interpol rtol=5.0e-7
    end
end

