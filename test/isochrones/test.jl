@testset "MetallicityTransformationsUsedInMIST" begin
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

@testset "MetallicityTransformationsUsedInPARSEC" begin
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

@testset "DownloadIsochroneEzMIST" begin
    family, age, metal, filter = :mist, 1.0e7, 0.0, "UBVRIplus"
    df_iso = get_isochrone(family, age, metal, filter)
    @test typeof(df_iso) == DataFrame
end

@testset "DownloadIsochroneEzPARSEC" begin
    family, age, metal, filter = :parsec, 1.0e7, 0.0, "YBC_hsc"
    df_iso = get_isochrone(family, age, metal, filter)
    @test typeof(df_iso) == DataFrame
end
