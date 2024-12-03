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

# @testset "DownloadIsochroneEzMIST" begin
#     n=2
#     exp_l = 6
#     exp_h = 12
#     metal_l = -3.0
#     metal_h = 1.0
#     f_exp(x) = (exp_h-exp_l)*x+exp_l
#     f_metal(x) = (metal_h-metal_l)*x+metal_l
#     for i ∈ 1:2
#         age = 10^f_exp(rand())
#         metal = f_metal(rand())
#         @show age, metal
#         family, age, metal, filter = :mist, age, metal, "UBVRIplus"
#         df_iso = get_isochrone(family, age, metal, filter)
#         @test typeof(df_iso) == DataFrame
#     end
# end

@testset "DownloadIsochroneEzPARSEC" begin
    n=2
    exp_l = 7
    exp_h = 11
    metal_l = -2.99
    metal_h = 1.0
    f_exp(x) = (exp_h-exp_l)*x+exp_l
    f_metal(x) = (metal_h-metal_l)*x+metal_l
    for i ∈ 1:5
        age = 10^f_exp(rand())
        metal = f_metal(rand())
        @show age, metal
        family, age, metal, filter = :parsec, age, metal, "YBC_hsc"
        df_iso = get_isochrone(family, age, metal, filter)
        @test typeof(df_iso) == DataFrame
    end
end
# @testset "DownloadIsochroneEzPARSEC" begin
#     family, age, metal, filter = :parsec, 1.0e7, -3.0, "YBC_hsc"
#     df_iso = get_isochrone(family, age, metal, filter)
#     @test typeof(df_iso) == DataFrame
# end
