# @testset "SSL certificates used by EzPadova" begin
#     python_exe = ENV["JULIA_PYTHONCALL_EXE"]
#     env_name = split(python_exe, "/")[end-2]  # El nombre del entorno está en la penúltima posición
#     ssl_cert_file = joinpath("/home/mmestre/.conda/envs", env_name, "ssl/cert.pem")
#     ssl_cert_dir = joinpath("/home/mmestre/.conda/envs", env_name, "ssl/certs")
#     @test ENV["SSL_CERT_FILE"] == ssl_cert_file
#     @test ENV["SSL_CERT_DIR"] == ssl_cert_dir
# end

# @testset "MetallicityTransformationsUsedInMist" begin
#     n = 50
#     for i ∈ 1:n
#         feh = -4 + rand()*(1-(-4))
#         z = feh2z_mist(feh)
#         feh₂ = z2feh_mist(z)
#         @test feh₂ ≈ feh rtol=5.0e-7
#         z₂ = feh2z_mist(feh)
#         @test z₂ ≈ z rtol=5.0e-7
#     end
# end

# @testset "MetallicityTransformationsUsedInParsec" begin
#     n = 50
#     for i ∈ 1:n
#         feh = -4 + rand()*(1 - (-4))
#         z = feh2z_parsec(feh)
#         feh₂ = z2feh_parsec(z)
#         @test feh₂ ≈ feh rtol=5.0e-5
#         z₂ = feh2z_parsec(feh)
#         @test z₂ ≈ z rtol=5.0e-7
#     end
# end

# @testset "DownloadIsochroneEzMist" begin
#     n=4
#     exp_l = 5
#     exp_h = 10.3
#     metal_l = -4
#     metal_h = 0.5
#     f_exp(x) = (exp_h-exp_l)*x+exp_l
#     f_metal(x) = (metal_h-metal_l)*x+metal_l
#     family, filter = :mist, "UBVRIplus"
#     for i ∈ 1:n
#         age = 10^f_exp(rand())
#         metal = f_metal(rand())
#         df_iso = get_isochrone(family, filter, age, metal)
#         @test typeof(df_iso) == DataFrame
#     end
# end

@testset "DownloadIsochroneEzParsec" begin
    n=4
    exp_l = 5
    exp_h = 10.3
    metal_l = -2.19999
    metal_h = 0.5
    f_exp(x) = (exp_h-exp_l)*x+exp_l
    f_metal(x) = (metal_h-metal_l)*x+metal_l
    family, filter = :parsec, "hsc"
    for i ∈ 1:n
        age = 10^f_exp(rand())
        metal = f_metal(rand())
        df_iso = get_isochrone(family, filter, age, metal)
        @test typeof(df_iso) == DataFrame
    end
end

# @testset "DownloadIsochroneEzParsec Grid dispatch" begin
#     n=8
#     exp_l = 5
#     exp_h = 10.3
#     age_l, age_h = 10.0.^(exp_l, exp_h)
#     metal_l = -2.19999
#     metal_h = 0.5
#     step_age = (age_h-age_l)/n
#     step_metal = (metal_h-metal_l)/n
#     age = (age_l, age_h, step_age)
#     metal = (metal_l, metal_h, step_metal)
#     family, filter = :parsec, "hsc"
#     df_iso = get_isochrone(family, filter, age, metal)
#     @test typeof(df_iso) == DataFrame
# end



@testset "InterpolateParsecIsochrone" begin
    n = 1
    exp_l = 5
    exp_h = 10.3
    metal_l = -2.19999
    metal_h = 0.5
    family, filter = :parsec, "hsc"
    for i ∈ 1:n
        age = 10^f_exp(rand())
        metal = f_metal(rand())
        df_iso = get_isochrone(family, filter, age, metal)
        df_intp = interpolate_isochrone(family, filter, age, metal)
        @test df_intp.Mini ≈ df_iso.Min rtol=1.e-5
        @test df_intp.MH ≈ df_iso.MH rtol=1.e-5
        @test df_intp.logAge ≈ df_iso.logAge rtol=1.e-5
        @test df_intp ≈ df_iso rtol=1.e-5 # This includes the above comparison per field and more.
    end
end

