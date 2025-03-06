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
#     n=2
#     exp_l = 9
#     exp_h = 10.3
#     metal_l = -4
#     metal_h = 0.5
#     f_exp(x) = (exp_h-exp_l)*x+exp_l
#     f_metal(x) = (metal_h-metal_l)*x+metal_l
#     for i ∈ 1:n
#         age = 10^f_exp(rand())
#         metal = f_metal(rand())
#         @show age, metal
#         family, age, metal, filter = :mist, age, metal, "UBVRIplus"
#         df_iso = get_isochrone(family, age, metal, filter)
#         @test typeof(df_iso) == DataFrame
#     end
# end

# @testset "DownloadIsochroneEzParsec" begin
#     n=4
#     exp_l = 9
#     exp_h = 10.3
#     metal_l = -2.19999
#     metal_h = 0.5
#     f_exp(x) = (exp_h-exp_l)*x+exp_l
#     f_metal(x) = (metal_h-metal_l)*x+metal_l
#     for i ∈ 1:n
#         age = 10^f_exp(rand())
#         metal = f_metal(rand())
#         @show age, metal
#         family, age, metal, filter = :parsec, age, metal, "YBC_hsc"
#         df_iso = get_isochrone(family, age, metal, filter)
#         @show typeof(df_iso) == DataFrame
#     end
# end

@testset "DownloadIsochroneEzParsec Grid dispatch" begin
    n=4
    exp_l = 10
    exp_h = 10.2
    metal_l = -2.
    metal_h = 0.3
    age = (10^exp_l, 10^exp_h, n)
    metal = (metal_l, metal_h, n)
    @show age  metal typeof(age) typeof(metal)
    family, age, metal, filter = :parsec, age, metal, "YBC_hsc"
    df_iso = get_isochrone(family, age, metal, filter)
    @show typeof(df_iso) == DataFrame
end







# @testset "InterpolateParsecIsochone" begin
#     n=2
#     exp_l = 9
#     exp_h = 10.3
#     metal_l = -2.19999
#     metal_h = 0.5
#     f_exp(x) = (exp_h-exp_l)*x+exp_l
#     f_metal(x) = (metal_h-metal_l)*x+metal_l
#     for i ∈ 1:n
#         age = 10^f_exp(rand())
#         metal = f_metal(rand())
#         @show age, metal
#         family, age, metal, filter = :parsec, age, metal, "YBC_hsc"
#         df_download = get_isochrone(family, age, metal, filter)
#         df_interpol = interpolate_isochrone(family, age, metal, filter)
#         @test df_download ≈ df_interpol rtol=5.0e-7
#     end
# end

