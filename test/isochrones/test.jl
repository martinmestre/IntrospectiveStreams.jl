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
#     family, photsys = :mist, "UBVRIplus"
#     for i ∈ 1:n
#         age = 10^f_exp(rand())
#         metal = f_metal(rand())
#         df_iso = download_isochrone(family, photsys, age, metal)
#         @test typeof(df_iso) == DataFrame
#     end
# end

# @testset "DownloadIsochroneEzParsec" begin
#     u_age = 1e9
#     n=4
#     exp_l = 5
#     exp_h = log10((14-1e-12)u_age)
#     metal_l = -2.19999
#     metal_h = 0.5
#     f_exp(x) = (exp_h-exp_l)*x+exp_l
#     f_metal(x) = (metal_h-metal_l)*x+metal_l
#     family, photsys = :parsec, :hsc
#     for i ∈ 1:n
#         age = 10^f_exp(rand())/u_age
#         metal = f_metal(rand())
#         df_iso = download_isochrone(family, photsys, age, metal)
#         @test typeof(df_iso) == DataFrame
#     end
# end

# @testset "DownloadIsochroneEzParsec Grid dispatch" begin
#     u_age = 1e9
#     n=8
#     exp_l = 5
#     exp_h = log10((14-1e-12)u_age)
#     age_l, age_h = 10.0.^(exp_l, exp_h) ./u_age
#     metal_l = -2.19
#     metal_h = 0.5
#     step_age = (age_h-age_l)/n
#     step_metal = (metal_h-metal_l)/n
#     age = (age_l, age_h, step_age)
#     metal = (metal_l, metal_h, step_metal)
#     family, photsys = :parsec, :hsc
#     @show age
#     df_iso = download_isochrone(family, photsys, age, metal)
#     @test typeof(df_iso) == DataFrame
# end



@testset "InterpolateParsecIsochrone -- very unprecise" begin
    u_age = 1e9
    n = 1
    exp_l = 8.5
    exp_h = log10((14-1e-12)u_age)
    metal_l = -2.19
    metal_h = 0.5
    family, photsys = :parsec, :hsc
    f_exp(x) = (exp_h-exp_l)*x+exp_l
    f_metal(x) = (metal_h-metal_l)*x+metal_l
    for i ∈ 1:n
        age = 10^f_exp(rand()) / u_age
        metal = f_metal(rand())
        df_iso = download_isochrone(family, photsys, age, metal)
        df_intp = interpolate_isochrone(family, photsys, age, metal)
        @show nrow(df_iso) nrow(df_intp)
        len = minimum([nrow(df_iso),nrow(df_intp)])
        @show minimum(df_iso.Mini[1:len]) maximum(df_iso.Mini[1:len])
        @show minimum(df_intp.Mini[1:len]) maximum(df_intp.Mini[1:len])
        @test df_intp.Mini[1:len] ≈ df_iso.Mini[1:len] rtol=1.e-1
        @test df_intp.MH[1:len] ≈ df_iso.MH[1:len] rtol=1.e-1
        @test df_intp.logAge[1:len] ≈ df_iso.logAge[1:len] rtol=1.e-1
    end
end

