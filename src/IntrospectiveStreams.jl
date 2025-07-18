module IntrospectiveStreams

    using Reexport
    @reexport using DataFrames, DataFramesMeta
    using StaticArrays
    using FITSIO
    @reexport using CSV, JLD2
    using Interpolations
    @reexport using Colors
    @reexport using CairoMakie, AlgebraOfGraphics
    # using ElectronDisplay
    @reexport using PythonCall
    @reexport using Printf
    @reexport using LaTeXStrings
    @reexport using Showoff
    @reexport using StatsBase
    using PolygonOps

    const aog = AlgebraOfGraphics

    const coord = PythonCall.pynew()
    const u = PythonCall.pynew()
    const at = PythonCall.pynew()
    const galstreams = PythonCall.pynew()
    const galacoord = PythonCall.pynew()
    const pyia = PythonCall.pynew()
    const ezmist = PythonCall.pynew()
    const ezpadova = PythonCall.pynew()
    const os = PythonCall.pynew()

    const PKGDIR = normpath(joinpath(@__DIR__, ".."))

    function __init__()
        PythonCall.pycopy!(coord,pyimport("astropy.coordinates"))
        PythonCall.pycopy!(u,pyimport("astropy.units"))
        PythonCall.pycopy!(at,pyimport("astropy.table"))
        PythonCall.pycopy!(galstreams,pyimport("galstreams"))
        PythonCall.pycopy!(galacoord,pyimport("gala.coordinates"))
        PythonCall.pycopy!(pyia,pyimport("pyia"))
        PythonCall.pycopy!(ezmist,pyimport("ezmist"))
        PythonCall.pycopy!(ezpadova,pyimport("ezpadova"))
        PythonCall.pycopy!(os,pyimport("os"))
        os.environ["SSL_CERT_FILE"] = "/etc/ssl/certs/ca-certificates.crt"
    end


    export  u,
            coord,
            at,
            galacoord,
            pyia,
            ezmist,
            ezpadova


    export  get_dataframes,
            correct_extinction_Gaia,
            correct_extinction_Gaia_loop,
            curation!,
            load_stream_track,
            compute_in_self_coords!,
            mask_gc!,
            filter_cmd!,
            filter_on_sky!,
            filter_on_μ_plane!,
            filter_along_ϕ₁,
            filter_along_ϕ₁!,
            filter_with_ϕ₂!,
            filter_PWB18!,
            filter_box_on_μ_plane!,
            reflex_correct!,
            reflex_correct_steps!,
            clean_xmatch!,
            dropinfinite!,
            dropinfinite,
            interpolate_distance!

    export  z2feh_mist, z2feh_parsec,
            feh2z_mist, feh2z_parsec,
            list_age_metal_keys,
            remove_age_group,
            colorear!,
            categorize_phases!,
            download_isochrone,
            build_isochrone_grid,
            find_nearest_isochrone,
            interpolate_isochrone

    export  plot_histog_on_sky,
            plot_histog_on_sky_with_gc,
            plot_histog_on_sky_self_frame,
            plot_scatter_on_sky_self_frame,
            plot_scatter_on_sky_μ_arrows_self_frame,
            plot_scatter_on_sky_μ_corr_arrows_self_frame,
            plot_histog_on_μ_plane,
            plot_histog_on_μ_plane_self_frame,
            plot_scatter_on_μ_plane_self_frame,
            plot_scatter_on_μ_corr_plane_self_frame,
            plot_track_on_μ_rc_plane_self_frame,
            plot_histog_cmd,
            plot_cmd,
            plot_mags_density,
            plot_mags_histogram,
            plot_scatter_on_sky

    export  wongcolors,
            wongcolors_ext,
            paulcolors

    export  print_counts_per_mag_bin
    export  rename_mag!,
            rename_magerr!


    export name_files_Gaia,
           name_files_photometry,
           name_files_isochrone,
           name_files_all,
           file_gc

    export basic_pipeline,
           basic_pipeline_loop

    export example_pipeline,
           example_interpolate_isochrone

    export mode


    include("methods/extinctions.jl")
    include("methods/transformations.jl")
    include("methods/filters.jl")
    include("methods/get_data.jl")
    include("methods/isochrones/isochrone_tools.jl")
    include("methods/isochrones/isochrones.jl")
    include("methods/naming_files.jl")
    include("methods/missing.jl")
    include("methods/colors.jl")
    include("methods/plots.jl")
    include("methods/printing.jl")
    include("methods/renaming.jl")
    include("methods/statistics.jl")
    include("pipelines.jl")
    include("../examples/example_pipeline.jl")
    include("../examples/example_isochrone.jl")
end
