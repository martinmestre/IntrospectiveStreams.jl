module IntrospectiveStreams

    using DataFrames, DataFramesMeta
    using PolygonOps
    using StaticArrays
    using FITSIO
    using CSV
    using Interpolations
    using CairoMakie, AlgebraOfGraphics
    using ElectronDisplay
    using PythonCall

    const coord = PythonCall.pynew()
    const u = PythonCall.pynew()
    const at = PythonCall.pynew()
    const galstreams = PythonCall.pynew()
    const galacoord = PythonCall.pynew()
    const pyia = PythonCall.pynew()
    const ezmist = PythonCall.pynew()

    function __init__()
        PythonCall.pycopy!(coord,pyimport("astropy.coordinates"))
        PythonCall.pycopy!(u,pyimport("astropy.units"))
        PythonCall.pycopy!(at,pyimport("astropy.table"))
        PythonCall.pycopy!(galstreams,pyimport("galstreams"))
        PythonCall.pycopy!(galacoord,pyimport("gala.coordinates"))
        PythonCall.pycopy!(pyia,pyimport("pyia"))
        PythonCall.pycopy!(ezmist,pyimport("ezmist"))
    end


    export  correct_extinction_Gaia,
            correct_extinction_Gaia_loop,
            curation!,
            load_stream_track,
            get_isochrone,
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
            clean_xmatch!

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
            plot_track_on_μ_corr_plane_self_frame,
            plot_histog_cmd,
            plot_isochrone_cmd

    export name_files_Gaia,
           name_files_photometry,
           name_files_isochrone,
           name_files_all,
           file_gc

    export basic_pipeline,
           basic_pipeline_loop

    export example

    include("data_methods.jl")
    include("plot_methods.jl")
    include("naming_files.jl")
    include("pipelines.jl")
    include("../examples/example.jl")

end
