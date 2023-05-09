module IntrospectiveStreams

    using DataFrames, DataFramesMeta
    using PolygonOps
    using StaticArrays
    using PythonCall
    using CairoMakie, AlgebraOfGraphics
    using ElectronDisplay

    @py begin
        import galstreams
        import astropy.coordinates as coord
        import astropy.units as u
        import astropy.table as at
        import ezmist
        import gala.coordinates as galacoord
        import pyia
    end

    export  correct_extinction_Gaia,
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

    export basic_pipeline

    include("data_methods.jl")
    include("plot_methods.jl")
    include("pipelines.jl")
end
