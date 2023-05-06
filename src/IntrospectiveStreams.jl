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

    export  extinction_corr,
            curation!,
            load_stream_track,
            get_isochrone,
            compute_in_self_coords!,
            mask_gc!,
            filter_cmd!,
            filter_stream_on_sky!,
            filter_stream_on_μ_plane!,
            filter_along_ϕ₁,
            filter_along_ϕ₁!,
            filter_with_ϕ₂!,
            filter_from_PWB18!,
            filter_box_μ,
            reflex_correct!,
            clean_xmatch!

    export  plot_sky_histo,
            plot_sky_histo_gc,
            plot_sky_histo_selfFrame,
            plot_sky_scatter_selfFrame,
            plot_sky_scatter_μ_arrows_selfFrame,
            plot_sky_scatter_μ_arrows_corr_selfFrame,
            plot_cmd_histo,
            plot_isochrone,
            plot_isochrone_data,
            plot_μ,
            plot_μ_window,
            plot_μ_selfFrame,
            plot_μ_selfFrame_window,
            plot_μ_scatter_selfFrame_window,
            plot_μ_corr_scatter_selfFrame,
            plot_μ_corr_scatter_selfFrame_window,
            plot_μ_corr_histo_selfFrame_window,
            plot_μ_corr_track_selfFrame

    include("data_methods.jl")
    include("plot_methods.jl")

end
