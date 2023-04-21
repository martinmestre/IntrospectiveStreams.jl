module IntrospectiveStreams

    using DataFrames, DataFramesMeta
    using PolygonOps
    using StaticArrays
    using PythonCall

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
            compute_in_selfCoords!,
            mask_gc!,
            filter_cmd!,
            filter_stream_on_sky!,
            filter_stream_μ_space!,
            filter_with_track,
            filter_with_track!,
            filter_with_ϕ₂!,
            filter_PWB18!,
            filter_box_μ,
            reflex_correct!


    include("DataMethods.jl")

end
