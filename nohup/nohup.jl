using Pkg
Pkg.activate(".")
using IntrospectiveStreams

build_isochrone_grid(:parsec,:hsc,(0.1, 14.0, 0.1), (-2.19, 0.5, 0.01), 10:10:140)