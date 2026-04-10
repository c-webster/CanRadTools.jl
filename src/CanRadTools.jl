module CanRadTools

using DelimitedFiles, NCDatasets, ArchGDAL, Plots

include("plot_tools.jl")
include("prep_tools.jl")

export
    check_var_name!,    
    create_gridpts,
    create_gridpts_from_coord,
    create_gridpts_from_file,
    create_gridpts_from_raster,
    plot_static_variable,
    plot_time_varying_gif

end
