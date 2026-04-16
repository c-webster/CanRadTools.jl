module CanRadTools

using DelimitedFiles, NCDatasets, ArchGDAL, Plots

include("plot_tools.jl")
include("prep_tools.jl")

export
    check_var_name!,
    create_grid_from_file,
    create_grid_from_multiple_files,
    create_gridpts,
    create_gridpts_from_coord,
    create_gridpts_from_file,
    create_gridpts_from_raster,
    environments_flag_local,
    plot_static_variable,
    plot_time_varying_gif,
    write_grid_to_netcdf

end
