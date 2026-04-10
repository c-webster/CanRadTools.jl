"""
    plot_static_variable(infile::String,vf_name::String="vf_planar",outfiletype::String="png",timestamp::String=nothing,outdir::String=nothing)

Creates a geotiff or .png file of the specifed variable. 
If the variable is time varying (tvt or swr), then the timestamp must be supplied

# Arguments
- `infile`: Path to the input NetCDF file created using CanRad.create_grid_from_file()
- `var_name`: Name of the variable in the NetCDF file (default is "vf_planar"). If the specified variable name is not found, the function will search for variables containing the specified name and use the first match.
- `outfiletype`: Type of the output file, either "png" for a heatmap plot or "tif" for a georeferenced TIFF file (default is "png").
- `timestamp`: Optional timestamp to specify which time slice to plot for time-varying variables (as string in datetime format e.g. "2020-01-01T00:00:00").
- `outdir`: Directory where the output file will be saved. If not specified, the output file will be saved in the same directory as the input file.

# Returns
- A georeferenced TIFF file or a heatmap PNG file of the specified variable.

# Requirements
- input netcdf file must be created using CanRad.create_grid_from_file() to ensure the correct structure and metadata for georeferencing.

"""
function plot_static_variable(infile::String;var_name::String="vf_planar",outfiletype::String="png",timestamp::String="nothing",outdir::String="nothing",epsgnum::Number=-1)

    # load the data
    ds = NCDataset(infile)

    (outdir === "nothing") && (outdir = dirname(infile))
    !isdir(outdir) && mkdir(outdir)
    if timestamp == "nothing"
        outname = replace(basename(infile), ".nc" => "_$(var_name).$outfiletype")
    else
        outname = replace(basename(infile), ".nc" => "_$(var_name)_$timestamp.$outfiletype")
    end

    var_name = check_var_name!(var_name,keys(ds))

    if ndims(ds[var_name]) == 2
        dat = Float64.(coalesce.(ds[var_name][:,:], NaN))
    elseif ndims(ds[var_name]) == 3
        timedx = findall(string.(ds["datetime"][:]) .== timestamp)[1]
        dat = Float64.(coalesce.(ds[var_name][timedx,:,:], NaN))
    end

    if any(contains(var_name, v) for v in ["tvt","svf"])
        dat = dat ./ 100
        clim = (0,1)
    else
        clim = (0, ceil(maximum(dat) / 50) * 50)
    end

    easting = ds["easting"][:]
    northing = ds["northing"][:]
    cellsize = diff(sort(unique(easting)))[1]   

    # create the plot
    if outfiletype == "png"

        # sort y in ascending order and reorder dat accordingly
        xlims = (minimum(easting)-cellsize/2, maximum(easting)+cellsize/2)
        ylims = (minimum(northing)-cellsize/2, maximum(northing)+cellsize/2)
        p = heatmap(easting,reverse(northing),reverse(dat,dims=1),aspect_ratio=1,title=var_name,clim=clim,xlim=xlims,ylim=ylims)
        savefig(p,joinpath(outdir,outname))

    elseif outfiletype == "tif"

        if epsgnum == -1
            try
                epsg = ds["easting"].attrib["epsg"]
            catch
               @error "EPSG code not found in dataset attributes. Please provide an EPSG code as a function argument." 
            end
        else
            epsg = epsgnum
        end
        crs = ArchGDAL.toWKT(ArchGDAL.importEPSG(epsg))
        cellsize = diff(sort(unique(easting)))[1]
        nrows,ncols = size(nomissing(dat))
        xllcorner   = minimum(easting)-(cellsize/2)
        yllcorner   = minimum(northing)-(cellsize/2)
        gt = Float64.([xllcorner,cellsize,0,yllcorner+(nrows*cellsize),0,-cellsize])
        outdat = Float64.(transpose(dat))
        outdat[isnan.(outdat)] .= -1
        outfname = joinpath(outdir,outname)
        ArchGDAL.create(
            outfname,driver = ArchGDAL.getdriver("GTiff"),
            width=Int(ncols),height=Int.(nrows),nbands=1,dtype=Float64
        )do dataset
            ArchGDAL.write!(dataset,outdat, 1)
            # band = ArchGDAL.getband(dataset, 1)
            # ArchGDAL.setnodatavalue!(band, -9999.0)
            ArchGDAL.setgeotransform!(dataset, gt)
            ArchGDAL.setproj!(dataset, crs)
        end

    end

    close(ds)

    
end

"""
    check_var_name!(var_name::String,list_vars::Vector{String})

Used in plot_static_variable and plot_time_varying_gif to check if the specified variable name exists in the dataset. If not, it searches for variables containing the specified name and uses the first match.
"""
function check_var_name!(var_name::String,list_vars::Vector{String})

    if !(var_name in list_vars)
        matching_vars = [var for var in list_vars if contains(var, var_name)]
        if !isempty(matching_vars)
            println("Found matching variable(s): ", matching_vars)
            println("Using variable: ", matching_vars[1])
            var_name = matching_vars[1]
        else
            error("Variable $var_name not found in dataset")
        end
    end

    return var_name

end

"""
    plot_time_varying_gif(infile::String,var_name::String,t1::String,t2::String,outdir::String=nothing)

creates a gif of the specified variable between the specified time range.

# Arguments
- `infile`: Path to the input NetCDF file created using CanRad.create_grid_from_file()
- `var_name`: Name of the variable in the NetCDF file (e.g. "tvt"). If the specified variable name is not found, the function will search for variables containing the specified name and use the first match.
- `t1`: Start time of the range to plot (as string in datetime format e.g. "2020-01-01T00:00:00")   
- `t2`: End time of the range to plot (as string in datetime format e.g. "2020-01-01T00:00:00")
- `fps`: Frames per second for the output gif (default is 1)
- `outdir`: Directory where the output file will be saved. If not specified, the output file will be saved in the same directory as the input file.

"""
function plot_time_varying_gif(infile::String,var_name::String,t1::String,t2::String,fps::Int8=1,outdir::String="nothing")

    ds = NCDataset(infile)

    (outdir === "nothing") && (outdir = dirname(infile))
    !isdir(outdir) && mkdir(outdir)
    outname = replace(basename(infile), ".nc" => "_$(var_name)_$(t1)_$(t2).gif")

    tdx1 = findall(string.(ds["datetime"][:]) .== t1)[1]
    tdx2 = findall(string.(ds["datetime"][:]) .== t2)[1]

    var_name = check_var_name!(var_name,keys(ds))

    dat = ds[var_name][tdx1:tdx2,:,:]
    time_axis = string.(ds["datetime"][tdx1:tdx2])

    if any(contains(var_name, v) for v in ["tvt"])
        dat = dat ./ 100
        clim = (0,1)
    else # variable is shortwave radiation
        clim = (0, ceil(maximum(dat) / 50) * 50)
    end

    easting = ds["easting"][:]
    northing = ds["northing"][:]

    cellsize = diff(sort(unique(easting)))[1]
    xlims = (minimum(easting)-cellsize/2, maximum(easting)+cellsize/2)
    ylims = (minimum(northing)-cellsize/2, maximum(northing)+cellsize/2)

    # plot the gif
    anim = @animate for i in axes(dat, 1)
        date_str = split(time_axis[i],"T")
        title_str = "$(date_str[1])\n$(date_str[2])"
        heatmap(easting, reverse(northing), reverse(dat[i,:,:], dims=1), 
                aspect_ratio=1, title=title_str, clim=clim, 
                xlim=xlims, ylim=ylims)
    end 
    gif(anim, joinpath(outdir, outname),fps=fps)

    close(ds)

end