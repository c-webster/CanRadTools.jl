"""
    create_gridpts(xmin::Float64, xmax::Float64, ymin::Float64, ymax::Float64, spacing::Float64, taskID::String="task", outdir::String=".")

Create a regular grid of points within a specified rectangular region and save them to a file.

# Arguments
- `xmin::Float64`: Minimum x-coordinate of the grid region
- `xmax::Float64`: Maximum x-coordinate of the grid region
- `ymin::Float64`: Minimum y-coordinate of the grid region
- `ymax::Float64`: Maximum y-coordinate of the grid region
- `spacing::Number`: Distance between grid points in both x and y directions
- `buffer::Number`: additional buffer to add (positive) or subtract (negative) from the limits (default: 0)
- `taskID::String`: Identifier for the output file (default: "task")
- `outdir::String`: Output directory path (default: current directory ".")

# Returns
Nothing. The function writes grid points to a delimited text file.
Output values are in the same coordinate system as the input limits.
Output values are whole numbers, so any decimal values in the limits are rounded down for minimums and rounded up for maximums.

# Requirements
DelimitedFiles.jl package is required to read the input file and write the output file.

# Details
Generates a regular grid of points by creating linearly spaced coordinates in x and y dimensions
with the specified spacing. All combinations of x and y coordinates are created, resulting in a
rectangular grid pattern. The grid points are saved to a file named `{taskID}_gridpts.txt`.

# Example

"""
function create_gridpts(xmin::Float64,xmax::Float64,ymin::Float64,ymax::Float64,spacing::Number,buffer::Number=0.0,taskID::String="task",outdir::String=".")

    xpts = collect(xmin:spacing:xmax)
    ypts = collect(ymin:spacing:ymax)
    gridpts = [(x,y) for x in xpts for y in ypts]

    outpath = joinpath(outdir,"$(taskID)_gridpts.txt")

    writedlm(outpath,gridpts)

end

"""
    create_gridpts_from_coord(xvalue::Float64, yvalue::Float64, spacing::Float64, buffer::Float64, taskID::String="task", outdir::String=".")

Create a regular grid of points centered around a specified coordinate and save them to a file.

# Arguments
- `xvalue::Float64`: x-coordinate of the center point
- `yvalue::Float64`: y-coordinate of the center point
- `spacing::Float64`: Distance between grid points in both x and y directions
- `buffer::Float64`: Distance from the center point to the outer edge of the grid
- `taskID::String`: Identifier for the output file (default: "task")
- `outdir::String`: Output directory path (default: current directory ".")
- `delim::String`: Delimiter for the output file (default: tab). Options are "tab" or "comma"

# Returns
Nothing. The function writes grid points to a delimited text file.
Output values are in the same coordinate system as the input limits.

# Requirements
DelimitedFiles.jl package is required to read the input file and write the output file.

"""
function create_gridpts_from_coord(xvalue::Float64,yvalue::Float64,spacing::Float64,buffer::Float64,taskID::String="task",outdir::String=".",delim::String="tab")

    xpts = collect(xvalue-buffer:spacing:xvalue+buffer)
    ypts = collect(yvalue-buffer:spacing:yvalue+buffer)
    gridpts = [(x,y) for x in xpts for y in ypts]

    outpath = joinpath(outdir,"$(taskID)_gridpts.txt")

    if delim == "tab"
        writedlm(outpath,gridpts,'\t')
    elseif delim == "comma"
        writedlm(outpath,gridpts,',')
    end

end


"""
    create_gridpts_from_file(infile::String, spacing::Float64, taskID::String="task", outdir::String=".", buffer::Float64=0)

Create a regular grid of points based on limits defined in an input file and save them to a file.

# Arguments
- `infile::String`: Path to the input file containing either point coordinates or limits
- `spacing::Float64`: Distance between grid points in both x and y directions
- `taskID::String`: Identifier for the output file (default: "task")
- `outdir::String`: Output directory path (default: current directory ".")
- `buffer::Float64`: Optional buffer to add (positive) or subtract (negative) from the limits (default: 0)

# Returns
Nothing. The function writes grid points to a delimited text file.
Output values are in the same coordinate system as the input limits.
Output values are whole numbers, so any decimal values in the limits are rounded down for minimums and rounded up for maximums.

# Requirements
DelimitedFiles.jl package is required to read the input file and write the output file.

# Details
The function reads an input file that can be in one of two formats:
1. A list of point coordinates (x, y) with two columns.
2. A single line with four values representing limits (xmin, xmax, ymin, ymax).
The input file can be delimited by either tabs or commas and may include a header line. 

"""
function create_gridpts_from_file(infile::String,spacing::Float64,taskID::String="task",outdir::String=".",buffer::Float64=0)

    # check format of infile
    lines = readlines(infile)
    first_data_line = nothing
    for line in lines
        s = strip(line)
        if !isempty(s)
            first_data_line = s
            break
        end
    end

    n_tabs = count(==("\t"[1]), first_data_line)
    n_commas = count(==(','), first_data_line)
    delim = n_tabs > n_commas ? '\t' : ','

    first_fields = split(first_data_line, delim)
    has_header = any(isnothing(tryparse(Float64, strip(field))) for field in first_fields)
    skiprows = has_header ? 1 : 0

    # load the data
    pts = readdlm(infile, delim, Float64, '\n'; skipstart=skiprows)

    # calculate the limits based on the input data format        
    if size(pts,2) == 2
        xmin = floor(minimum(pts[:,1]) + buffer)
        xmax = ceil(maximum(pts[:,1]) - buffer)
        ymin = floor(minimum(pts[:,2]) + buffer)
        ymax = ceil(maximum(pts[:,2]) - buffer)
    elseif size(pts,2) == 4
        xmin = pts[1,1] + buffer
        xmax = pts[1,2] - buffer
        ymin = pts[1,3] + buffer
        ymax = pts[1,4] - buffer
    end

    xpts = collect(xmin:spacing:xmax)
    ypts = collect(ymin:spacing:ymax)
    gridpts = [(x,y) for x in xpts for y in ypts]

    outpath = joinpath(outdir,"$(taskID)_gridpts.txt")
    writedlm(outpath,gridpts)

end

"""
    create_gridpts_from_raster(infile::String, spacing::Float64, taskID::String="task", outdir::String=".", buffer::Float64=0)

Create a regular grid of points based on the spatial extent of a raster file and save them to a file.

# Arguments
- `infile::String`: Path to the input raster file (e.g., GeoTIFF)
- `spacing::Float64`: Distance between grid points in both x and y directions
- `taskID::String`: Identifier for the output file (default: "task")
- `outdir::String`: Output directory path (default: current directory ".")
- `buffer::Float64`: Optional buffer to add (positive) or subtract (negative) from the raster extent (default: 0)

# Returns
Nothing. The function writes grid points to a delimited text file.
Output values are in the same coordinate system as the input limits.
Output values are whole numbers, so any decimal values in the limits are rounded down for minimums and rounded up for maximums.

# Requirements
ArchGDAL.jl package is required to read the raster file and extract its spatial extent.
DelimitedFiles.jl package is required to write the output file.

# Details
The function reads the spatial extent of the input raster file using ArchGDAL, applies an optional buffer to the limits.
Generates a regular grid of points within that extent using the specified spacing.

"""
function create_gridpts_from_raster(infile::String, spacing::Float64, taskID::String="task", outdir::String=".", buffer::Float64=0)

    # load the raster and get the limits
    ArchGDAL.read(infile) do ds
        gt = ArchGDAL.getgeotransform(ds)
        ncols = ArchGDAL.width(ds)
        nrows = ArchGDAL.height(ds)
        xmin = gt[1]
        ymax = gt[4]
        xmax = gt[1] + gt[2]*ncols
        ymin = gt[4] + gt[6]*nrows

        xpts = collect(floor(xmin-buffer):spacing:ceil(xmax+buffer))
        ypts = collect(floor(ymin-buffer):spacing:ceil(ymax+buffer))
        gridpts = [(x,y) for x in xpts for y in ypts]

        outpath = joinpath(outdir,"$(taskID)_gridpts.txt")
        writedlm(outpath,gridpts)
    end

end