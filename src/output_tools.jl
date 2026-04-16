"""
    create_grid_from_file(infile::String)

Takes an output file from CanRad (created using `createfiles`) and creates a grid of points based on the coordinates in the file.
Uses the minimum and maximum coordinates and the spacing between points to create a complete grid of points.

Requirements:
Grid must be complete (will fail if there are missing points)

Assumptions:
- point spacing is consistent across the grid
- the grid has the same x and y dimensions (i.e. square cellsize)

"""
function create_grid_from_file(infile::String)

    ds = NCDataset(infile,"r")
    xvals = ds["easting"][:]
    yvals = ds["northing"][:]

    easting = sort(unique(xvals))
    northing = sort(unique(yvals))

    @assert size(unique(diff(sort(unique(xvals)))))[1] == 1 "spacing between points inconsistent in x dimension"
    @assert size(unique(diff(sort(unique(yvals)))))[1] == 1 "spacing between points inconsistent in y dimension"

    pt_spacing_x = diff(sort(unique(xvals)))[1]
    pt_spacing_y = diff(sort(unique(yvals)))[1]

    @assert pt_spacing_x == pt_spacing_y "cellsize not consistent across x and y dimensions"

    xmin = minimum(xvals) - pt_spacing_x/2
    xmax = maximum(xvals) + pt_spacing_x/2
    ymin = minimum(yvals) - pt_spacing_y/2
    ymax = maximum(yvals) + pt_spacing_y/2

    fulltilesize = Int64(max(xmax - xmin, ymax - ymin))

    @assert fulltilesize == (xmax - xmin) == (ymax - ymin) "Grid is not square"

    numptsgrid = Int64(div(fulltilesize,pt_spacing_x))
    numptsvec = numptsgrid^2

    # run some checks
    @assert numptsvec == size(xvals,1) == size(yvals,1) "Number of points in grid does not match number of points in input file"
    
    B = hcat(xvals,yvals)
    p = sortperm(view.(Ref(B), 1:size(B,1), :));

    # get a list of the variable names in ds
    list_vars = keys(ds)

    numptstime = any(contains.(list_vars, "datetime")) ? size(ds["datetime"][:], 1) : 1

    exdir = splitpath(infile)[1:end-1] |> joinpath
    outname = "Output_" * splitpath(infile)[end][8:end-3] * "_grid.nc"
    ds_grid = NCDataset(joinpath(exdir,outname),"c",format=:netcdf4_classic)
    defDim(ds_grid,"locX",size(easting,1)); defDim(ds_grid,"locY",size(northing,1))
    epsg_code = haskey(ds["easting"].attrib,"epsg") ? ds["easting"].attrib["epsg"] : missing
    defVar(ds_grid,"easting",Float64.(easting),("locX",),attrib=[
            "long_name" => "easting coordinate of center of pixel",
            "epsg"  => epsg_code,]); 
    defVar(ds_grid,"northing",Float64.(reverse(northing)),("locY",),attrib=[
            "long_name" => "northing coordinate of center of pixel",
            "epsg"  => epsg_code,])

    if any(contains.(list_vars, "datetime"))
        defVar(ds_grid,"datetime",ds["datetime"][:],("datetime",),attrib=ds["datetime"].attrib)
    end

    list_vars = filter(var -> any(startswith(string(var), prefix) for prefix in ["svf", "tvt", "swr"]), list_vars)
    
    for var in list_vars
        println("Processing variable: ", var)

        if contains(string(var),"svf")

            dat = nomissing(ds[string(var)][:])
            dat_new = Int8.(reverse(reshape(dat[p],(numptsgrid,numptsgrid)),dims=1))
            if haskey(ds[string(var)].attrib, "_FillValue")
                src_attrib = ds[string(var)].attrib
                attrib_no_fill = Dict{String,Any}(
                    String(k) => v for (k, v) in pairs(src_attrib) if String(k) != "_FillValue"
                )
                defVar(ds_grid,string(var),dat_new,("locY","locX",),deflatelevel=5,attrib=attrib_no_fill)
            else
                defVar(ds_grid,string(var),dat_new,("locY","locX",),deflatelevel=5,fillvalue=Int8(-1),attrib=ds[string(var)].attrib)
            end

        elseif contains(string(var),"tvt") || contains(string(var),"swr")

            is_swr = contains(String(var), "swr")
            T = is_swr ? Int16 : Int8
            fillv = T(-1)

            dat = nomissing(ds[string(var)][:,:])
            dat_new = fill(fillv,numptstime,numptsgrid,numptsgrid)
            for x in 1:1:size(dat,1)
                dat_new[x,:,:] = reverse(reshape(dat[x,p],numptsgrid,numptsgrid),dims=1)
            end
            if haskey(ds[string(var)].attrib, "_FillValue")
                src_attrib = ds[string(var)].attrib
                attrib_no_fill = Dict{String,Any}(
                    String(k) => v for (k, v) in pairs(src_attrib) if String(k) != "_FillValue"
                )
                defVar(ds_grid,string(var),T.(dat_new),("datetime","locY","locX",),deflatelevel=5,attrib=attrib_no_fill)
            else
                defVar(ds_grid,string(var),T.(dat_new),("datetime","locY","locX",),deflatelevel=5,fillvalue = T(-1),attrib=ds[string(var)].attrib)
            end
        end

    end

    close(ds)
    close(ds_grid)


end


"""
    create_grid_from_multiple_files(infile::String)


Takes a list of output files from CanRad (created using `createfiles`) and creates a grid of points based on the coordinates in all the files.
Uses the minimum and maximum coordinates and the spacing between points to create a complete grid of points.
Missing values in the grid are replaced with nan. 

Requirements:
All files must have the same variables

Assumptions:
- point spacing is consistent across the grid

Arguments:
- `allfiles::Vector{String}`: A vector of file paths to the input NetCDF files.
    e.g. obtained from `glob(joinpath(exdir,"*","*","*.nc"))` where `exdir` is the directory containing the output from CanRad.
- `par_in::Dict{String, Any}`: The parameters dictionary used in the CanRad runs. This is used to determine which variables are present and to set attributes in the output file.
- `outfile::String`: Optional. The path to the output NetCDF file. If not provided, a default name is generated.

Note that depending on the number of total points and/or the variables, this can take some time to run and can use a lot of memory. 

Usage:
See examples/gridded in CanRad repository of example usage

"""

function create_grid_from_multiple_files(allfiles::Vector{String},par_in::Dict{String, Any},outfile::String="")

    if outfile == ""
        outputfile = split(allfiles[1],"/")[1]*"_gridded.nc"
    else
        outputfile = outfile
    end

    # load in the metadata
    ds0 = NCDataset(allfiles[1],"r")

    datetime = ds0["datetime"][:]
    numptstime = length(datetime)

    # get attributes
    var_attrib = Dict{String, Dict{String, Any}}()
    for vname in keys(ds0)
        d = Dict{String, Any}()
        for (k, v) in pairs(ds0[vname].attrib)
            string(k) == "_FillValue" && continue
            d[string(k)] = v
        end
        var_attrib[vname] = d
    end

    # initialise all possible variables
    east = Float64[]
    north = Float64[]

    envstrings = environments_flag_local(par_in)

    evergreen = false
    terrain = false
    leafoff = false
    leafon = false
    tvt = false
    swr = false

    for envname in envstrings

        if envname == "evergreen"
            evergreen = true
            global svf_planar_evergreen = []
            global svf_hemi_evergreen = []
            global tvt_evergreen = []
            global swr_total_evergreen = []
            global swr_direct_evergreen = []

            (sum(contains.(keys(ds0),"tvt")) > 0) && (tvt = true)
            (sum(contains.(keys(ds0),"swr")) > 0) && (swr = true)
        end

        if envname == "terrain"
            terrain = true
            global svf_planar_terrain = []
            global svf_hemi_terrain = []
            global tvt_terrain = []
            global swr_total_terrain = []
            global swr_direct_terrain = []
            (sum(contains.(keys(ds0),"tvt")) > 0) && (tvt = true)
            (sum(contains.(keys(ds0),"swr")) > 0) && (swr = true)
        end

        if envname == "leafoff"
            leafoff = true
            global svf_planar_leafoff = []
            global svf_hemi_leafoff = []   
            global tvt_leafoff = []
            global swr_total_leafoff = []
            global swr_direct_leafoff = []  
            (sum(contains.(keys(ds0),"tvt")) > 0) && (tvt = true)
            (sum(contains.(keys(ds0),"swr")) > 0) && (swr = true)
        end 

        if envname == "leafon"
            leafon = true
            global svf_planar_leafon = []
            global svf_hemi_leafon = []    
            global tvt_leafon = []
            global swr_total_leafon = []
            global swr_direct_leafon = []
            (sum(contains.(keys(ds0),"tvt")) > 0) && (tvt = true)
            (sum(contains.(keys(ds0),"swr")) > 0) && (swr = true)
        end

    end

    close(ds0)


    for t in allfiles[1:end-3]

            tds = NCDataset(t)

            append!(east, tds["easting"][:])
            append!(north, tds["northing"][:])

            if terrain
                append!(svf_planar_terrain,tds["svf_planar_terrain"][:])
                append!(svf_hemi_terrain,tds["svf_hemi_terrain"][:])
                tvt && append!(tvt_terrain,vec(tds["tvt_terrain"][:,:]))
                swr && append!(swr_total_terrain,vec(tds["swr_total_terrain"][:,:]))
                swr && append!(swr_direct_terrain,vec(tds["swr_direct_terrain"][:,:]))
            end

            if evergreen
                append!(svf_planar_evergreen,tds["svf_planar_evergreen"][:])
                append!(svf_hemi_evergreen,tds["svf_hemi_evergreen"][:])
                tvt && append!(tvt_evergreen,vec(tds["tvt_evergreen"][:,:]))
                swr && append!(swr_total_evergreen,vec(tds["swr_total_evergreen"][:,:]))
                swr && append!(swr_direct_evergreen,vec(tds["swr_direct_evergreen"][:,:]))
            end

            if leafoff
                append!(svf_planar_leafoff,tds["svf_planar_leafoff"][:])
                append!(svf_hemi_leafoff,tds["svf_hemi_leafoff"][:])
                tvt && append!(tvt_leafoff,vec(tds["tvt_leafoff"][:,:]))
                swr && append!(swr_total_leafoff,vec(tds["swr_total_leafoff"][:,:]))
                swr && append!(swr_direct_leafoff,vec(tds["swr_direct_leafoff"][:,:]))
            end

            if leafon
                append!(svf_planar_leafon,tds["svf_planar_leafon"][:])
                append!(svf_hemi_leafon,tds["svf_hemi_leafon"][:])
                tvt && append!(tvt_leafon,vec(tds["tvt_leafon"][:,:]))
                swr && append!(swr_total_leafon,vec(tds["swr_total_leafon"][:,:]))
                swr && append!(swr_direct_leafon,vec(tds["swr_direct_leafon"][:,:]))
            end

            close(tds)
            println("loaded data from $(t)")
    end

    # get information for creating the new grid structure
    
    B = hcat(east,north)
    p = sortperm(view.(Ref(B), 1:size(B,1), :))
    # p = sortperm(Tuple.(eachrow(B)))

    # now reorder to grid format
    @assert size(unique(diff(sort(unique(east)))))[1] == 1 "spacing between points inconsistent in x dimension"
    @assert size(unique(diff(sort(unique(north)))))[1] == 1 "spacing between points inconsistent in y dimension"

    pt_spacing_x = diff(sort(unique(east)))[1]
    pt_spacing_y = diff(sort(unique(north)))[1]

    xvals = sort(unique(east))
    yvals = sort(unique(north))

    numpts_x = length(xvals)
    numpts_y = length(yvals)

    full_grid_size = length(xvals) * length(yvals)

    if full_grid_size != size(B,1)
        full_grid = reduce(vcat, [[x y] for x in minimum(east):pt_spacing_x:maximum(east) for y in minimum(north):pt_spacing_y:maximum(north)])
        bset = Set(Tuple.(eachrow(B)))
        locs = [Tuple(r) in bset for r in eachrow(full_grid)]
    else
        locs = trues(size(B,1))
    end

    numptsvec = size(locs,1)

    # rearrange the data and save to the output file

    println("creating gridded file at "*outputfile)

    ds_out = NCDataset(outputfile,"c",format=:netcdf4_classic)
    defDim(ds_out,"locX",numpts_x); defDim(ds_out,"locY",numpts_y)
    defVar(ds_out,"easting",Float64.(unique(sort(east))),("locX",),attrib=[
            "long_name" => "easting coordinate of centre of pixel",
            "epsg"  => par_in["epsg_code"],]); 
    defVar(ds_out,"northing",Float64.(reverse(unique(sort(north)))),("locY",),attrib=[
            "long_name" => "northing coordinate of centre of pixel",
            "epsg"  => par_in["epsg_code"],])
    defDim(ds_out,"DateTime",size(datetime,1))
    defVar(ds_out,"datetime",datetime,("DateTime",),attrib=var_attrib["datetime"])

    
    if evergreen
        ds_out = write_grid_to_netcdf(ds_out,
                            svf_planar_evergreen,svf_hemi_evergreen,
                            tvt,tvt_evergreen,
                            swr,swr_total_evergreen,swr_direct_evergreen,
                            numptstime,numptsvec,numpts_x,numpts_y,
                            "evergreen",var_attrib,locs,p)
    end

    if terrain
        ds_out = write_grid_to_netcdf(ds_out,
                            svf_planar_terrain,svf_hemi_terrain,
                            tvt,tvt_terrain,
                            swr,swr_total_terrain,swr_direct_terrain,
                            numptstime,numptsvec,numpts_x,numpts_y,
                            "terrain",var_attrib,locs,p)
    end

    if leafoff
        ds_out = write_grid_to_netcdf(ds_out,
                            svf_planar_leafoff,svf_hemi_leafoff,
                            tvt,tvt_leafoff,
                            swr,swr_total_leafoff,swr_direct_leafoff,
                            numptstime,numptsvec,numpts_x,numpts_y,
                            "leafoff",var_attrib,locs,p)
    end

    if leafon
        ds_out = write_grid_to_netcdf(ds_out,
                            svf_planar_leafon,svf_hemi_leafon,
                            tvt,tvt_leafon,
                            swr,swr_total_leafon,swr_direct_leafon,
                            numptstime,numptsvec,numpts_x,numpts_y,
                            "leafon",var_attrib,locs,p)
    end

    close(ds_out)

    println("created gridded file at "*outputfile)


end

"""
    write_grid_to_netcdf(ds_out::NCDataset,
                            svf_planar,svf_hemi,
                            tvt,tvt_dat,
                            swr,swr_total,swr_direct,
                            numptstime,numptsvec,numpts_x,numpts_y,
                            envname,var_attrib,locs,p)

Helper function to write the data for a single environment to the output netcdf file in the correct format.

Called by `create_grid_from_multiple_files`.
"""
function write_grid_to_netcdf(ds_out::NCDataset,
                            svf_planar,svf_hemi,
                            tvt,tvt_dat,
                            swr,swr_total,swr_direct,
                            numptstime,numptsvec,numpts_x,numpts_y,
                            envname,var_attrib,locs,p)

        tvt && (tvt_dat = permutedims(reshape(tvt_dat, numptstime, :)))
        swr && (swr_total = permutedims(reshape(swr_total, numptstime, :)))
        swr && (swr_direct = permutedims(reshape(swr_direct, numptstime, :)))

        # lazy re-allocating array even if full_grid = true
        svf_planar_out = fill(-1,numptsvec,1)
        svf_planar_out[locs] .= svf_planar[p]
        defVar(ds_out,"svf_planar_$envname",Int8.(reverse(reshape(svf_planar_out,numpts_y,numpts_x),dims=1)),("locY","locX",),deflatelevel=5,attrib=var_attrib["svf_planar_$envname"])

        svf_hemi_out = fill(-1,numptsvec,1)
        svf_hemi_out[locs] .= svf_hemi[p]
        defVar(ds_out,"svf_hemi_$envname",Int8.(reverse(reshape(svf_hemi_out,numpts_y,numpts_x),dims=1)),("locY","locX",),deflatelevel=5,attrib=var_attrib["svf_hemi_$envname"])

        if tvt
            tvt_dat_out = fill(-1,numptsvec,numptstime)
            tvt_dat_out[locs,:] .= tvt_dat[p,:]
            tvt_dat_newdim = reverse(reshape(tvt_dat_out', numptstime, numpts_y, numpts_x),dims=2)
            defVar(ds_out,"tvt_$envname",Int8.(tvt_dat_newdim),("DateTime","locY","locX",),fillvalue = Int8(-1),attrib=var_attrib["tvt_$envname"])
        end

        if swr
            swr_total_out = fill(-1,numptsvec,numptstime)
            swr_total_out[locs,:] .= swr_total[p,:]
            swr_total_newdim = reverse(reshape(swr_total_out', numptstime, numpts_y, numpts_x),dims=2)
            defVar(ds_out,"swr_total_$envname",Int32.(swr_total_newdim),("DateTime","locY","locX",),fillvalue = Int32(-1),attrib=var_attrib["swr_total_$envname"])

            swr_direct_out = fill(-1,numptsvec,numptstime)
            swr_direct_out[locs,:] .= swr_direct[p,:]
            swr_direct_newdim = reverse(reshape(swr_direct_out', numptstime, numpts_y, numpts_x),dims=2)
            defVar(ds_out,"swr_direct_$envname",Int32.(swr_direct_newdim),("DateTime","locY","locX",),fillvalue = Int32(-1),attrib=var_attrib["swr_direct_$envname"])
        end

    return ds_out

end

"""
    environments_flag_local(par_in::Dict{String, Any})

Helper function called by `create_grid_from_multiple_files` to determine which variables to read in and write to the output file.
"""
function environments_flag_local(par_in)

    envstrings = String[]
    par_in["calc_terrain"] && push!(envstrings, "terrain")
    (par_in["phenology"] == "leafon" || par_in["phenology"] == "both") && push!(envstrings, "leafon")
    (par_in["phenology"] == "leafoff" || par_in["phenology"] == "both") && push!(envstrings, "leafoff")
    par_in["forest_type"] == "evergreen" && push!(envstrings, "evergreen")

    return envstrings

end

