# get basinwide mean values for the Indian Ocean

# activate environment with DrWatson
include("intro.jl")

using Revise
using HistoricalIndianOcean, DrWatson, TMI, PyPlot, CSV, DataFrames

# What temperature field?
delta = ["T","T_tait","T_5564"]

# input parameters
LxyT = 450_000 # m
LzT = 450 # m

# try the same thing for watermass gradients.
LxyS = 2_000_000 # xylengthscale = meters
LzS = 1000 # zlengthscale  = meters

sratio = (1/5.)^2 
σobs = 0.14  # deg C # obs error from HMS Challenger

LzAVG = [500, 1000] # meters

tratio = [1.0,2.0] # amplify variability expected by a decadal average

σS = [0.5,1.0] # first-guess size of anomalies, deg C

## Next set the fixed variables ################

# For the purposes of the line plot in the current version of the manuscript I binned the data into the following depth ranges, which give a reasonably balanced set of obs in each range
zgridedge = [-1, 1, 100, 300, 500, 1000, 2000, 3000]
zgrid = [0, 50, 200, 400, 750, 1500, 2500]

# Several parameter containers
allparams = @strdict delta σobs tratio sratio LxyT LzT LxyS LzS σS LzAVG zgrid
allparams["zgrid"] = [zgrid]
accessvars = ["delta","LzAVG","tratio","σS"]

dicts = dict_list(allparams)

for (i, d) in enumerate(dicts)

    output = basinwide_avg(d)

    # output full state of analysis to jld2
    @tagsave(datadir("all",savename("DTbar",d,"jld2",accesses=accessvars)), output)

    # output profile information to csv
    xax = "ΔT̄ [°C]"
    yax = "z [m]"
    zoutput = Dict(yax => output["zgrid"], xax => output["T̄"], "σΔT̄ [°C]" => output["σT̄"])
    df = DataFrame(zoutput)
    println(df)
    #!isdir(datadir("csv")) && mkdir(datadir("csv"))
    CSV.write(datadir("all",savename("DTbar",d,"csv",accesses=accessvars)),df)

    # make profile figure
    figure(i)
    clf()
    plot(output["T̄"],zgrid)
    errorbar(output["T̄"],zgrid,xerr=output["σT̄"])
    xlabel(xax)
    ylabel(yax)
    grid()
    gca().invert_yaxis()

    titlelabel = replace(savename(d,accesses=accessvars),"_" => " ")
    title(titlelabel,fontsize=10)

    !isdir(plotsdir("all")) && mkdir(plotsdir("all"))
    figname = plotsdir("all",savename("DTbar",d,"pdf",accesses=accessvars))
    savefig(figname)

end

# tests
#readdir(datadir("csv"))
#readdir(plotsdir())

# Left to do: compare to bin-averaged ΔT
#ΔTmean = mean(skipmissing(ΔT))
#ΔTstd = std(skipmissing(ΔT))
#err_naive = ΔTstd/sqrt(count(!ismissing,ΔT))

#=
# testing Named Tuples
# locs = NamedTuple{(:lon,:lat,:depth),T}()
# locs = (lon = lon[1], lat=lat[1], depth=depth[1])
# names = (:lon,:lat,:depth)
# loc = NamedTuple{names,(T,T,T)}()
# for i in eachindex(lon)
#     loc[i] = (lon = lon[i], lat=lat[i], depth=depth[i])
# end

# this works, but doesn't specify types
#NamedTuple{(:lon,:lat,:depth)}(loc1)
#T = Float64
# loc = Vector{@NamedTuple{lon::T,lat::T,depth::T}}(undef,2)
# loc[1] = @NamedTuple{lon::T,lat::T,depth::T}((loc1))
# loc[2] = @NamedTuple{lon::T,lat::T,depth::T}((loc2))
=#
