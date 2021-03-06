# get basinwide mean values for the Indian Ocean

# activate environment with DrWatson
include("intro.jl")

using Revise
using HistoricalIndianOcean, DrWatson, TMI, CSV, DataFrames

# to make plots
using Plots
plotlyjs() # gr() is simpler but has problems on Wayland/Ubuntu 21.10

# What temperature field?
delta = ["T","T_tait","T_5564"]

# input parameters
LxyT = 450_000 # m
LzT = 450 # m

# try the same thing for watermass gradients.
LxyS = 2_000_000 # xylengthscale = meters
LzS = 1000 # zlengthscale  = meters

sratio = (1/5.0)^2

σobs = 0.14  # °C: obs error from HMS Challenger

LzAVG = 500 # meters

tratio = 3.0 # amplify variability expected by a decadal average

σS = 1.0 # first-guess size of anomalies, deg C

latsouth = -50

## Next set the fixed variables ################

# For the purposes of the line plot in the current version of the manuscript I binned the data into the following depth ranges, which give a reasonably balanced set of obs in each range
zgridedge = [-1, 1, 100, 300, 500, 1000, 2000, 3000]
zgridmid = [0, 50, 200, 400, 750, 1500, 2500]
zgrid = Vector{Any}(undef,1)
zgrid[1] = zgridmid 

zstar = 700 # meters, depth over which heat content difference is calculated

# save computation by loading TMI grid as a fixed variable.
# only γ is needed
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# Several parameter containers
bestparams = @strdict delta σobs tratio sratio LxyT LzT LxyS LzS σS LzAVG latsouth zgrid zstar γ TMIversion 

dicts = dict_list(bestparams)

for (i, d) in enumerate(dicts)

    output = basinwide_avg(d)

    println("H=",output["H"])
    println("σH=",output["σH"])
    
    # output full state of analysis to jld2
    !isdir(datadir("best")) && mkdir(datadir("best"))
    @tagsave(datadir("best",savename("DTbar",d,"jld2",accesses=["delta"])), output)

    # output profile information to csv
    xax = "ΔT̄ [°C]"
    yax = "z [m]"
    zoutput = Dict(yax => output["zgrid"], xax => output["T̄"], "σΔT̄ [°C]" => output["σT̄"])
    df = DataFrame(zoutput)
    println(df)
    CSV.write(datadir("best",savename("DTbar",d,"csv",accesses=["delta"])),df)

    dfscalar = DataFrame(Dict("ΔH [ZJ]" => output["H"],"σ(ΔH) [ZJ]" => output["σH"], "z⋆ [m]" => zstar))
    println("save to ",datadir("best",savename("DH",d,"csv",accesses=["delta"])))
    CSV.write(datadir("best",savename("DH",d,"csv",accesses=["delta"])),dfscalar)

    # make profile figure
    titlelabel = replace(savename(d,accesses=["delta"]),"_" => " ")
    plot(output["T̄"],zgrid[1],xerr=output["σT̄"],
         xlabel=xax,ylabel=yax,title=titlelabel,yflip=true,legend=false)

    # # make whole path if necessary
    !isdir(plotsdir("best")) && mkpath(plotsdir("best"))
    figname = plotsdir("best",savename("DTbar",d,"pdf",accesses=["delta"]))
    savefig(figname)

end

# tests
#readdir(datadir("csv"))
#readdir(plotsdir())

# Left to do: compare to bin-averaged ΔT
#ΔTmean = mean(skipmissing(ΔT))
#ΔTstd = std(skipmissing(ΔT))
#err_naive = ΔTstd/sqrt(count(!ismissing,ΔT))
