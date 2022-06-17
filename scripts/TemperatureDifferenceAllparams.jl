# get basinwide mean values for the Indian Ocean

# activate environment with DrWatson
include("intro.jl")

using Revise
using HistoricalIndianOcean, DrWatson, TMI, CSV, DataFrames

# to make figures
using Plots
plotlyjs()

# What temperature field?
delta = ["T","T_tait","T_5564"]

# input parameters
LxyT = 450_000 # m
LzT = 450 # m

# try the same thing for watermass gradients.
LxyS = 2_000_000 # xylengthscale = meters
LzS = 1000 # zlengthscale  = meters

sratio = [1/25, 1/4] 
σobs = 0.14  # deg C # obs error from HMS Challenger

LzAVG = [500] # meters

tratio = [2.0,3.0,4.0] # amplify variability expected by a decadal average

σS = [1.0] # first-guess size of anomalies, deg C

#latsouth = [-50,-40,-30]
latsouth = [-50]

## Next set the fixed variables ################

# For the purposes of the line plot in the current version of the manuscript I binned the data into the following depth ranges, which give a reasonably balanced set of obs in each range
zgridedge = [-1, 1, 100, 300, 500, 1000, 2000, 3000]
zgridmid = [0, 50, 200, 400, 750, 1500, 2500]
zgrid = Vector{Any}(undef,1)
zgrid[1] = zgridmid 

zstar = 700

# save computation by loading TMI grid as a fixed variable.
# only γ is needed
TMIversion = "modern_90x45x33_GH10_GH12"
A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)

# Several parameter containers
allparams = @strdict delta σobs tratio sratio LxyT LzT LxyS LzS σS LzAVG latsouth zgrid zstar γ TMIversion

accessvars = ["delta","tratio","sratio","latsouth"]

dicts = dict_list(allparams)

for (i, d) in enumerate(dicts)

    output = basinwide_avg(d)

    println("H=",output["H"])
    println("σH=",output["σH"])

    # output full state of analysis to jld2
    !isdir(datadir("all")) && mkdir(datadir("all"))
    @tagsave(datadir("all",savename("DTbar",d,"jld2",accesses=accessvars)), output)

    # output profile information to csv
    xax = "ΔT̄ [°C]"
    yax = "z [m]"
    zoutput = Dict(yax => output["zgrid"], xax => output["T̄"], "σΔT̄ [°C]" => output["σT̄"])
    df = DataFrame(zoutput)
    println(df)
    #!isdir(datadir("csv")) && mkdir(datadir("csv"))
    CSV.write(datadir("all",savename("DTbar",d,"csv",accesses=accessvars)),df)

    dfscalar = DataFrame(Dict("ΔH [ZJ]" => output["H"],"σ(ΔH) [ZJ]" => output["σH"],"z⋆ [m]" => zstar ))
    CSV.write(datadir("all",savename("DH",d,"csv",accesses=accessvars)),dfscalar)

    titlelabel = replace(savename(d,accesses=["delta"]),"_" => " ")
    plot(output["T̄"],zgrid[1],xerr=output["σT̄"],
         xlabel=xax,ylabel=yax,title=titlelabel,yflip=true,legend=false)

    !isdir(plotsdir("all")) && mkpath(plotsdir("all"))
    figname = plotsdir("all",savename("DTbar",d,"pdf",accesses=accessvars))
    savefig(figname)

end


# Left to do: compare to bin-averaged ΔT
#ΔTmean = mean(skipmissing(ΔT))
#ΔTstd = std(skipmissing(ΔT))
#err_naive = ΔTstd/sqrt(count(!ismissing,ΔT))
