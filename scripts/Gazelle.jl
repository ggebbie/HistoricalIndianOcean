# get basinwide mean values for the Indian Ocean
using Revise
using HistoricalIndianOcean, DrWatson, GoogleDrive, NCDatasets, TMI, PyPlot

# input parameters
Lxy_decadal = 450_000 #m
Lz_decadal = 450 #m

# try the same thing for watermass gradients.
Lxy_watermass = 2_000_000 # xylengthscale = 3e3 m
Lz_watermass = 1000 # zlengthscale  = 1.5e3 m  

watermassvar = (1/5.)^2 
σobs = 0.54  # deg C # obs error from HMS Challenger

Lz_basinwideavg = 1000 # meters

# For the purposes of the line plot in the current version of the manuscript I binned the data into the following depth ranges, which give a reasonably balanced set of obs in each range
zgridedge = [0, 100, 300, 500, 1000, 2000, 3000]
zgrid = [0, 50, 200, 400, 750, 1500, 2500]

## Download TMI grid and error data from GoogleDrive
#inputfile = datadir("TMI_"*TMIversion*".nc")
#TMIversion = "modern_90x45x33_GH10_GH12"
#A, Alu, γ, TMIfile, L, B = TMI.config_from_nc(TMIversion)

# HistoricalIndianOcean.nc: data from Gazelle, etc.
# https://drive.google.com/file/d/1-mOo6dwHVwv0TJMoFiYR5QWxykJhNwWK/view?usp=sharing
file_id = "1-mOo6dwHVwv0TJMoFiYR5QWxykJhNwWK"
url = "https://docs.google.com/uc?export=download&id="*file_id
filename = google_download(url,datadir())

nc = NCDataset(filename)

igood = findall(!ismissing,nc["delta_T"])
nobs = count(!ismissing,nc["delta_T"])

locs = Vector{Loc}(undef,nobs)

# lon, lat, depth
#[locs[r] = Loc(nc["lon"][r],nc["lat"][r],nc["depth"][r]) for r in eachindex(nc["delta_T"])]
[locs[r] = Loc(nc["lon"][igood[r]],nc["lat"][igood[r]],nc["depth"][igood[r]]) for r in eachindex(igood)]

Rqq = error_covariance(locs,σobs,watermassvar,Lxy_decadal,Lz_decadal, Lxy_watermass, Lz_watermass)

zobs = collect(nc["depth"][igood])
S = vertical_smoothness(zgrid,Lz_basinwideavg)

# make H matrix, vertical map onto grid

H = vertical_map(zobs,zgrid)

#################################
# not good because type includes missing
#ΔT = nc["delta_T"][igood]

# here ΔT is type vector float (no missing)
ΔT = zeros(nobs)
for ii in eachindex(igood)
    println(ii)
    ΔT[ii] = nc["delta_T"][igood[ii]]
end

#figure()
#plot(ΔT,-nc["depth"],marker="o")

ΔT̄,σT̄ = basinwide_avg(ΔT,Rqq,H,S)

figure()
plot(ΔT̄,-zgrid)
plot(ΔT̄.+σT̄,-zgrid)
plot(ΔT̄.-σT̄,-zgrid)
grid()

ΔTmean = mean(skipmissing(ΔT))
ΔTstd = std(skipmissing(ΔT))

err_naive = ΔTstd/sqrt(count(!ismissing,ΔT))


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
