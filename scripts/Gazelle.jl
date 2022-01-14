# get basinwide mean values for the Indian Ocean
using Revise, HistoricalIndianOcean, DrWatson, GoogleDrive, NCDatasets

# Download data from GoogleDrive
#inputfile = datadir("TMI_"*TMIversion*".nc")

# HistoricalIndianOcean.nc
# https://drive.google.com/file/d/1-mOo6dwHVwv0TJMoFiYR5QWxykJhNwWK/view?usp=sharing
file_id = "1-mOo6dwHVwv0TJMoFiYR5QWxykJhNwWK"
#url = "https://drive.google.com/file/d/"*file_id*"/view?usp=sharing"
url = "https://docs.google.com/uc?export=download&id="*file_id
filename_tmp = google_download(url,datadir())
#filename=datadir("HistoricalIndianOcean.nc")
#mv(filename_tmp,filename,force=true)
nc = NCDataset(filename)


locs = Vector{Loc}(undef,2)

# lon, lat, depth
loc1 = Loc(90, 0, 1200)
loc2 = Loc(89, -1, 2200)

# this works, but doesn't specify types
#NamedTuple{(:lon,:lat,:depth)}(loc1)
#T = Float64
# loc = Vector{@NamedTuple{lon::T,lat::T,depth::T}}(undef,2)
# loc[1] = @NamedTuple{lon::T,lat::T,depth::T}((loc1))
# loc[2] = @NamedTuple{lon::T,lat::T,depth::T}((loc2))

# assign the Vector
locs[1] = loc1
locs[2] = loc2

Lxy_decadal = 450_000 #km
Lz_decadal = 450 #m

# try the same thing for watermass gradients.
Lxy_watermass = 2_000_000 # m
#                      %xylengthscale = 3e3;
Lz_watermass = 1000 # m
#                      %zlengthscale  = 1.5e3;% 
TTerr = 0.14^2
    
Rqq = error_covariance(locs,TTerr,Lxy_decadal,Lz_decadal, Lxy_watermass, Lz_watermass)

Lxy = 450_000 # meters 
Lz = 450 # meters 


# take synthetic observations
# get observational uncertainty
σθ = readtracer(inputfile,"σ"*variable)

# testing Named Tuples
# locs = NamedTuple{(:lon,:lat,:depth),T}()
# locs = (lon = lon[1], lat=lat[1], depth=depth[1])
# names = (:lon,:lat,:depth)
# loc = NamedTuple{names,(T,T,T)}()
# for i in eachindex(lon)
#     loc[i] = (lon = lon[i], lat=lat[i], depth=depth[i])
# end
