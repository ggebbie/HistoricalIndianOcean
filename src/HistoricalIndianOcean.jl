module HistoricalIndianOcean

using DrWatson, Distances, TMI, Interpolations, LinearAlgebra

export R2_covariance, Loc, observational_covariance,
    weighted_covariance, error_covariance, woce_error

import TMI.interpindex

println("new module for functions")

struct Loc
    lon::Float64
    lat::Float64
    depth::Float64
end

"""
function error_covariance(locs,watermassvar,Lxy_decadal,Lz_decadal, Lxy_watermass, Lz_watermass)

    locs: locations of data
    σobs: observational covariance
    watermassvar: fraction of water-mass variability to total (heaving, waves) variability 
"""
function error_covariance(locs,σobs,watermassvar,Lxy_decadal,Lz_decadal, Lxy_watermass, Lz_watermass)

    # hard-wire the TMI version
    TMIversion = "modern_90x45x33_GH10_GH12"
    A, Alu, γ, TMIfile, L, B = TMI.config_from_nc(TMIversion)
    ση = woce_error(TMIversion,locs,γ)

    # multiply expected error by water-mass variability fraction
    Rss = 2*weighted_covariance(locs,watermassvar*ση,Lxy_watermass,Lz_watermass)
    Rtt =  2*weighted_covariance(locs,ση,Lxy_decadal,Lz_decadal)

    Rmm = observational_covariance(σobs)
    Rqq  = Rss + Rtt + Rmm
    
  # %% use 0.2 K for uncertainty. Instead of 0.14 K.
  # Rqq = (0.14.^2 + 0.03.^2).*eye(Nobs,Nobs) + 2.*RSxx + 2.*RTxx;
  # iRqq = inv(Rqq);
    return Rqq
end

# function haversined(loc1::Loc,loc2::Loc)
#     return distance = haversine((loc1.lon,loc1.lat),(loc2.lon,loc2.lat))
# end
haversine(loc1::Loc,loc2::Loc) = Distances.haversine((loc1.lon,loc1.lat),(loc2.lon,loc2.lat))

function weighted_covariance(locs,ση,Lxy,Lz)

    Rρ = R2_covariance(locs,Lxy,Lz)
    Rηη = Rρ.* (ση * ση')

    return Rηη
    
end

observational_covariance(σobs) = σobs^2 * I

function R2_covariance(locs,Lxy,Lz)

    # get haversine distance
    C = Matrix{Float64}(undef,length(locs),length(locs))

    for i in eachindex(locs)
        for j in eachindex(locs)
            rxy = (haversine(locs[i],locs[j])/Lxy)^2
            rz  = ((locs[i].depth-locs[j].depth)/Lz)^2
            #   exponent = (xydist./xylengthscale).^2 + (zdist./zlengthscale).^2;
            C[i,j] = exp(-(rxy+rz))
        end
    end
    return C 
end

# Get WOCE error as a decadal proxy at locations
function woce_error(TMIversion,locs,γ)

    # read WOCE error
    inputfile = datadir("TMI_"*TMIversion*".nc")

    # take synthetic observations
    # get observational uncertainty
    
    σθ = TMI.readtracer(inputfile,"σθ")

    # get weighted interpolation indices
    # equivalent to TMI E matrix
    N = length(locs)
    wis= Vector{Tuple{Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}, Interpolations.WeightedAdjIndex{2, Float64}}}(undef,N)
    [wis[i] = interpindex(locs[i],γ) for i in 1:N]

    # get WOCE error at this location
    σsample =  observe(σθ,wis,γ)

    return σsample

end

"""
    function interpindex(loc,γ)

    add method to TMI.interpindex
    all Loc type to work with interpindex
"""
function interpindex(loc::Loc,γ::TMI.grid)

    loctuple = (loc.lon,loc.lat,loc.depth)
    wis = interpindex(loctuple,γ)
        
end

end # module HistoricalIndianOcean
