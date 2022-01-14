module HistoricalIndianOcean

using DrWatson, Distances

export R2_covariance, Loc, decadal_covariance,
    watermass_covariance, error_covariance

println("new module for functions")

struct Loc
    lon::Float64
    lat::Float64
    depth::Float64
end

function watermass_covariance(locs,TTerr,Lxy,Lz)
    Rxx = R2_covariance(locs,Lxy,Lz)
    return RSxx = TTerr*Rxx/25.
end

function error_covariance(locs,TTerr,Lxy_decadal,Lz_decadal, Lxy_watermass, Lz_watermass)

    RSxx = watermass_covariance(locs,TTerr,Lxy_watermass,Lz_watermass)
    RTxx =  decadal_covariance(locs,Lxy_decadal,Lz_decadal)
    Rqq  = RSxx + RTxx
  # %% use 0.2 K for uncertainty. Instead of 0.14 K.
  # Rqq = (0.14.^2 + 0.03.^2).*eye(Nobs,Nobs) + 2.*RSxx + 2.*RTxx;
  # iRqq = inv(Rqq);
    return Rqq
end

# function haversined(loc1::Loc,loc2::Loc)
#     return distance = haversine((loc1.lon,loc1.lat),(loc2.lon,loc2.lat))
# end
haversine(loc1::Loc,loc2::Loc) = Distances.haversine((loc1.lon,loc1.lat),(loc2.lon,loc2.lat))

decadal_covariance(locs,Lxy,Lz) = R2_covariance(locs,Lxy,Lz)

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
function woce_error()

    # read WOCE error
    
    # use TMI to get E matrix

    # get WOCE error at this location
    ETerr = E*Terr;

    # more stuff related to covariance?
    TTerr = ETerr*ETerr';
    RTxx = TTerr.*Rxx;
end


end # module
