module HistoricalIndianOcean

using DrWatson, Distances, TMI,
    Interpolations, LinearAlgebra,
    GoogleDrive, NCDatasets

export R2_covariance, Loc, observational_covariance,
    weighted_covariance, error_covariance, woce_error,
    vertical_smoothness, vertical_map, basinwide_avg,
    read_historical_data, degree_meters, indianarea

import TMI.interpindex

println("new module for functions")

struct Loc
    lon::Float64
    lat::Float64
    depth::Float64
end

"""
function error_covariance(locs,sratio,LxyT,LzT, LxyS, LzS)

    locs: locations of data
    σobs: observational covariance
    sratio: fraction of water-mass variability to total (heaving, waves) variability 
"""
function error_covariance(locs,σobs,tratio,sratio,LxyT,LzT, LxyS, LzS)

    # hard-wire the TMI version
    TMIversion = "modern_90x45x33_GH10_GH12"
    A, Alu, γ, TMIfile, L, B = TMI.config_from_nc(TMIversion)

    ση = woce_error(TMIversion,locs,γ)

    # multiply expected error by water-mass variability fraction
    # tratio = how much bigger is sqrt of variability for this time interval
    # relative to the woce time interval.
    # For longer time intervals, tratio is larger.
    #Rtt =  weighted_covariance(locs,ση,LxyT,LzT)

    # assume independence between two time periods
    # add covariance together
    #Rss =  weighted_covariance(locs,√sratio*ση,LxyS,LzS) + weighted_covariance(locs,√sratio*√tratio*ση,LxyS,LzS)

    Rss = (tratio + 1) * sratio * weighted_covariance(locs,ση,LxyS,LzS)
    
    Rtt = (tratio + 1) * weighted_covariance(locs,ση,LxyT,LzT)
    
    Rmm = observational_covariance(σobs)
    
    Rqq  = Rss + Rtt + Rmm
    
    return Rqq
end

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

    # if NaN, replace with big value. Better solution: use 2 x 2 degree woce error
    replace!(σsample, NaN => 1.0)
    
    return σsample

end

"""
    function interpindex(loc,γ)

    add method to TMI.interpindex
    all Loc type to work with interpindex
"""
function interpindex(loc::Loc,γ::Grid)

    loctuple = (loc.lon,loc.lat,loc.depth)
    wis = interpindex(loctuple,γ)
        
end

function vertical_smoothness(zgrid,σ,Lz)
    
    # some smoothness in vertical profile.
    ngrid = length(zgrid)
    S = zeros(ngrid,ngrid)

    for nz in eachindex(zgrid)
        dz = zgrid[nz] .- zgrid
        exponent = (dz./Lz).^2
        S[nz,:] = σ^2 * exp.(-exponent);
    end

    # a small perturbation to the diagonal for
    # stability 
    S += 1e-8 * I

    return S
end

function vertical_map(zobs,zgrid)

    E = zeros(length(zobs),length(zgrid))
    
    for zz in eachindex(zobs)

        Δz = zgrid .- zobs[zz]

        # last (smallest abs) negative value
        ineg = count(x -> x < 0, Δz)
        
        # first (smallest abs) positive value: ineg + 1

        if iszero(ineg)

            # offscale low obs depth, pick smallest grid depth
            E[zz,begin] = 1
            
        elseif ineg == length(zgrid)

            # offscale high obs depth, pick deepest grid depth
            E[zz,end] = 1
            
        else
            
            wneg = Δz[ineg+1] / (-Δz[ineg] + Δz[ineg+1])
            wpos = -Δz[ineg] / (-Δz[ineg] + Δz[ineg+1])
            E[zz,ineg] = wneg
            E[zz,ineg+1] = wpos
            
        end
    end

    return E
end

"""
    function basinwide_avg(T,Rqq,V,S)
    Short version of arguments
"""
function basinwide_avg(T,Rqq,V,S)
    
    # solve for T̄.
    iRqqV = Rqq\V;
    inside = iRqqV'*V + S\I
    rhs = iRqqV'*T
    T̄ = inside\rhs

    Ctt = inv(inside)
    σT̄ = zeros(length(T̄))
    for ii in eachindex(T̄)
        σT̄[ii] = sqrt(Ctt[ii,ii]);
    end
    
    return T̄,σT̄

end

"""
    function basinwide_avg(T,Rqq,V,S)
    Long version of arguments
    T, H, locs, zgrid are fixed
    params is a dictionary of variable parameters
"""
function basinwide_avg(params)

    @unpack delta, σobs, tratio, sratio, LxyT, LzT, LxyS, LzS, σS, LzAVG, zgrid, zstar = params
    
    T,locs = read_historical_data(delta)

    zobs = zeros(length(locs))
    [zobs[r] = locs[r].depth for r in eachindex(locs)]

    # make H matrix, vertical map onto grid
    V = vertical_map(zobs,zgrid)

    Rqq = error_covariance(locs,σobs,tratio,sratio,LxyT,LzT, LxyS, LzS)
    S = vertical_smoothness(zgrid,σS,LzAVG)

        # solve for T̄.
    iRqqV = Rqq\V;
    inside = iRqqV'*V + S\I
    rhs = iRqqV'*T
    T̄ = inside\rhs

    CT̄ = inv(inside)
    σT̄ = zeros(length(T̄))
    for ii in eachindex(T̄)
        σT̄[ii] = sqrt(CT̄[ii,ii]);
    end

    # compute Heat Content down to depth 
    Km =  degree_meters(T̄,zgrid,zstar)

    # Indian Ocean area divided by 10^21.
    # 15% of global ocean area
    indian_area_m21 = 0.15*0.71*4π*6300_000^2*1e-21

    factor = 3996.0 * 1035.0 * indian_area_m21 # ZJ/(Km)
    # heat content
    H = factor * Km

    # uncertainty in heat content
    σH = 0.0
    nbins = size(CT̄,2)
    tmp = zeros(nbins)
    for nn = 1:nbins
        tmp[nn] = factor*degree_meters(CT̄[:,nn],zgrid,zstar)
    end
    σH = √(factor*degree_meters(tmp,zgrid,zstar))
    
    output = copy(params)
    output["T̄"]= T̄
    output["σT̄"] = σT̄
    output["H"] = H # ZJ
    output["σH"] = σH # error covariance of basinwide avg
    return output
    
end

"""
   Read data regarding historical Indian Ocean cruises
"""
function read_historical_data(delta)
    
    # HistoricalIndianOcean.nc: data from Gazelle, etc.
    # https://drive.google.com/file/d/1-mOo6dwHVwv0TJMoFiYR5QWxykJhNwWK/view?usp=sharing
    # file_id = "1-mOo6dwHVwv0TJMoFiYR5QWxykJhNwWK" old file ID
    file_id = "1_D3XYORYj3J5Ud6Ul6dbu1XjKQnTzEAa"
    url = "https://docs.google.com/uc?export=download&id="*file_id
    filename = datadir("HistoricalIndianOcean.nc")

    # make data directory
    !isdir(datadir()) && mkdir(datadir())

    # download if not already done
    !isfile(filename) &&  google_download(url,datadir())

    nc = NCDataset(filename)

    Tname = "delta_"*delta
    igood = findall(!ismissing,nc[Tname])
    nobs = count(!ismissing,nc[Tname])

    locs = Vector{Loc}(undef,nobs)

    # lon, lat, depth
    #[locs[r] = Loc(nc["lon"][r],nc["lat"][r],nc["depth"][r]) for r in eachindex(nc["delta_T"])]
    [locs[r] = Loc(nc["lon"][igood[r]],nc["lat"][igood[r]],nc["depth"][igood[r]]) for r in eachindex(igood)]

    #################################
    # not good because type includes missing
    #ΔT = nc["delta_T"][igood]

    # here ΔT is type vector float (no missing)
    ΔT = zeros(nobs)
    for ii in eachindex(igood)
        ΔT[ii] = nc[Tname][igood[ii]]
    end

    return ΔT, locs

end

"""
    function degree_meters(T,z,zstar)
"""
function degree_meters(T,z,zstar)

    # layer thickness
    # Δ::Vector{Float64} = z[2:end] - z[1:end-1]
    Km = 0.0
    nbins = length(T) - 1
    # layer fraction

    for i in 1:nbins
        Δ = z[i+1] - z[i]
        if z[i] < zstar < z[i+1]
            δ = (zstar - z[i])/Δ
        elseif zstar ≥ z[i+1]
            δ = 1.0
        elseif zstar ≤ z[i]
            δ = 0.0
        end

        Km += (T[i] * (2-δ) + T[i+1] * δ) * δ * Δ / 2 
    end
    
    return Km
end

"""
function indianarea(latboundary,γ)

calculate Indian Ocean area [m²] given a southern latitudinal boundary
"""
function indianarea(latboundary,TMIversion ="modern_90x45x33_GH10_GH12")
    
    A, Alu, γ, TMIfile, L, B = config_from_nc(TMIversion)
    area = cellarea(γ)

    # now read indian ocean d values.
    bindian = zerosurfaceboundary(γ)

    list = ("ADEL","SUBANTIND","TROPIND")
    for region in list
        bindian += TMI.surfaceregion(TMIversion,region,γ)
    end

    # make a mask for the latitudinal boundary
    for i in eachindex(γ.lon)
        for j in eachindex(γ.lat)
            if γ.lat[j] < latboundary && !isnan(bindian.tracer[i,j])
                bindian.tracer[i,j] =  0.0
            end
        end
    end

    # now multiply mask by area
    return sum((bindian.tracer .* area)[γ.wet[:,:,1]])
end

end # module HistoricalIndianOcean
