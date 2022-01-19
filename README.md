# HistoricalIndianOcean

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> HistoricalIndianOcean

It is authored by G Jake Gebbie.

# Code structure

The script at `scripts/GazelleValdiviaPlanet.jl`
solves for the Indian Ocean average temperature change from the time of three historical ocean voyages (from the Gazelle, Valdivia, and Planet, respectively) up to the most recent climatology from the World Ocean Atlas. The scripts runs 8 sets of parameters (and creates 8 solutions) for temperature change and its standard error as a function of depth. 

Basinwide-temperature change is solved using code in `src/HistoricalIndianOcean.jl` that represents an inverse method that accounts for three distinct sources of contamination in temperature observations. The list of available functions is stated in the `export` statement at the top of this source code.

The most "objective" choice of parameters is probably:

`Tbar_LxyS=2000000` m \\
`LxyT=450000` m \\
`LzAVG=500` m \\
`LzS=1000` m \\
`LzT=450` m \\
`sratio=0.04` \\
`tratio=2.0` \\
 `σS=1.0` degrees C \\
 `σobs=0.14` degrees C .

All error bars are 1 sigma. Unicode symbols are used in some places; please submit a GitHub Issue if they cause problems.

# Code output

Output of the code is not included with this GitHub repository. Figures are created in a `plots` directory. CSV output for vertical profile quantities is created in `data/csv`. The complete state of the analysis is saved in Julia output at `data/jld2`. Output in NetCDF or Excel formats could also be created with a little extra work.

# How to run the code

Install the latest stable version of Julia (currently 1.7+). Figure are made with PyPlot which requires Julia to have PyCall installed. Instructions are available at https://github.com/ggebbie/TMI.jl for this step. The code can be run at the command line by changing to the repository base directory and running:

`julia scripts/GazelleValdiviaPlanet.jl` .

For an interactive session, it is possible to run the lines of `GazelleValdiviaPlanet.jl` line by line. The statement `include("intro.jl")` will automatically activate the project environment and download and install the proper Julia packages. 

# Additional steps

1. explanation of all parameters
2. copy and update relevant parts of Gebbie & Huybers 2019 supplementary material into GitHub.
3. check that bin averages reproduce your previous results.
4. Figures are rudimentary. With the csv output, perhaps it is easiest for you to make them publication ready.

# Reproducibility

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
