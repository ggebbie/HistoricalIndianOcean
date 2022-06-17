# HistoricalIndianOcean

[(https://zenodo.org/badge/doi/10.5281/zenodo.6654779.svg)]

These codes for calculating Indian Ocean basinwide temperature change accompany the paper:

Wenegrat, J.O., E. Bonanno, U. Rack, and G. Gebbie, 2022: A century of observed temperature change in the Indian Ocean. Geophys. Res. Lett., in press, 2022.

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> HistoricalIndianOcean

It is authored by G Jake Gebbie.

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

# Code structure

The script at `scripts/TemperatureDifference_Allparams.jl`
solves for the Indian Ocean average temperature change from the time of two different historical time periods up to the most recent climatology from the World Ocean Atlas. The scripts runs 24 sets of parameters (and creates 24 solutions) for temperature change and its standard error as a function of depth. 

The script at `scripts/TemperatureDifference_Bestparams.jl`
solves for the Indian Ocean average temperature change using the best algorithmic parameters: \
`Tbar_LxyS=2000000` m \
`LxyT=450000` m \
`LzAVG=500` m \
`LzS=1000` m \
`LzT=450` m \
`sratio=0.04` \
`tratio=2.0` \
 `σS=1.0` degrees C \
 `σobs=0.14` degrees C .

3 cases are solved: 
1. from the time of the historical cruises Valdivia, Gazelle, and Planet  up to the most recent climatology from the World Ocean Atlas,
2. from this same time but with the Tait pressure correction applied to historical data,
3. from the time of the historical cruises Valdivia, Gazelle, and Planet  up to the World Ocean Atlas 1955-1964 time interval.

Basinwide-temperature change is solved using code in `src/HistoricalIndianOcean.jl` that represents an inverse method that accounts for three distinct sources of contamination in temperature observations. The list of available functions is stated in the `export` statement at the top of this source code.


All error bars are 1 sigma. Unicode symbols are used in some places; please submit a GitHub Issue if they cause problems.

# Code output

Output of the code is not included with this GitHub repository. Figures are created in a `plots` directory. The `all` directory is output from the sensitivity study using all parameters. The `best` directory contains just the plots from the best parameter choises.  CSV output for vertical profile quantities is created in `data/*/*csv`. The complete state of the analysis is saved in Julia output at `data/*/*jld2`. Output in NetCDF or Excel formats could also be created with a little extra work.

# How to run the code

Install the latest stable version of Julia (currently 1.7+). Figures are made with PyPlot which requires Julia to have PyCall installed. Instructions are available at https://github.com/ggebbie/TMI.jl for this step. Once the package has been instantiated following #1 above, the code can be run at the command line by changing to the repository base directory and running:

`julia scripts/TemperatureDifferenceBestparams.jl` .

For an interactive session, it is possible to run the lines of `TemperatureDifferenceBestparams.jl` line by line. The statement `include("intro.jl")` will automatically activate the project environment and download and install the proper Julia packages. 

# LaTeX documentation

Documentation is found in a LaTeX-generated document in `papers` directory. 

How to set up a LaTeX file for compilation/bibliography:
- open *tex file to enter AucTeX
- `C-c e` to open ebib, should open main.bib, If it doesn't, then toggle LaTeX mode with `M-x LaTeX-mode` (not just latex-mode, not sure why there are two modes with different capitalization)
  - first time creating dependent bib database, use `M c` or similar
    - later times, be sure to `o` open the dependent database in ebib
 - `z` to minimize ebib
   - [associate a database with tex buffer](http://joostkremers.github.io/ebib/ebib-manual.html#associating-a-database-with-a-text-buffer), doesn't seem to work automagically 
     use elisp to set `ebib-local-bibfiles` or `C-h v` to manually customize
    `M-:` then  ` (setq ebib-local-bibfiles '("HistoricalIndianOcean.bib")) `
     `C-c b` will insert new citation from main.bib, add it to the local bib file
     - takes a couple of compiles to get all references/bib correct in pdf file
- [Biblio.el to look up references from internet](https://github.com/cpitclaudel/biblio.el/blob/master/README.md) \
 `M-x biblio-lookup` or `M-x crossref-lookup`
  First open ebib in the main database before doing a crossref lookup

 how to import new citations 
 - `B` in main.bib, give doi.
   or 
  - `M-x biblio-lookup`, choose CrossRef, use `e` to import to ebib
    - must switch back to dependent bib file before citing in tex, use `o`
      - `C-c b` to cite new key\
Otherwise the dependent bib will be overwritten with full main.bib. 

how to compile LaTeX and view PDF 
 - `C-c C-a` ;; i.e., do all
