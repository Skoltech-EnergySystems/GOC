# GOC
Grid Optimization Competition_Challenge 1.
To be able work with code you can install "JuliaPro Personal" from https://juliacomputing.com/products/juliapro.html. After installation the most optimal is to use "Juno for JuliaPro".
The input data is presented Matlab compatible pscopf.m, generator.csv, contingency.csv, powersystem.raw. 
.m file lacks column headings, so the meaning of values were derived from .raw file. RAW file can be opened with PSS 34 Xplore program.
Now the input operation from .m and .csv files is fully automated. As a result two tables: BLGS (Branch-Load-Generator-Shunt) and Br (Branch) are created. They consist all necessary information for calculations in proper format.
The current problem with the code is that it solves to infeasible solution.
