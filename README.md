
#The files in this repo are for supporting the review of manuscript submitted to J. Appl. Electrochem

### The files in the folder ecm\_model contains equivalent circuit fitting model and the corresponding parameters for each experimental case. This fit is done using zview software. It can be seen that the zview software gives only a few number of data points in its trial version
### The files in the folder exp\_data contains the actual experimental data points of EIS measurements. 
### The file eis\_gsa.jl is the julia code for simulating the equivalent circuit model using the parameters estimated using zview. This simulation can be done by calling the function ecm\_fit(). There are two input files are required for this simulation. "ecm.csv" is a file that specifies the parameter values and "freqdata.csv" specifying the maximum and minimum frequency limit and the number of points to be considered per decade. The function eis\_gsa() performs the global sensitivity analysis using Sobol sampling. 

