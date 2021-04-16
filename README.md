# ZikaKinetics
Filenames

Files with names containing HighMOI/one-step/sigle-cycle/OS are for one-step infection kinetics, while LowMOI/multi-step/multi-cycle/MS are for infection started by fewer infectous virus and propagate for multiple cycles. African/VR stands for data of the Afrcican strain Zika virus (ATCC VR1838), and Asian/PR stands for the Puerto Rico isolate that belongs to the Asian lineage of Zika virus. Exponential/ED means the file is about the exponetial-distributed delay models, gamma/GD represnets the gamma-distributed delay models, and Fixed/FD stands for fixed delay models.

Running the scripts

Packages deSolve and minpack.lm need are necessary for running the scripts. Excute the following lines before running the scripts.

  install.packages('deSolve')
  
  install.packages('minpack.lm')
  
The directory of the data files is set by the variable DVdir, and can be customized. The sripts will read data files from the directory and write output to it.
To achieve the minimal pretion error, usually multiple starting values of parameters need to be tested. This can be done by for loops, and the scripts also provide the option of defining the starting values by an input and running the srcipts in batch.

Raw data

Degradation.csv
