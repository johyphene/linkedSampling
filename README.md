# linkedSampling
These R files are for running a linked emulator following the equations outlined in
Linked Gaussian Process Emulation for Systems of Computer Models Using Matern Kernels and Adaptive Design by Deyu Ming and Serge Guillas.

The files in the burgers folder run the code necessary to build this linked emulator for Burgers equation.
It is necessary to download the data files in the folder along with the R script.

The files in the richards folder run the code necessary to build this linked emulator for Richards equation.
It is necessary to download the data files in the folder along with the R script.
The output soil column graphs in the file show the water content of the 3m tall soil column in cells of 30mm each.

The files in the main folder "GPGraphGenerator.R" and "surfaceSampling.R" run the linked emulator on toy functions.
We then sample from the linked emulator to get spatially correlated samples. 
GPGraphGenerator graphs slices/cross-sections of the surfaces so we can see how things are working at a slice near the middle.
surfaceSampling graphs the entire surface. 
In each case the samples are shown inside of the 95% confidence intervals given by the linked emulator.

Once the data files are downloaded and the paths are adjusted accordingly in the corresponding files, all of the files should run without issue.
If you would like to use and edit the files, you may run into more issues in their current states.
All of the files run with no trend on the second level of emulation which then simplifies many of the calculations in the linked emulator and sampling.
They are not equipped for general cases currently. 
Future updates will focus on making sure any hardcoded dimensions are removed if they persist along with making these general cases rather than just no trend cases.
