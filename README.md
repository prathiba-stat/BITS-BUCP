# BITS-BUCP
This repository contains files that can be used to fit BITS and BUCP models for single case design data
                                                 Instructions
Download all the files on to a given folder. For instance, lets say it is “c:/foo”. I will use this path to demonstrate everything here. Please change this path to whatever folder you create. Note that the path is denoted using a forward slash / and not a backward slash.
The files you will download are:
1.	BITS.R – This program is used to run a Bayesian Interrupted Time-Series model for one single case dataset that is on column 2 of the input data which is in the form of a csv file. In this case that file is called Au_Figure2_SelfBlame.csv. The single case plots are saved with a suffix -SSDplot.jpg, the ROPE of the effect size is saved with a suffix -effect size ROPE.jpg, and the traceplots and histograms from the Bayesian analysis are stored with the suffix -BITS.jpg. All plots are saved in the same folder as the syntax files. The results are saved as a csv file with the suffix -BITSresults.csv.
2.	Au_Figure2_SelfBlame.csv – This CSV file contains 6 columns with the identification of the subject on the top row (header). The first column contains indicators of baseline vs treatment phases (i.e. 0 or 1, respectively), the second column contains the outcome variable of the first person pertaining to these phases in chronological order, the third and the fourth contain data for the second person (phase indicator, outcome variable, respectively), and so on. Thus, this file contains 6 columns of data for 3 subjects, P7, P8, and P9 taken from Au et al. (2007).
3.	BUCP.R – This program is used to run a Bayesian Unknown Change-point model for one single case dataset that is on column 2 of the input data which is in the form of a csv file. In this case that file is called Au_Figure2_SelfBlame.csv. he single case plots are saved with a suffix -SSDplot.jpg, and the traceplots and histograms of the change-point from the Bayesian analysis are stored with the suffix -Change-Point.jpg. All plots are saved in the same folder as the syntax files. The results are saved as a csv file with the suffix -BUCPresults.csv.
4.	BITS-looped.R – This program is an extension of BITS.R that automates running the BITS syntax for all participants.
5.	BUCP-looped.R – This program is an extension of BUCP.R that automates running the BUCP syntax for all participants.
6.	plotPost.R – This is a file that facilitates drawing the ROPE plot. DO NOT CHANGE ANYTHING IN THIS FILE
7.	openGraphSaveGraph.R – This is a file that facilitates drawing the ROPE plot. DO NOT CHANGE ANYTHING IN THIS FILE
8.	HDIofMCMC.R – This is a file that facilitates drawing the ROPE plot. DO NOT CHANGE ANYTHING IN THIS FILE
9.	plots.R – This is a file that facilitates drawing the ROPE plot. DO NOT CHANGE ANYTHING IN THIS FILE
10.	plot_SSD.R – This is a file that facilitates drawing the SSD plot. DO NOT CHANGE ANYTHING IN THIS FILE
11. plot_SSD_line.R - This is a file that facilititates drawing the SSD plot with lines of best fit. DO NOT CHANGE ANYTHING IN THIS FILE

In order to run a syntax, for instance BITS.R, open the file in R or Rstudio. You can initially run the syntax line by line to understand it. All comments are marked with # symbol which means it is a line to explain the syntax that appears below it. 
Some places where you may consider changing parameters to suit your data in BITS.R, BITS-looped.R, BUCP.R, andBUCP-looped.R (using a control+f function) are:
setwd
filename (if you have your own csv file please place that name within quotes)
Set priors according to what you believe would be the scale of the parameter. For instance if you think the mean of the outcome variable may range between 2 and 7, you could use a distribution of dnorm(5, .05). This is the equivalent of a variance of 20 allowing us to be more agnostic about our prior beliefs.



