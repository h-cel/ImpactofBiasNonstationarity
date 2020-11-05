# Impact of bias nonstationarity on the performance of uni- and multivariate bias-adjusting methods

Cite as 

## General overview

This is the code used for the calculations and evaluation in Van de Velde et al. (to be submitted). 
In this paper, 5 bias-adjusting methods are evaluated under climate change circumstances and compared with the change in bias between calibration and validation period. These 5 methods are QDM, mQDM, MBCn, MRQNBC and dOTC.

Note that although this code is publicly available, the observed data used in the paper cannot be shared and has to be requested from the Royal Meteorological Institute in Uccle, Belgium.

## Structure

In this code, you can find the following files. A documentation is included in each file.

Bias adjustment and evaluation:

* a_loadClimateData
* b_configurationBiasAdjustment
  * prepareBiasdata
  * BiasAdjustment
      * occAdj_SSR.m
      * occAdj_TDA.m
        * T
      * occAdj_Threshold.m
      * QDM
      * mQDM
      * MBCn
        * RandMatrix
        * EnDist
      * MRQNBC
	* QDM_all
¨	* coeff
   	* coeffPeriodic
      * dOTC
	* DistEucl
	* NearestFrobenius
	* OptTransPlan
	* OTC 
      * postprocessingSSR
* c_BAEvaluation
  * TruncateObs
  * BA_Evaluation
    * case0
    * CDD
    * Discharge
    * fullfig
    * matload
    * PDMPieter
    * PRCPTOT
    * R10
    * R20
    * RMSE
    * RX1day
    * RX5day
    * SDII
    * SpellDist
    * TransProb
    * val2OyRP

Visualisation:

* Visualisation

Data analysis:

* BiasChangeAll
  * BiasChangeCalc
* ClimaticChanges


## Detailed overview of main files

Six main (or configuration) files are included in this code:
* a_loadClimateData: loads the climate data
* b_configurationBiasCorrection: launches the bias correction/bias adjustment
* c_BCEvaluation: calculates all indices and the RB_O and RB_MB values as used in the article
* Visualisation: uses the RB_O and RB_MB values calculated to make the graphs used in the article
* BiasChangeAll: calculates the R index described in Section 2.4
* ClimaticChanges: makes the plot shown in Section 2.1

Last update on 05/11/'20 by J. Van de Velde
