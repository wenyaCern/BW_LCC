## Blastwave model + Local charge conservation (Only for Pion , used to try the dragon model with PID, cannot add LCC into this model sucessfully. Need much more tests.)

- Add signal : CMW signal and CME siganl
Switch spatial point correponding to shematic of CMW/CME. Tune the fCMW/fCME = 0 to close the signal part. Note: adding signal needs 
to change setup of parameters.

- Run the code : 
The code can be run directly with "root simpleBW.C" and switch "calcOnline=true" --> to calculate the final results/observables "online".
Also the code can be run on the PC farm with several sub-jobs. Turn calcOnline equals false, and store all of the TProfiles. The final
results can be obtained with the coresponding code "CalcOffline.C". 

- To submite jobs on PC farm ... 
