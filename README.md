# Radar-Interference-JSTSP
Author: Sian Jin (sianjin@uw.edu)

This document describes MATLAB files for obtaining the figures in the JSTSP (Special Issue on Recent Advances in Automotive Radar Signal Processing) paper: 
"FMCW Radar Network: Multiple Access and Interference Mitigation". 


Prerequisites: MATLAB 2020 

Folders:

1)  PHY: PHY layer simulation scripts for the 4 interference mitigation schemes: UA-RFDM, UA-FH, UA-RFDM-PC, UA-FH-PC.

	 PHY layer performance of UA-RFDM is obtained by the file twoDimFFTUA-RFDM.m.  
	 PHY layer performance of UA-FH is obtained by the file twoDimFFTUA-FH.m.  
	 PHY layer performance of UA-RFDM-PC is obtained by the file twoDimFFTUA-RFDM-PC.m.  
	 PHY layer performance of UA-FH-PC is obtained by the file twoDimFFTUA-FH-PC.m.    
	 
	 The PC related parameter p0 is obtained through UARFDMPCCalculatep0.m and UAFHPCCalculatep0.m for UA-RFDM-PC and UA-FH-PC, respectively.
  
2)  MAC: MAC layer simulation scripts for the 4 interference mitigation schemes: UA-RFDM, UA-FH, UA-RFDM-PC, UA-FH-PC.

	 Multiple access capacity VS density performance of the 4 interference mitigation schemes is obtained by the file capacity.m. 
	 Probability of target misdetection VS density performance of the 4 interference mitigation schemes is obtained by the file pmd.m. 
	 Multiple access capacity VS probability of target misdetection performance of the 4 interference mitigation schemes is obtained by the file capacityVspmd.m. 
  
