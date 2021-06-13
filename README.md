# Hybrid_beamforming_mmWave_Massive_MIMO

These codes are written in MATLAB and are simulated for SIC (Successive Interference Cancelation) based hybrid precoding 
using mmWave channel in 3D scenario where both azimuth and elevation angles are taken into account 
and antennas at the Tx and Rx are arranged in 2D i.e. UPA or URA, in other words full-dimentional
massive MIMO is considered here. These codes are simulated for multi-user and multi-cell cases and also 
compared their performance with different parameters. ALso codes compares different factorization methods
(SVD, EVD, GMD) for capacity optimization. 

Below described the functions/operations performed by different files -
1) SIC_based_HP_3D_MU_FC_SC.m - It simulates SIC based hybrid precoding and optimal precoding schemes
using sub connected and fully connected structures separately for multi-user case in 3D scenario.
2) SIC_based_HP_3D_MU_FC_SC_comparisions.m - In this code the achievable rate of different schemes 
i.e. SIC based hybrid precoding and optimal precoding schemes using sub connected and fully connected 
structures separately for multi-user case in 3D scenario, are compared with the different parameters like 
number of Tx/Rx antennas, number of RF chains at Tx/Rx, number of paths/rays and with the number of users. 
3) SIC_based_HP_multicell.m - Simulates SIC based hybrid precoding and optimal precoding schemes using sub 
connected and fully connected structures separately for multi-cell case in 3D scenario.
4) SIC_based_HP_multicell_comparisions.m - Here the achievable rate of different schemes i.e. SIC based 
hybrid precoding and optimal precoding schemes using sub connected and fully connected structures separately 
for multi-cell case in 3D scenario, are compared with the different parameters like number of Tx/Rx antennas, 
number of RF chains at Tx/Rx, number of paths/rays, number of users per cell and with the number of cells. 
5) SIC_based_HP_3D_MU_FC_SC_using_EVD.m - Compares between SVD and EVD for capacity optimization for 
multi-user case in 3D scenario, where capacity optimization is carried on the SIC based and optimal 
precoding algorithms using sub connected and fully connected structures separately.
6) test_2.m - It is a simulation code for GMD based hybrid precoding using mmWave channel in 2D scenario.
which uses OMP algorithm to find analog precoder and uses GMD for calculating digital precoder.
Acually this is not a complete code, whose work is currently in progress (Code will be updated when 
simulation completes). 
