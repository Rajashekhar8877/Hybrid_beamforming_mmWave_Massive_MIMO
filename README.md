# Hybrid_beamforming_mmWave_Massive_MIMO

MIMO uses 2 to 4 antennas at the transmitter and receiver, which requires individual RF units for each 
antenna. Massive MIMO uses further large number of antennas (at least 64 antennas) to improve throughput and 
spectral efficiency, so it is one of the important technology in 5G. Beamforming is a signal processing 
technique which sends the signal in a particular direction rather than in a broadcast manner, resulting into 
increased dirctivity, energy efficiency, system gain and signal quality. Hence, it is called as one of the key 
technology for 5G and adopted by massive MIMO, but it can support only single data stream with one user. Precoding
is a generalization of beamforming to support multiple data streams which superposes multiple beams. There are 2
types of precoding i.e. analog and digital precoding, whose combination of these 2 gives hybrid precoding which
gives better performance with less complexity. In Successive Interference Cancelation (SIC) based hybrid precoding
method, each antenna array is optimized separately (instead of jointly) such that while optimizing the specific
array the contribution of optimized array is removed from the total capacity. 

These codes are written in MATLAB and are simulated for SIC (Successive Interference Cancelation) based hybrid 
precoding using mmWave channel in 3D scenario where both azimuth and elevation angles are taken into account 
and antennas at the Tx and Rx are arranged in 2D i.e. UPA or URA, in other words full-dimentional
massive MIMO is considered here. These codes are simulated for multi-user and multi-cell cases and also 
compared their performance with different parameters. Also codes compares different factorization methods
(like SVD, EVD, GMD) for capacity optimization in hybrid precoding. 

Below described the functions/operations performed by different files -
1) SIC_based_HP_3D_MU_FC_SC.m - It simulates SIC based hybrid precoding and optimal precoding schemes
using sub connected and fully connected structures separately for multi-user case in 3D scenario.
2) SIC_based_HP_3D_MU_FC_SC_comparisions.m - In this code the achievable rate of different schemes 
i.e. SIC based hybrid precoding and optimal precoding schemes using sub connected and fully connected 
structures separately for multi-user case in 3D scenario, are compared with the different parameters like 
number of Tx/Rx antennas, number of RF chains at Tx/Rx, number of paths/rays and with the number of users.
Here to compare with a specific parameter, you have uncomment that specific section for that parameter in the
code (comment all other sections related to other parameters) and run it.
3) SIC_based_HP_multicell.m - Simulates SIC based hybrid precoding and optimal precoding schemes using sub 
connected and fully connected structures separately for multi-cell case in 3D scenario.
4) SIC_based_HP_multicell_comparisions.m - Here the achievable rate of different schemes i.e. SIC based 
hybrid precoding and optimal precoding schemes using sub connected and fully connected structures separately 
for multi-cell case in 3D scenario, are compared with the different parameters like number of Tx/Rx antennas, 
number of RF chains at Tx/Rx, number of paths/rays, number of users per cell and with the number of cells. 
Here to compare with a specific parameter, you have uncomment that specific section for that parameter in the
code (comment all other sections related to other parameters) and run it.
5) SIC_based_HP_3D_MU_FC_SC_using_EVD.m - Compares between SVD and EVD for capacity optimization for 
multi-user case in 3D scenario, where capacity optimization is carried on the SIC based and optimal 
precoding algorithms using sub connected and fully connected structures separately.
6) test_2.m - It is a simulation code for GMD based hybrid precoding using mmWave channel in 2D scenario.
which uses OMP algorithm to find analog precoder and uses GMD for calculating digital precoder.
Acually this is not a complete code, whose work is currently in progress (Code will be updated when 
simulation completes). 

Note that, in 2D scenario antennas will be arranged in linear i.e. ULA and only azimuth angle is considered
in the channel. 
