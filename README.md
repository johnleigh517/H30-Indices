# H30-Indices
H30-Indices, an adapted form of K-indices developed by Yamamazi et al., 2022. Here, an adapted set of python scripts used to derive K-index can also derive H30 and H60-Indices. The original K-indices script was developed by the TCD Solar Physics Group, https://github.com/TCDSolar/K-IndexCalculator (Blake, 2017). These scripts were updated to operation capacity to derive K-Indices for the MagIE website by the DIAS/TCD Solar Physics group (Malone-Leigh, 2025). Here, slight tweeks were made to allow the scripts to generate H30-Indices.


In essence, the adaptation is simple. The original function to calculate K-indices was adapted by:
  a) removing the K = 9 limit for H60 and H30 and replacing it with an unbounded limit following Yamamazi et al. 2022. 
  b) Allowing the K-index module to run at a shorter 60-minute and 30-minute cadence compared to the normal 3-hour
