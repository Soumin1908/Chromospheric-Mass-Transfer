## This page contains the codes and data to reproduce the results in the paper titled "Enhanced Mass Transfer from Extended Shock Supported Stellar Chromospheres" (Bhattacharjee+, 2026)

1) chromosphere_paper - contains the MESA runs. Contains all the inlists, as well as the history.data outputs that are used to generate the plots in the paper. Each folder inside contains one run. 'full_chrom', 'red_chrom' and 'ritter' stand for Full Chromosphere case, Reduced Chromosphere case, and R88 case, respectively (refer paper). Note p___ signifies the period. Note that p2000 is actually the 2250-day period case discussed in the paper.

To generate the MESA plots in the paper, run the 'chromosphere_analyze.ipynb' and 'chromosphere_analyze_new_wind.ipynb' notebooks, in the same directory as uploaded. They produce Figures 5, 6, 7, 10, 11

2) paper_notebooks - Contains the files and codes to generate the remaining figures in the paper -- the mass transfer/wind calculations, as well as the numerical approximation codes. They are commented and explained.

All kinds of comments welcome! Contact: sbhatta2@caltech.edu

Disclaimer: GPT 6.1 Sol have been used at multiple instances of this project, including helping in writing several section of the code, polishing the notebooks, and adding comments. 
