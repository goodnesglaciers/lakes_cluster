# lakes_cluster (_Stevens et al.,_ In progress)
Recreate figures and movies in Stevens et al. (202?): _Supraglacial lake hydro-fracture not advanced inland by lower-elevation lake drainages in Kalaallit Nunaat/Greenland._ 

## System requirements
This code requires MATLAB R2023B to run. To install, download this code repository to your computer (~150 MB). Also download the FigShare data deposition corresponding to this manuscript (~33 GB). A typical install time on a "normal" desktop computer should be <20 minutes, depending on internet speed. 

Next, follow the instructions in each subsection below for running data-processing and figure-creation scripts. The expected output are the paper figures and reproductions of all quantitative results presented in the manuscript. A typical run time on a "normal" desktop computer will be a few minutes for each figure-producing script. The **nevis** subglacial-hydrology model (_Hewitt,_ 2013) takes ~12 hours to run over a full melt season.

## Probabilistic temporal-cluster analysis
Within [homogeneous_poisson](homogeneous_poisson/), scripts for counting up clusters of lake-drainage events, organized by drainage mechanism: 
    
+ hydro-fracture events (Fig. 3); 
+ moulin events (Fig. 4a–d); and 
+ overspill events (Fig. 4e–h).
  
Daily runoff accumulated at ice-sheet surface elevations within study region of interest courtesy of _Noël et al._ (2019). 

## GNSS-derived quantities
Within [GNSS_derived](GNSS_derived/), scripts for figures containing GNSS-observed and/or GNSS-derived estimates of:

+ station horizontal and vertical positions (Supplementary Information);
+ station horizontal velocities (Supplementary Information);
+ station bed-separation rates (Supplementary Information);
+ inter-station strain rates (Figs. 5 and 6; Supplementary Information); and,
+ map of lake-drainage mechanisms, lake-to-lake supraglacial connections, and GNSS stations (Fig. 1).

## Physical plausibility of hydro-fracture event clusters
Within [physical_plausibility](physical_plausibility/), scripts for physical plausibility of hydro-fracture event clusters C1–C6 (Figs. 7 and 8; Supplementary Information). If the user prefers not to run the Network Inversion Filter (NIF) code, three additional .mat files are needed from the _Stevens et al._ (2024) [lake_strain](https://github.com/goodnesglaciers/lake_strain) repository; alternatively, these large (>25 MB) .mat files may be downloaded on Dropbox [here](https://www.dropbox.com/scl/fo/bbcwl97e8qg0yb9jehusz/h?rlkey=dfclznsaqxzcxtidm9w6gim5j&dl=0).

## Subglacial-hydrology model
Within [nevis_lakes_cluster](nevis_lakes_cluster/), model-run and plotting scripts for the **nevis** subglacial-hydrology model (_Hewitt,_ 2013), which is equivalent to the model version used in _Stevens et al._ (2022) [nevis_helheim](https://github.com/goodnesglaciers/nevis_helheim), save for parameter-value choices and the Central-West-Greenland model domain.  

The model is forced by estimated rates of daily runoff courtesy of _Noël et al._ (2019). Model domain requires BedMachine Greenland v.5 (_Morlighem et al.,_ 2017; 2022). Ice-sheet basal velocities of the model domain are set to surface-velocity values observed by the 2022 MEaSUREs Annual Velocity Mosaic (_Joughin et al.,_ 2015). 

+ [nevis](nevis_lakes_cluster/nevis/): model code.
+ [nevis_lakesix_noSK_300m_ub](nevis_lakes_cluster/nevis_lakesix_noSK_300m_ub/): daily output files for 2022.
+ [nevis_lakesix_2023_noSK_300m_ub](nevis_lakes_cluster/nevis_lakesix_2023_noSK_300m_ub/): daily output files for 2023.
+ Surface runoff, subglacial discharge, and lake-drainage events in 2022 (Supplementary Information Movie M1).
+ Surface runoff, subglacial discharge, and lake-drainage events in 2023 (Supplementary Information Movie M1).

### Example subglacial-hydrology model output alongside mechanistic lake-drainage catalogue: ###
(The mechanistic lake-drainage catalogues are located in these two repositories: [mechanistic_drainage_catalogue_2022](https://github.com/goodnesglaciers/mechanistic_drainage_catalogue_2022) and [mechanistic_drainage_catalogue_2023](https://github.com/goodnesglaciers/mechanistic_drainage_catalogue_2023). The lake-drainage catalogues include figures of every image available for each lake identified by FASTER (_Williams et al.,_ 2018a; 2018b).)

![classifier_2022_HFpossible_E_q_Q_250424_214](https://github.com/user-attachments/assets/d76372ab-fe37-482e-ba6f-2b9fdb42dee0)
Repository Figure 1. Surface runoff, subglacial discharge, and lake-drainage events on 2022/214. (**upper panel**) (blue colormap) Daily surface runoff $E$ for study region courtesy of _Noël et al._ (2019). Day of year of movie frame plotted as map title. (**middle panel**) (blue colormap) Modelled subglacial discharge $q$ forced by daily, distributed surface runoff inputs shown in upper panel. In both map-view panels, ice-sheet surface elevation shown with 100-m-elevation contours in grey (_Morlighem et al.,_ 2017; 2022). GNSS stations shown with solid black triangles. Black-outlined symbols show location and drainage mechanism of supraglacial lakes, with the symbol colour indicating whether the lake is (blue) filling; (goldenrod) initially draining; (green) continuing to drain; or (grey) an empty, dry, or frozen lake basin. Lake-drainage mechanisms are: (stars) hydro-fracture, (overturned triangles) moulin, (circles) overspill, and (square) no-exit, frozen. The location and drainage timing of (diamonds) water-filled crevasses are also shown. Lakes and water-filled crevasses are first plotted in time on the day of year in which they attain $v_{crit}$, their critical volume required to hydro-fracture to the ice-sheet bed. (**bottom panel**) (solid-blue line) Surface-runoff and (dashed-blue line) basal-melt inputs summed across the model domain; proglacial discharge exiting the model domain shown in yellow. Purple vertical bar tracks time. 

## References

Hewitt, I. J. Seasonal changes in ice sheet motion due to melt water lubrication. Earth and Planetary Science Letters 371–372, 16–25 (2013).

Joughin, I., Smith, B., Howat, I. & Scambos, T. MEaSUREs Greenland Ice Sheet Velocity Map from InSAR Data, Version 2. NASA National Snow and Ice Data Center Distributed Active Archive Center https://doi.org/10.5067/OC7B04ZM9G6Q (2015).

Noël, B., Van De Berg, W. J., Lhermitte, S. & Van Den Broeke, M. R. Rapid ablation zone expansion amplifies north Greenland mass loss. Sci. Adv. 5, eaaw0123 (2019).

Morlighem, M. et al. BedMachine v3: Complete Bed Topography and Ocean Bathymetry Mapping of Greenland From Multibeam Echo Sounding Combined With Mass Conservation. Geophysical Research Letters 44, (2017).

Morlighem, M. et al. IceBridge BedMachine Greenland, Version 5. NASA National Snow and Ice Data Center Distributed Active Archive Center https://doi.org/10.5067/GMEVBWFLWA7X (2022).

Stevens, L. A. et al. Tidewater-glacier response to supraglacial lake drainage. Nat Commun 13, 6065 (2022).

Stevens, L. A. Tidewater-glacier response to supraglacial lake drainage (v1.0). Zenodo. https://doi.org/10.5281/zenodo.7023662 (2022).

Stevens, L. A. et al. Elastic Stress Coupling Between Supraglacial Lakes. JGR Earth Surface 129, e2023JF007481 (2024).

Stevens, L. A., & S. Larochelle. Elastic stress coupling between supraglacial lakes (v1.2). Zenodo. https://doi.org/10.5281/zenodo.10650188 (2024).

Williamson, A. G., Banwell, A. F., Willis, I. C. & Arnold, N. S. Dual-satellite (Sentinel-2 and Landsat 8) remote sensing of supraglacial lakes in Greenland. The Cryosphere 12, 3045–3065 (2018a).

Williamson, A. Full source code for the Fully Automated Supraglacial lake Tracking at Enhanced Resolution ("FASTER") algorithm. Apollo - University of Cambridge Repository. https://doi.org/10.17863/CAM.25769 (2018b).

## Correspondence 
Have questions? Please address correspondence to L.A.S. (laura.stevens@earth.ox.ac.uk).
