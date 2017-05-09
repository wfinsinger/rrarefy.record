This R function resamples randomly and without replacement counts (for instance pollen counts), and calculates rarefied richness (N0), diversity (N1, N2) and evenness (N2/N0) (after Hill (1973) and Birks et al. (2016)).

For additional details, please refer to Finsinger et al. (2017), or contact me directly.

----------------------------------------------------------------------------------------
The R function requires a data.frame and input parameters:
---
One data.frame formatted as:
--
  - Column 1: sample depths
  - Column 2: sample ages
  - Columns 3-n: counts for n taxa

Input parameters:
--
  - counts.file = name of the data.frame
  - site.name = site name
  - zones = zone boundaries. By default zones=NULL
  - n.random = number of samples taken randomly and without replacement. By default n.random=1000
  
Output:
---
the function creates an output directory in the working directory. The new directory will contain:
  - a .csv file with the richness, diversity and evenness estimates,
  - a pdf file with an output figure.

Dependencies
---
The function uses the rrarefy() function from the 'vegan' package (Oksanen et al. 2016) to draw random samples without replacement.

----------------------------------------------------------------------------------------
References:
---
Birks, H.J.B., Felde, V.A., Bjune, A.E., Grytnes, J.-A., Seppä, H., Giesecke, T., 2016. Does pollen-assemblage richness reflect       floristic richness? A review of recent developments and future challenges. Review of Palaeobotany and Palynology 228, 1–25. http://dx.doi.org/10.1016/j.revpalbo.2015.12.011

Finsinger, W., Morales-Molino, C., Gałka, M., Valsecchi, V., Bojovic, S., Tinner, W., 2017. Holocene vegetation and fire dynamics at Crveni Potok, a small mire in the Dinaric Alps (Tara National Park, Serbia). Quaternary Science Reviews 167, 63–77. http://dx.doi.org/10.1016/j.quascirev.2017.04.032

Hill, M.O., 1973. Diversity and Evenness: A Unifying Notation and Its Consequences. Ecology 54, 427–432.

Oksanen, J., Blanchet, G.B., Friendly, M., Kindt, R., Legendre, P., McGlinn, D., Minchin, P.R., O’Hara, R.B., Simpson, G.L., Solymos, P., Stevens, M.H.H., Szoecs, E., Wagner, H., 2016. vegan: Community Ecology Package. R package version 2.4-1.
