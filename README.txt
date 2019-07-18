This script uses the MSnbase package (https://github.com/lgatto/MSnbase) to find specific MS2 scans in an LC-MS/MS file in the mzXML format, to export them in order to produce MS/MS databases. This script was written to avoid tedious manual selection of scans when building databases with hundreds of LC-MS runs. 
As it uses MSnbase, cite the following article in case of publication : Gatto L, Lilley KS. MSnbase - an R/Bioconductor package for isobaric tagged mass spectrometry data visualization, processing and quantitation. Bioinformatics. 2012 Jan 15;28(2):288-9. doi:10.1093/bioinformatics/btr645. PubMed PMID:22113085.

##################################################
# Variables :
MasterList : contains the names of the files to be analysed in the first column, along with the masses to be searched in the next columns : adducts, from the most likely ie +/-H to the less likely, 2M-2H+TFA for example. For now, only the first adduct found will be processed. Names of the columns should be explicit as to which adduct is being considered and without special characters (ie "H", "Na") because these will be used in the naming of the output files (see adductsFull variable). An example is present in the repository along with five mzXML files.
adductsFull : names of the columns in MasterList, referring to the names of the adducts being processed. These names will be used for the names of the output files
topN : number of spectra to be merged (from the most intense). Can be set to 1 if only the first scan should be exported, or Inf if all scans are to me merged. If you are using multiple collision energies, topN should be a multiple of the number of collision energies use so as to have the same number of spectra for each energy (3 collision energies : topN =3, 6, 9, 12, 15 etc...). Filtering by collision energy is not yet implemented but can be done by using the collisionenergy() function from MSnbase.
logFile : contains data for each file. Was_detected : "y" if any of the adducts could be found and "n" otherwise. Adduct : first adduct found and processed from the MasterList. isomers : number of isomers / different peaks for the same mz found. spectra_merged : number of spectra used to produce the consensus/merged spectrum, outFiles : names of the output files (same as the input but with the adduct and the RT in the case of isomers).
interval : mass error in ppm to search for the precursor m/z (ie MS/MS spectra with that precursor m/z in the precursorMz tag of the mzXML). 15 ppm is recommended for most LC-MS machines, but the display for precursorMz in Waters machines are often very far from the target mz, and we recommend 3000 ppm, which has been used to a good effect (this has nothing to do with the calibration / measure accuracy of the machine).
writeMGF = FALSE #TRUE : export spectra to MGF, FALSE : data is not exported
mzd = 0.02 # Paremeter for the consensusSpectrum function from MSnbase (see consensusSpectrum)
minProp = 0.11 # Paremeter for the consensusSpectrum function from MSnbase (see consensusSpectrum)

# Varialbes for rudimentary peak picking (to differentiate isomers with different RTs) :
noiseLevel : MS2 scan TIC under which an MS2 scan should not be considered. This helps finding different peaks for a same mass
peakDist : time in seconds between two consecutive MS2 scans, above which they will be considered as belonging to two different peaks
minScans : minimum number of MS2 scans a peak must have to be processed. Set to 1 to pick all detected peaks. Usefull when there are several peaks because of rare MS2 scans with the same precursor mz as the target metabolite (usually 1-2 scans). 

##################################################
Recommended use :
# We recommend doing several exploration runs for your data with tolerant parameters :

writeMGF = FALSE (to avoid writing MGF files during test runs)
topN : Inf (you will be able to see all MS2 scans with the target mz precursor)
noiseLevel : 0 (will take all scans with the targeted mz precursor)
peakDist : Inf (isomers will not be searched)
minScans : 1 (no restriction on the number of scans to be considered)
mzd = 0.02 #consensusSpectrum parameter from MSnbase
minProp = 0.11 #consensusSpectrum parameter from MSnbase

This first step will help you see how many of your files cannot be processed as they probably don't have MS/MS spectra for you molecule (or that it is present as an adduct missing from the MasterList, make sure to check for these files if the run actually is "empty" using a software to visualise your data, like MZmine(http://mzmine.github.io/)). You will also be able to see the total number of MS/MS scans for each ions in the logFile (spectra_merged column).
You can evaluate the quality of the merged scans, as their are all displayed once they've been produced. Check the actual spectra to spot if any peaks are being deleted, or if too much noise is being introduced. You can note the "j" values of each spectrum to make a list of unsatisfactory spectra to rerun them with other parameters. Change the mzd and minProp values to improve the results.

# If you are looking to differentiate peaks with the same precursor m/z, start probing with different parameters while verifying the fitness of the parameters using a software to visualise your data.

topN : 10 (the 10 most intense MS/MS spectra will be merged)
noiseLevel : 1E6 (example of noise level for an orbitrap)
peakDist : 40 (a gap of 40 seconds between two consecutive MS/MS scans will be considered as a gap between two peaks)
minScans : 1 (no restriction on the number of scans to be considered)

Once the correct parameters have been found, set writeMGF to TRUE and run the script to merge and export your spectra to MGF. Rerun it with other parameters on the "j" values of the remaining spectra that gave unsatisfactory results with the first parameters.
