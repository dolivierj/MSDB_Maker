library(MSnbase)
setwd("path")
masterList= as.matrix(read.csv("MasterList.csv")) #Parameter file with compound names and m/z for several adducts
adductsFull = colnames(masterList) # Adducts provided by the masterList
topN = 10 #number of spectra to be merged (from the most intense). Set to Inf if all spectra should be merged
interval = 15 #mass error in ppm : 15 for thermo Orbitrap is fine, 3000 for waters QToF
writeMGF = FALSE #TRUE : export spectra to MGF, FALSE : data is not exported
mzd = 0.02 # Paremeter for the consensusSpectrum function from MSnbase
minProp = 0.11 # Paremeter for the consensusSpectrum function from MSnbase
noiseLevel = 0 #useful for isomer differentiation, otherwise a a very low noiseLevel or 0 would do no harm
peakDist = Inf #minimum distance between two MS/MS scans to consider different peaks. set Inf if no isomer differentiation is required
minScans = 1 #set to 1 to avoid restrictions
logFile = matrix(data = NA, nrow = length(masterList[,2]), ncol = 6)
colnames(logFile) = c("files", "Was_detected", "Adduct", "isomers","spectra_merged", "outFiles")
logFile[,1] = masterList[,1]

#voidpeaks = which(logFile[,2] == "n") # Used to review molecules that were not detected
#logFile[voidpeaks] # Used to review molecules that were not detected

#length(masterList[,2])
for (j in 1:length(masterList[,2])) {
  
  #Setting target parameters
  k = 2 #current adduct counter. 2 = H, if nothing is found in H, will sequentially try adducts in the masterList
  file = masterList[[j,1]] #Name of the file selected in the masterList
  targetMass = as.numeric(masterList[[j,k]]) #Expected measured mass for the compound
  
  #Readfile
  msmsData = readMSData(file, pdata = NULL, msLevel. = 2,
                       verbose = isMSnbaseVerbose(), centroided. = NA, smoothed. = NA,
                       cache. = 1L, mode = c("onDisk"))
  
  
  #Filter scans using the noiseLevel
  msmsData = msmsData[as.numeric(which(tic(msmsData) > noiseLevel))]
  
  
  #find the MS/MS spectra with the targeted precursor, switch to next adduct if nothing found. If absolutely no spectra are found, skip to next compound
  while ((length(filterPrecursorMz(msmsData, targetMass, ppm = interval)) == 0) && (k < length(adductsFull))) {
    k = k + 1
    targetMass = as.numeric(masterList[[j,k]])
  }
  if (length(filterPrecursorMz(msmsData, targetMass, ppm = interval)) == 0) {
    logFile[j,2] = "n"
    next
  }else{
    msmsData = filterPrecursorMz(msmsData, targetMass, ppm = interval)
    logFile[j,2] = "y"
    logFile[j,3] = adductsFull[k]
  }
  #Find time gaps between scans superior to peakDist to differentiate isomers / peaks
  rtVector = as.numeric(MSnbase::rtime(msmsData))
  
  if (length(rtVector)==1 && length(rtVector) >= minScans) { #if only one scan is found in the filtered data, no merge, show only the one scan
    outFile = paste0(gsub(".mzXML", "_", file), adductsFull[k], "_single.mgf", collapse = NULL)
    plot(mz(msmsData[[1]]), intensity(msmsData[[1]]), type = "h", main = paste0("j=",j,", ",outFile))
    logFile[j,4] = 0
    logFile[j,5] = 1
    logFile[j,6] = outFile
    if(writeMGF == TRUE) {MSnbase::writeMgfData(msmsData[[1]],con = outFile)}
    #MSnbase::writeMgfData(msmsData[[1]],con = outFile)
    next
  }else{}
  
  rtGaps = as.numeric(diag(as.matrix(dist(rtVector))[2:length(msmsData),]))
  cutlist = which(rtGaps >= peakDist) #put back peakDist
  
  ### HERE : escape route for files without groups / isomers (short route)
  if (!length(cutlist) == 0) {#long route
    #Add here the number of isomers
    logFile[j,4] = length(cutlist)
    groupList = list()
    l = 1
    for (i in 1:length(cutlist)) {
      groupList[[i]] = l:cutlist[i]
      l = cutlist[i] + 1
    }
    groupList[[length(cutlist)+1]] = l:length(rtVector)
    
    #Selects which groups are valid (number of spectra >= minScans)
    selectedGroups = c()
    for (i in 1:length(groupList)) {
      if (length(groupList[[i]]) >= minScans) {
        selectedGroups = c(selectedGroups, i)
      }
    }
    #Merge scans for each group
    nbspectra = ""
    outFileList = ""
    for (i in 1:length(selectedGroups)) {
      data_03 = msmsData[groupList[[selectedGroups[i]]]] #select from msmsData the scans for each group
      if (length(data_03) > topN) { ### Check wether there are more than topN scans. If not, skip as selecting scans is useless
        msmsIntensities = as.numeric(tic(data_03)) #Find the tic of each scan in data_03
        msmslowestInt = min(msmsIntensities[rank(-msmsIntensities) <= topN]) # find the least intense scan among the topN most intense scans
        scans2merge = which(msmsIntensities >= msmslowestInt) #here are the topN most intense spectra
        data_03 = data_03[scans2merge] #reduce data_03 to the topN number of scans
      }else{} #else : there are less than topN scans, no need to reduce data_03
      rtval = as.character(data_03[[which(as.numeric(tic(data_03)) == max(as.numeric(tic(data_03))))[1]]]@rt) #get the RT of the group / isomer
      rtval = gsub("\\.", "s", rtval) #replace the dot in the RT by "s" (seconds)
      #merge the spectra in data_03 : 
      mergedSpectra = consensusSpectrum(Spectra(spectra(data_03), elementMetadata = DataFrame(idx = 1:length(spectra(data_03)))), mzd, minProp)
      nbspectra = paste0(nbspectra,"|", length(data_03))
      outFile = paste0(gsub(".mzXML", "_", file), adductsFull[k], "_", rtval, "_merged.mgf", collapse = NULL)
      outFileList = paste0(outFileList, "|", outFile)
      plot(mz(mergedSpectra), intensity(mergedSpectra), type = "h", main = paste0("j=",j,", ",outFile))
      #write the file for each isomer with the RT to differentiate them
      #Report in the log file the number of isomers detected, their RT and their file name (make an outFile variable)
      if(writeMGF == TRUE){MSnbase::writeMgfData(mergedSpectra,con = outFile)}
      #MSnbase::writeMgfData(mergedSpectra,con = outFile)
    }
    nbspectra = substring(nbspectra, 2, nchar(nbspectra))
    outFileList = substring(outFileList, 2, nchar(outFileList))
    logFile[j,5] = nbspectra
    logFile[j,6] = outFileList
  }else{#short route
    #Select the top N spectra and merge
    if (length(msmsData) < minScans) {
      next
    }else{}
    msmsIntensities = as.numeric(tic(msmsData))
    msmslowestInt = min(msmsIntensities[rank(-msmsIntensities) <= topN])
    scans2merge = which(msmsIntensities >= msmslowestInt) #here are the topN most intense spectra
    msmsData = msmsData[scans2merge]
    mergedSpectra = consensusSpectrum(Spectra(spectra(msmsData), elementMetadata = DataFrame(idx = 1:length(spectra(msmsData)))), mzd, minProp)
    outFile = paste0(gsub(".mzXML", "_", file), adductsFull[k], "_merged.mgf", collapse = NULL)
    plot(mz(mergedSpectra), intensity(mergedSpectra), type = "h", main = paste0("j=",j,", ",outFile))
    logFile[j,4] = 0
    logFile[j,5] = length(msmsData)
    logFile[j,6] = outFile
    if(writeMGF == TRUE) {MSnbase::writeMgfData(mergedSpectra,con = outFile)}
    #MSnbase::writeMgfData(mergedSpectra,con = outFile)
  }
}
