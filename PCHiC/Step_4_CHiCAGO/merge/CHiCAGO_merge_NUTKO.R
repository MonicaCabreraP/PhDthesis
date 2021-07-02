library(Chicago)
outputDirectory <- ("/gpfs/projects/bsc08/bsc08471/p53/Step2_CHiCAGO/downsampled/NUTWT_recalibration/NUTKO")
weightsPath <- file.path("/gpfs/home/bsc08/bsc08471/softwares/CHiCAGO/weights_recalibration")
settingsFile <- file.path(weightsPath, "NUTWT.settings")
outprefix <- ("NUTKO")
files <- c(
    file.path("/gpfs/projects/bsc08/bsc08471/p53/Step2_CHiCAGO/downsampled/chinputs/NUT1_KO_downsampled/NUT1_KO_downsampled.chinput"),
    file.path("/gpfs/projects/bsc08/bsc08471/p53/Step2_CHiCAGO/downsampled/chinputs/NUT2_KO_downsampled/NUT2_KO_downsampled.chinput"))
DesignDir <- ("/gpfs/home/bsc08/bsc08471/softwares/CHiCAGO/designDir")
cd <- setExperiment(designDir = DesignDir, settingsFile = settingsFile)
cd@settings$brownianNoise.seed <- 103
cd <- readAndMerge(files=files, cd=cd)
cd <- chicagoPipeline(cd, outprefix = outprefix, printMemory = T)
saveRDS(cd, paste0(outputDirectory,"/",outprefix, ".Rds"))
exportResults(cd, file.path(outputDirectory, paste0(outprefix,"_cutoff_0")), cutoff = 0)
exportResults(cd, file.path(outputDirectory, paste0(outprefix,"_cutoff_5")), cutoff = 5)
