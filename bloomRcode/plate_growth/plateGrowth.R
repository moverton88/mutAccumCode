# 03/19/13
# R Code for detecting quantifying results of plate growth assays from
# Finding the sources of missing heritability in a yeast cross
# http://www.nature.com/nature/journal/v494/n7436/full/nature11867.html
# Joshua Bloom (jbloom@princeton.edu)
# This work is licensed under GPL-2.

# Designed for EBImage version 3.6
# note, does not work with EBImage version 4 due to unresolved issue with rotate function
require(EBImage)
require(gdata)
require(bit)
require(splus2R)
require(doMC)
require(foreach)
require(locfit)
require(reshape)
require(nlme)
require(lme4)
require(pedigree)
require(robustbase)
require(MASS)
registerDoMC(cores=4)

# root directory for data
KSERV='/media/kserver/kruglyak/raid1/home/jbloom'
#KSERV=Sys.getenv('KSERV')

# File and directory setup 

# Make a folder to contain the experiment data.
# Make two subfolders called 'scans' and 'keys'
# In 'scans' put all of your image files (jpg or tiff)
# In 'keys' put all of your key files, including the main key file, the coordinates file, and the strain layout files

#### DESCRIPTION OF KEYS ###################################################################################################


    #####'Key.txt'#################################################################################################
    # a tab-delimited text file with the following 7 columns (column names are not case sensitive):
    #
    # StrainLayout          The filename of a tab-delimited text file with the same dimensions as the plate scanned 
    #                       with strain names as text, blank cells are tolerated and will be converted to 'BLANK'
    #                       by the program.
    # PermutationGroup      An indicator for the the permutation set for each plate, can be any number or string
    #                       Also can be used as an indicator for matched batches
    # Filename              Name of the image file containing plate scanned (short file name, not the full path)
    # PlatePosition         Position on scanner of the plate (order goes top left to bottom left to top right to 
    #                       bottom right A,B,C,(D),(E))
    # Condition             Drug on the plate or some other condition identifier (string or number)
    # Concentration         Concentration on plate or some other identifier (string, number, or blank)
    # Control               [TRUE or FALSE] Flag to indicate that a plate is a control plate 
    ################################################################################################################


    ####'Coordinates.txt'###########################################################################################
    # a tab-delimited text file containing the centers of the corner colonies for each plate
    # the easiest way to make this file is to start ImageJ
    # go to File->Import->Image Sequence
    #       check 'Sort names nnumerically' and 'Use virtual stack' (if you have a lot of images and are RAM)
    # if it is difficult to see the corner colonies go to Image->Adjust->Brightness/Contrast, usually increasing
    # brightness and contrast helps
    # select the multipoint tool
    # From the top right to bottom right to top left to bottom left
    # click on the four colonies on the corners of each plate
    # within each plate it doesn't matter what order you pick the corners
    # you do not have to 
    # with the following 4 columns
    #
    # Filename              Name of the image file containing plate scanned (short file name, not the full path)
    # X                     X-coordinate of grid point
    # Y                     Y-coorinate of grid point
    # PlatePosition         Position on scanner of the plate (order goes top left to bottom left to top right to 
    #                       bottom right A,B,C,(D),(E))
    #
    ################################################################################################################

###########################################################################################################################


# parameter descriptions 
# <fx.file>          = full file path to location of coplementary code containing the functions for this program
# <batch>            = string or number indicating batch, (default =1)
# <working.dir>      = full directory path to directory containing keys and scans
# <key.file>         = filename of main key (whatever you've called it) (default = "Key.csv")
# <coordinates.file> = filename of coordinates key (whatever you've called it) (default = "Coordinates.csv")
# <faceup>           = [TRUE or FALSE]  if plates were scanned faceup then TRUE (default TRUE)
# <rep.4>            = [TRUE or FALSE]  if 4 technical replicates of 96 were pinned in a 384 well format then TRUE (default FALSE)
# <inches>           = [TRUE or FALSE]  are the grid coordinates in Inches or Pixels (if inches then TRUE, if pixels then FALSE) (default FALSE)
# <hp>               = [TRUE or FALSE]  if scanner was HP scanner then TRUE otherwise FALSE (current setup is EPSON scanner, so default is FALSE)
# <bit.depth>        = [8,16] bit depth of scan (default 8)
# <border.offset>    = [INTEGER] distance in pixels from center of colonies on edges to edge of plate (default 26)
# <diagnose.grids>   = [TRUE or FALSE] display overlays of grids on plates while running code, for debugging, slow, might not work anymore (default FALSE)
# <saveImg>          = [TRUE or FALSE] save images 
# <platetype>        = [96,384,1536] 

# Parameters start here -----------------------------------------------
fx.file             <- paste(KSERV, '/source/BYxRM_031713/plateGrowthFx.R', sep='/')
source(fx.file)
batch               <- 1

working.dir = paste(KSERV, '1000BYxRM/phenotyping/XW_20130316/', sep='/')

#check these carefully
key.file            <- 'Key.txt'
coordinates.file    <- 'Coordinates.txt'
faceup              <- TRUE
rep.4               <- FALSE
inches              <- TRUE
hp                  <- FALSE
bit.depth           <- 8
border.offset       <- 34
diagnose.grids      <- FALSE
saveImg             <- FALSE
platetype           <- 384
#-------------------------------------------------------------------


images.dir           <-  paste(working.dir, 'scans/', sep='')
keys.dir             <-  paste(working.dir, 'keys/', sep='')
outfile              <-  paste(working.dir, 'Results.bin', sep='')
outfile.processed    <-  paste(working.dir, 'ResultsProcessed.bin', sep='')
outDirectoryTXTfiles <-  paste(working.dir, 'ResultsTXT/',sep='')
phenoData.bin        <-  paste(working.dir, 'phenoData.bin', sep='')
averageSegData.file  <-  paste(working.dir, 'SegData.txt', sep='')

params = list( outfile          = outfile,    
               images.dir       = images.dir,
               keys.dir         = keys.dir,  
               faceup           = faceup,    
               plate.type       = unlist(plate.types[as.character(platetype)]),
               rep.4            = rep.4,     
               inches           = inches,    
               border.offset    = border.offset,
               diagnose.grids   = diagnose.grids,
               hp               = hp,        
               bit.depth        = bit.depth, 
               saveImg          = saveImg,
               outfile.processed= outfile.processed )
#----------------------------------------------------------------------


#----------------------------------------------------------------------
Key = makeKey(params)
AllStrainLayouts = makeStrainLayouts(Key, params) 
AllColonyGrids   = makeColonyGrids(params)
SKey = split(Key, Key$FILENAME)
#----------------------------------------------------------------------

# Multicore version
#Results = foreach( n=1:length(SKey))  %dopar%   runPlateGrowth(n, SKey,AllStrainLayouts, AllColonyGrids, params)

Results = list()
for(i in 1:length(SKey)) {
    Results[[i]]=runPlateGrowth(i, SKey,AllStrainLayouts, AllColonyGrids, params)
}

Results = unlistify(Results, batch)
save(Results, file=outfile)
Results.Processed=processResults(Results)
save(Results.Processed, file=outfile.processed)
writeTXTfiles(outDirectoryTXTfiles, Results.Processed)

# End core functions  --------------------------------------------



















# Various Filters and post-processing starts here




# post-processing starts here
all.pxs=sapply(Results, function(x){x[,'g.effr']})

# remove plate if top 5% quantile of data has radius <9
bs =apply(all.pxs, 2, function(x) {quantile(x, .95, na.rm=T, names=FALSE) } )
Results[which(bs<9)]=NULL

# Filter on features 
for( i in 1:length(Results) ){
    pa = Results[[i]]
    badstuff = (pa[,'g.edge']>17) | (pa[,'m.pxs']>3600) |  (pa[,'g.sf']>1 & pa[,'g.effr']>20) 
    # default |     (pa[,'g.pdm']>33.5 ) | (pa[,'g.p']>195)  |  (pa[,'g.pdsd']>3.3 ) |  (pa[,'g.sf']>2) | (pa[,'g.acirc']>.4)
    #badstuff = (pa[,'g.edge']>30) | (pa[,'m.pxs']>3900) 
    #|  (pa[,'g.sf']>1 & pa[,'g.effr']>20) 
    # (pa[,'g.pdm']>33.5 ) | (pa[,'g.p']>195)  |  (pa[,'g.pdsd']>3.3 ) |  (pa[,'g.sf']>2) | (pa[,'g.acirc']>.4)
    Results[[i]][badstuff,-c(3,4)] = NA
    # no growth on control plate is NA
}

Na.cnt =sapply(Results, function(x) { sum(is.na(x[,'g.effr'])) })
Results[Na.cnt>57]=NULL

Results.annot    = getAnnotation(Results)
for( i in 1:length(Results) ){
     pa = Results[[i]]
     if(Results.annot[i, 'control']=='TRUE' ) { Results[[i]][pa[,'g.effr']<11,-c(3,4)]=NA}
}

Results.Processed=processResults(Results)

col.row = expand.grid(seq(1,24), seq(1,16))
edges = list(    
    top.row = col.row[,2]==1,
    bottom.row = col.row[,2]==16,
    left.col = col.row[,1]==1,
    right.col =col.row[,1]==24)

for( i in 1:length(Results.Processed) ){
    pa = Results.Processed[[i]]
    pagn = pa[,'g.effr.norm']
    e.ps= sapply(edges, function(x) { tryCatch( {t.test(pagn[x], pagn)$p.value}, error=function(e) {return(1) } )  } )
    for(j in 1:4){   if(e.ps[j]<.05) {Results.Processed[[i]][edges[[j]], -c(3,4)]=NA;   } } }

Results.annot    = getAnnotation(Results.Processed)
Results.controls = Results.Processed[Results.annot$control=='TRUE']
Results.exp      = Results.Processed[Results.annot$control=='FALSE']

R.c.a = getAnnotation(Results.controls)
R.e.a = getAnnotation(Results.exp)

R.c.YPD = Results.controls[!grepl('YNB|ynb', R.c.a$condition)]
R.e.YPD = Results.exp[!grepl('YNB|ynb', R.e.a$condition)]
R.c.YPD.a = getAnnotation(R.c.YPD)
R.e.YPD.a = getAnnotation(R.e.YPD)

R.C.P.s =split(R.c.YPD, paste(R.c.YPD.a$batch, R.c.YPD.a$strainlayout,sep='_'))
R.E.P.s =split(R.e.YPD, paste(R.e.YPD.a$batch, R.e.YPD.a$strainlayout, sep='_'))

YPD.avg=lapply(R.C.P.s, avgPlateControls)
R.E.P.s=R.E.P.s[names(R.E.P.s) %in% names(YPD.avg)]
YPD.avg = YPD.avg[names(YPD.avg) %in% names(R.E.P.s)]

YPD.norm = mapply( makedtf, R.E.P.s, YPD.avg, SIMPLIFY=FALSE)

data1=unlistify2(YPD.norm)
data3=unlistify2(R.C.P.s)
# to facilitate post-processing repeat g.effr columns for controls with name of 'controls removed'
data3=lapply(data3, makedtf_include_control_data)

sResults = c(data1, data3)
sResults.annot = getAnnotation(sResults)
sResults.annot$concentration = as.numeric(sResults.annot$concentration)

##writeTXTfiles(outDirectoryTXTfiles, sResults)
#
sby  = paste(sResults.annot$batch, sResults.annot$condition, sResults.annot$concentration, sep='_')
dResults = split(sResults, sby)

phenoData= lapply(dResults, groupStrainPhenos, pheno='g.effr.norm.ypdc') 

save(file ='/media/kserver/kruglyak/raid1/home/jbloom/1000BYxRM/phenotyping/ST_JB_Haploid_Phenotyping_090512/pAB_103012.bin', list=c('phenoData.a', 'phenoData.b'))
#phenoData = reducePhenoData(phenoData, strains_used)
#seg.data = getTable(phenoData,'notsegs.mean')
#phenoData.sa = lapply(phenoData, function(x){x$notsegs.all})
#phenoData.a = phenoData
#BroadH2      = sapply(phenoData.sa, calc.BroadH2, jackknife=FALSE)
#write.table(BroadH2, file='/media/kserver/kruglyak/raid1/home/jbloom/1000BYxRM/phenotyping/090612-YPSYJM/BroadH2.txt', sep='\t', quote=F)

seg.data = getTable(phenoData,'notsegs.mean')
phenoData.sa = lapply(phenoData, function(x){x$notsegs.all})

