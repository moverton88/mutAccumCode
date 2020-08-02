# 03/19/13
# R script for heritability and QTL analysis in:
# Finding the sources of missing heritability in a yeast cross
# http://www.nature.com/nature/journal/v494/n7436/full/nature11867.html
# Joshua Bloom (jbloom@princeton.edu)
# This work is licensed under GPL-2.

library(lme4)
library(foreach)
library(doMC)
library(rrBLUP)
library(qtl)


registerDoMC(cores=8)

# location of QTL_mappingFx.R
src.directory= '/media/kserver/kruglyak/raid1/home/jbloom/source/BYxRM_031713/'
source((paste(src.directory, "QTL_mappingFx.R", sep="")))

# Loads 'pheno_raw' list (contains phenotype measurements for each trait)
load(url("http://genomics-pubs.princeton.edu/YeastCross_BYxRM/data/pheno_raw.Rdata"))

# Loads R/QTL object 'cross' (contains average phenotype for each segregant for each trait,
# markers for QTL mapping and genetic map
load(url("http://genomics-pubs.princeton.edu/YeastCross_BYxRM/data/cross.Rdata"))



#Calculate broad-sense heritability ####################################################################################
BroadH2      = sapply(pheno_raw, calc.BroadH2, jackknife=FALSE)
# Calculate SE for broad-sense heritability
BroadH2.jkSE = sapply(pheno_raw, calc.BroadH2, jackknife=TRUE)
#########################################################################################################################



#Calculate narrow-sense heritability#####################################################################################
#pdata.01 = mean for each segregant (for Figure 2, same scale as used for QTL mapping)
#pdata.02 = one value for each segregant (for Figure 1, same scale as broad-sense)

gdata     = extractGenotype(cross)
n.pheno   = countStrainsPerTrait(cross$pheno) 
pdata.01  = extractScaledPhenotype(cross, TRUE)

#first element of each replicate that isn't NA
pdata.02 =  getOnePhenoValue(pheno_raw, pdata.01)
n.pheno.02 = apply(pdata.02, 2, function(x) {sum(!is.na(x))})


# use rrBLUP or regress package
# mixed.solve is in rrBLUP ... constraint is that there can only be one random effect term and random effect term covariance matrix
# regress allows multiple random effects term (slower?) but doesn't compute BLUP estimates
# can also use emma.MLE
# Calculate segregant relatedness
A = A.mat(gdata, shrink=FALSE)/2

# narrow-sense heritability for avergage segregant phenotype
narrowH2 = foreach(i=1:ncol(pdata.01)) %dopar% {print(i); mixed.solve(pdata.01[,i], K=A, method='REML') }
narrowH2 = sapply(narrowH2, function(x) {x$Vu/(x$Ve+x$Vu) })

# narrow-sense heritability on same scale as broad-sense heritability
narrowH2.fig1 = foreach(i=1:ncol(pdata.02)) %dopar% {print(i); mixed.solve(pdata.02[,i], K=A, method='REML') }
narrowH2.fig1 = sapply(narrowH2.fig1, function(x) {x$Vu/(x$Ve+x$Vu) })

# get standard error for narrowH2 (repeat for narrowH2.fig1) 
jackh2s = calc.mixed.jack.pseudovalues(pdata.01, A)
narrowH2.se = rep(NA, length(narrowH2))
for (t in 1:length(narrowH2) ) { narrowH2.se[t]=calc.jkSE(n.pheno[t], jackh2s[[t]], narrowH2[t])}
narrowH2.bias = rep(NA, length(narrowH2))
for (t in 1:length(.A) ) { narrowH2.bias[t]=calc.jkbias(n.pheno[t], jackh2s[[t]], NarrowH2[t])}


# calculate narrow h2 per chromosome
mis = getMarkerIndexSplit(cross)
# build relationship matrices for subsets of gdata
chrRelMats = lapply(mis, function(x) { A.mat(gdata[,x])/2 })

nh2.chr = foreach(i=1:ncol(pdata.01)) %dopar% {print(i); 
    lapply(chrRelMats, function(x) {
                mixed.solve(pdata.01[,i], K=x, method='REML')})}
pchr.mat= sapply(nh2.chr, function(y) {sapply(y, function(x) {x$Vu/(x$Ve+x$Vu) } ) } )

# end Narrow-sense heritability code ######################################################################################






# QTL Mapping ################################################################################################# 

mindex.split = getMarkerIndexSplit(cross)
# get chromosome offsets  (to convert from marker index per chromosome to marker index across genome)
chr.mindex.offset = sapply(mindex.split, min)-1

######extract phenotypes, genotypes, and number of individuals phenotyped per trait 
gdata     = extractGenotype(cross)
n.pheno   = countStrainsPerTrait(cross$pheno) 
pdata.01      = extractScaledPhenotype(cross, TRUE)
LODS.01       = get.LOD.by.COR(n.pheno, pdata.01, gdata)
LODS.01s      = LODmatrix.2.scanone(LODS.01, cross)
peaklist.01   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.01) 
#LODS.01.FDR   = getPeakFDR(peaklist.01$chr.peaks.lod, pdata.01, gdata, 1000)
peakArray.01  = getPeakArray(peaklist.01, 2.69)

pdata.02      = getPhenoResids(pdata.01, gdata, peakArray.01) 
LODS.02       = get.LOD.by.COR(n.pheno, pdata.02, gdata)
LODS.02s      = LODmatrix.2.scanone(LODS.02, cross, LODS.01s)
peaklist.02   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.02) 
#LODS.02.FDR   = getPeakFDR(peaklist.02$chr.peaks.lod, pdata.02, gdata, 1000)
peakArray.02  = getPeakArray(peaklist.02, 2.88)
#peakArray.02  = rbind(peakArray.01, peakArray.02)

pdata.03      = getPhenoResids(pdata.02, gdata, peakArray.02) 
LODS.03       = get.LOD.by.COR(n.pheno, pdata.03, gdata)
LODS.03s      = LODmatrix.2.scanone(LODS.03, cross, LODS.01s)
peaklist.03   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.03) 
#LODS.03.FDR   = getPeakFDR(peaklist.03$chr.peaks.lod, pdata.03, gdata, 1000)
peakArray.03  = getPeakArray(peaklist.03, 3.75)
#peakArray.03  = rbind(peakArray.02, peakArray.03)
 
pdata.04      = getPhenoResids(pdata.03, gdata, peakArray.03) 
LODS.04       = get.LOD.by.COR(n.pheno, pdata.04, gdata)
LODS.04s      = LODmatrix.2.scanone(LODS.04, cross, LODS.01s)
peaklist.04   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.04) 
#LODS.04.FDR   = getPeakFDR(peaklist.04$chr.peaks.lod, pdata.04, gdata, 1000)
peakArray.04  = getPeakArray(peaklist.04, 4.63)
#peakArray.04  = rbind(peakArray.03, peakArray.04)

pA1=cbind(peakArray.01, 'J1')
pA2=cbind(peakArray.02, 'J2')
pA3=cbind(peakArray.03, 'J3')
pA4=cbind(peakArray.04, 'J4')

names(pA1)[3]='jump'
names(pA2)[3]='jump'
names(pA3)[3]='jump'
names(pA4)[3]='jump'

peak.index = data.frame(rbind(pA1,pA2,pA3, pA4))
peak.index = peak.index[order(peak.index$trait, peak.index$markerIndex),]
peak.index = split(peak.index, peak.index$trait)

#save(peak.index, file='~/1000BYxRM/QTL/peakindex.bin')
# also ran with refine = TRUE 
#save(fQTLs_FDR05r, file = '~/1000BYxRM/QTL/fQTLs_FDR05r.bin')

fQTLs_FDR05 = foreach( i=1:length(peak.index) ) %dopar% doQTLModel(i,peak.index,cross, LODS.01s, refine=TRUE)
names(fQTLs_FDR05) = names(peak.index)
fQTLs_FDR05 =lapply(fQTLs_FDR05, function(x) {
        x$CIs = data.frame((as.vector(x$fqtl$ests$ests[-1])), x$CIs)
        names(x$CIs)[1]='eff.size'
        return(x)
        })

# variance explained without cross-validation
VarExpQTLsAdditive=sapply(fQTLs_FDR05, function(x) { 
    y=NA
    tryCatch( 
        {
         y= as.numeric(x$fqtl$result.full[,'%var'][1]) },
           error= function(e){y=NA})
   return(y/100)} )
VarExpQTLsAdditive=t(data.frame(t(VarExpQTLsAdditive)))

QTLcnt =sapply(
    fQTLs_FDR05, function(x) {
    if(!is.null(nrow(x$qtlMarkersCP))){ return(nrow(x$qtlMarkersCP)) }
    else{return(NA)}
})
QTLcnt = t(data.frame(t(QTLcnt)))


# variance explained with cross-validation----------------------------------------------------------------------------------
cv10 = doCV(cross, LODS.01s, mindex.split, chr.mindex.offset,cv.fold=10)

# summarize results from cross validation
cv_ve10 = sapply(cv10, function(x) { x$ve })
cv_se10 = apply(cv_ve10, 1, function(x) {sqrt(sum((x-mean(x))^2)) } )
cv_se10 = cv_se10/sqrt(10)

Additive_10xCV.nqtl   = sapply(cv10, function(x){nrow(do.call('rbind',x$peak.index))})
Additive_10xCV.mean   = apply(cv_ve10, 1, mean)
Additive_10xCV.median = apply(cv_ve10, 1, median)
Additive_10xCV.min    = apply(cv_ve10, 1,quantile, .25)
Additive_10xCV.max    = apply(cv_ve10, 1,quantile, .75)
Additive_10xCV.se     = cv_se10
cv10.H2               = data.frame(Additive_10xCV.mean, Additive_10xCV.median, Additive_10xCV.min,  Additive_10xCV.max, Additive_10xCV.se)
rownames(cv10.H2)     = nH2$Trait
rm(cv10)
# end QTL mapping ##############################################################################################################################



# calculate LOD scores (as in Figure 4) ------------------------------------------------------------------------------
# last column in each data frame is the LOD score from composite mapping
cimMOD= foreach (p=1:length(peak.index)) %dopar% {
    print(p)
    so.out = LODS.01s[,c(1,2,p+2)]
    acv = peak.index[[p]]$markerIndex
    if(is.null(acv)){ retun(so.out) }
    else{ 
        josh.cim = sapply(1:length(acv), function(i) {
                    scanone(cross, pheno.col=p, method='mr', addcovar=gdata[,acv[-i]])$lod })
        colnames(josh.cim)=acv
        so.out = data.frame(so.out, josh.cim)
        class(so.out)=c('scanone', 'data.frame')
        josh.cim.max = apply(so.out[-c(1,2)], 1, max)
        so.out$cim=josh.cim.max
           return(so.out)
    }
}
names(cimMOD) = names(fQTLs_FDR05)
save(cimMOD, file= '~/1000BYxRM/QTL/cimMOD.bin')
#--------------------------------------------------------------------------------------------------------------------






# Detecting QTL-QTL  interactions###############################################################################################
library(eqtl)
# setup for full 2D scan on reduced marker set------------------------------------------------------------------------
    
    pm = pull.map(cross)
    allm = as.character(do.call('c', sapply(pm, names)))
    p.1cm = as.character(do.call('c', sapply(pm, pickMarkerSubset, 1)))

    pcross = drop.markers(cross, allm[!(allm %in% p.1cm)])
    pcross = calc.genoprob(pcross)

    pAF =do.call('rbind', peak.index)[,c(1,2)]
    #pAF = rbind(peakArray.01, peakArray.02, peakArray.03, peakArray.04)
    pdata.05      = getPhenoResids(pdata.01, gdata, pAF, intercept=FALSE)
    pcross2       = pcross
    pcross2$pheno = pdata.05
    g.pcross.data = extractGenotype(pcross2)
    # NOTE .... scantwo was run on cluster 


    for(ff in 1:100) {
        pcross3 = pcross2
        pcross3$pheno=pcross3$pheno[sample(1:nrow(pcross3$pheno)),]

        iLODshp = array(0, dim=c(46,ncol(g.pcross.data), ncol(g.pcross.data)))
        iLODsp=foreach(i = 1:(ncol(g.pcross.data)-25) ) %dopar% {
            print(i)
            nullInt = seq(1, i+24)
            gInt = g.pcross.data*g.pcross.data[,i]
            gInt[,nullInt]=0
            get.LOD.by.COR(n.pheno, pcross3$pheno, gInt)
        }
        for(i in 1:(ncol(g.pcross.data)-25)) { print(i); iLODshp[,i,]=iLODsp[[i]]    }
        save(iLODshp, file =paste('~/Desktop/s2/exp', ff, '.bin', sep=''))
    }

    # results from R/qtl scantwo run on cluster
    load('~/Desktop/s2/obs.bin')
    load('~/Desktop/s2/exp_peaks.bin')
    
    iexp.peaks
    iobs.peaks =apply(iLODsh, 1, peakfinder2D, threshold=3)
    names(iobs.peaks)=names(iexp.peaks[[1]])

    tiobs.cnt = sapply(iobs.peaks, function(x) {     sapply(seq(2,7, .1), function(tt) {sum(x[,'lod']>tt)}) })
    rownames(tiobs.cnt)=seq(2,7,.1)
    tiexp.peaks = list()
    for (i in 1:46) { tiexp.peaks[[i]] =sapply(iexp.peaks, function(x) x[[i]][,'lod'])  }
    tiexp.cnt=sapply(tiexp.peaks, function(slist) {
        iecnt= sapply(slist, function(x) {     sapply(seq(2,7, .1), function(tt) {sum(x>tt)}) })
        rownames(iecnt)=seq(2,7,.1)
        apply(iecnt, 1,mean)
    })
    colnames(tiexp.cnt)=names(iobs.peaks)
    full2D.FDR.trait = tiexp.cnt/tiobs.cnt

    fiobs.cnt = apply(tiobs.cnt,1,mean)
    fiexp.cnt=lapply(tiexp.peaks, function(slist) {
        iecnt= sapply(slist, function(x) {     sapply(seq(2,7, .1), function(tt) {sum(x>tt)}) })
        rownames(iecnt)=seq(2,7,.1)
        return(iecnt)
    })

    fiexp.cnt=do.call('cbind', fiexp.cnt)
    fiexp.cnt=apply(fiexp.cnt, 1, mean)

    full2D.FDR= fiexp.cnt/fiobs.cnt
    # FDR 10% = 6.2
    # FDR 5% = 6.5


    g.pcross.data

    iobs.peaks.g = lapply(iobs.peaks, function(z) {
           zx= colnames(g.pcross.data)[ z[,'x'] ]
           zy= colnames(g.pcross.data)[ z[,'y'] ]
           xg = match(zx, colnames(gdata))
           yg = match(zy, colnames(gdata))
           cbind(z, xg,yg)
            })
    save(iobs.peaks.g, file='~/1000BYxRM/QTL/iobspeakg.bin')

    iobs.peaks.sig = sapply(iobs.peaks, function(x) {x[x[,3]>6.2,]})


# transform x and y to marker indices from cross and gdata
#-------------------------------------------------------------------------------------------------------------------------------
#plot(ipeaks[,'x'], ipeaks[,'y'], pch='.', xlim=c(0,2800), ylim=c(0,2800), type='n')
#text(ipeaks[,'x'], ipeaks[,'y'], round(ipeaks[,'lod'],1), pch='.', xlim=c(0,2800), ylim=c(0,2800))
#abline(0,1)

# For marginal 2D scan (as in figure 5) --------------------------------------------------------------------------------
    pm = pull.map(cross)
    allm = as.character(do.call('c', sapply(pm, names)))
    p.1cm = as.character(do.call('c', sapply(pm, pickMarkerSubset, .5)))

    pcross = drop.markers(cross, allm[!(allm %in% p.1cm)])
    pcross = calc.genoprob(pcross)

    pAF =do.call('rbind', peak.index)[,c(1,2)]
    #pAF = rbind(peakArray.01, peakArray.02, peakArray.03, peakArray.04)
    pdata.05      = getPhenoResids(pdata.01, gdata, pAF, intercept=FALSE)
    pcross2       = pcross
    pcross2$pheno = pdata.05
    g.pcross.data = extractGenotype(pcross2)



imindex.split = getMarkerIndexSplit(pcross2)
# get chromosome offsets  (to convert from marker index per chromosome to marker index across genome)
ichr.mindex.offset = sapply(imindex.split, min)-1

intScans =scanIntwMain(fQTLs_FDR05, n.pheno, pcross2, g.pcross.data, imindex.split, ichr.mindex.offset, doFDR, doGWER)

#save(intScans, file='~/1000BYxRM/QTL/intScans.bin')
load('~/1000BYxRM/QTL/intScans.bin')
#save.image('~/Desktop/040112.bin')

int_FDR10=sapply(intScans, function(x) {
    pFDR = attr(x, 'pFDR') 
    plot(pFDR, as.numeric(names(pFDR))  )
    if(min(pFDR, na.rm=T)>.25) { return(10) } else { ss = smooth.spline(pFDR, as.numeric(names(pFDR)))
        predict(ss, .10)$y   }
})
names(int_FDR10)=names(fQTLs_FDR05)
 

intScanPeaks= foreach (t = 1:length(intScans) ) %dopar% {
    ipeaks.FDR10      =  define.peak2(intScans[[t]], lodcolumn='all', th=int_FDR10[t],  si=1.5, window.size =30, round=3)
    ipeaks.FDR10      =  cleanPeaks(ipeaks.FDR10)
    ipeaks.FDR10      =  map.peak(ipeaks.FDR10)
    ipeaks.FDR10      = tryCatch({    
                intPeakMarker = find.marker(cross, ipeaks.FDR10$chr, ipeaks.FDR10$cM) 
                mainPeakLOC  = t(sapply(ipeaks.FDR10$trait, find.markerpos, cross=cross))
                ipeaks.FDR10     =  data.frame(ipeaks.FDR10[,1], mainPeakLOC, intPeakMarker, ipeaks.FDR10[,c(2,3)])
                names(ipeaks.FDR10) =c('main.Marker', 'main.chr', 'main.pos', 'int.Marker', 'int.chr', 'int.pos')
              
                ipeaks.FDR10$main.Marker = as.character(ipeaks.FDR10$main.Marker )
                ipeaks.FDR10$int.Marker  = as.character(ipeaks.FDR10$int.Marker )
                ipeaks.FDR10$main.chr = as.numeric(as.character(ipeaks.FDR10$main.chr ))
                ipeaks.FDR10$main.pos = as.numeric(ipeaks.FDR10$main.pos )
                ipeaks.FDR10$int.chr = as.numeric(as.character(ipeaks.FDR10$int.chr ))
                ipeaks.FDR10$int.pos = as.numeric( ipeaks.FDR10$int.pos )

                ipeaks.FDR10.m = (unique(ipeaks.FDR10$main.Marker))
                ipeaks.FDR10.m = ipeaks.FDR10[match(ipeaks.FDR10.m, ipeaks.FDR10$main.Marker), c(1,2,3)]
        
                temp = ipeaks.FDR10
                for(j in 1:nrow(ipeaks.FDR10) ) {
                        cc = as.numeric(ipeaks.FDR10$int.chr[j])
                        pp = ipeaks.FDR10$int.pos[j]
                        cmat = cc-ipeaks.FDR10.m$main.chr
                        pmat = abs(pp-ipeaks.FDR10.m$main.pos)
                        matchme = which(cmat==0 & pmat<20)
                        if(length(matchme)==1) {
                           temp[j, 4:6] = ipeaks.FDR10.m[matchme,] 
                        }
                }
             return(temp)
            },error=function(e) {   return(NULL)        }) 
    return(ipeaks.FDR10)
}
names(intScanPeaks)=names(fQTLs_FDR05)
save(intScanPeaks, file='~/1000BYxRM/QTL/intScanPeaks.bin')

intcount = sapply(intScanPeaks, function(x) {nrow(x) })
intcount = sapply(intcount, function(x) { if(is.null(x)){0} else{x} })

vInt = foreach( t = 1:length(intScans)) %dopar% {
    print(t)
    additive.markers = rownames(fQTLs_FDR05[[t]]$qtlMarkersCP)
    additive.pos = find.markerpos(cross, additive.markers)
    rqtl = makeqtl(cross, chr=additive.pos$chr, pos = additive.pos$pos, what='prob')
    rfit = fitqtl(cross, pheno.col=t, qtl=rqtl)
    rvexp = rfit$result.full[1,'%var']/100
    if(!is.null(intScanPeaks[[t]])) {

        interactions     = intScanPeaks[[t]]
        uniqueQTLs       = unique((c(additive.markers, interactions$main.Marker, interactions$int.Marker)))
        names(uniqueQTLs) = paste('Q', 1:length(uniqueQTLs), sep='')

        interactions      = data.frame(interactions, names(uniqueQTLs[match(interactions$main.Marker, uniqueQTLs)]), 
                                                 names(uniqueQTLs[match(interactions$int.Marker, uniqueQTLs)]))
        names(interactions)[c(7,8)]=c('int1', 'int2')

        iqs = paste(as.character(interactions$int1), as.character(interactions$int2), sep='*')
        iqs = paste(iqs, collapse = ' + ')

        qtlpos =find.markerpos(cross, uniqueQTLs) 
        mqtl = makeqtl(cross, chr=qtlpos$chr, pos=qtlpos$pos, what='prob')
        afit = fitqtl(cross, pheno.col=t, qtl=mqtl)
        imodel = paste( attr(afit, 'formula'), iqs ,sep=' + ')
        ifit = fitqtl(cross, pheno.col=t, qtl=mqtl, formula= imodel)
        vexp = ifit$result.full[1,'%var']/100
    } else{   vexp =NA; ifit=NA  }
   return(list(rvexp = rvexp, vexp=vexp, ifit=ifit) ) }
names(vInt)=names(intScans)


via=sapply(vInt, function(x) {tryCatch( {x$ifit$result.drop }, error =function(e) {NULL})   } )

via.df = do.call('rbind', via)
via.l = sapply(via, nrow)
uvl=unlist(via.l)
traitN=rep(names(uvl), uvl)
qtlN = rownames(via.df)
via.df=data.frame(traitN, qtlN, via.df)
i.index=grep('@.*@', via.df$qtlN)
via.df = via.df[i.index,]
iH2[iH2$Trait == via.df[,1], 'BroadH2'] 
via.df$BH2= iH2$BroadH2[match(via.df[,1], iH2$Trait)]
via.df$X.var/(via.df$BH2*100)

pdf(file= '~/Desktop/NaVarExp.pdf' ,width =8 ,height =8)   

naexp =sapply(vInt , function(x) {x$vexp-x$rvexp})
aexp =sapply(vInt , function(x) {x$rvexp})
fexp =sapply(vInt , function(x) {x$vexp})

naexp.z = naexp
naexp.z[is.na(naexp.z)]=0

# end Detecting QTL-QTL  interactions###############################################################################################
