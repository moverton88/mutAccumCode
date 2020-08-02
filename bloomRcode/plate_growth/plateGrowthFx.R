# Functions
# makeGrid()
# segmentImages()
# binaryImageFilter()
# getColonies()

runPlateGrowth = function(n, SKey, AllStrainLayouts, AllColonyGrids, params) {
    attach(params, warn.conflicts=FALSE)
    scann = names(SKey)[n]
    print(scann) 
    thisScan = SKey[[scann]]
    ufile   = unique(thisScan$FullFileName)
    sfile   = unique(thisScan$FILENAME)
    
    allsuper  = (with(thisScan, {paste(STRAINLAYOUT,PERMUTATIONGROUP, PLATEPOSITION, CONDITION, CONCENTRATION, CONTROL) }) )
    usuper    = unique(allsuper)
    uthisScan=  thisScan[match(usuper, allsuper),]

    scanImg = readScan(ufile, bit.depth=bit.depth, hp=hp)
    ColonyGrids   = AllColonyGrids[[sfile]]
    
    StrainLayouts = AllStrainLayouts[uthisScan$STRAINLAYOUT] 
    plateIndices =  1:length(StrainLayouts)
    if(rep.4==TRUE) {StrainLayouts = lapply(StrainLayouts, rep4) }
    
    extracted.plates = extractPlateImages(scanImg, ColonyGrids, borderoffset=border.offset, 
                                          diagnose=diagnose.grids, faceup=faceup)
    print('Done Extracting Plates')
    plates.rotated = extracted.plates$platesRotated
    spot.seeds     = extracted.plates$spot.seeds
   
    rm(extracted.plates)
    plate.masks = mapply(segmentPlate, plates.rotated, spot.seeds, SIMPLIFY=FALSE)
    
    print('Done Segmenting Plates')
    features = mapply(moments, plate.masks, plates.rotated, SIMPLIFY=FALSE) 
    hfeatures = lapply(plate.masks, hullFeatures)
    afeatures = mapply(cbind, features,hfeatures, SIMPLIFY=FALSE)

    nfeatures = mapply(nameSpots,afeatures[plateIndices], StrainLayouts[plateIndices],   SIMPLIFY=FALSE)

    print('Done Extracting Features')
    names(nfeatures)=as.character(apply(uthisScan[,1:(ncol(uthisScan)-1)], 1, paste, collapse='::'))
   
    #print(plateIndices)
    #saveImages(plates.rotated[plateIndices], plate.masks[plateIndices], nfeatures, scann, images.dir)     
    gc()
    detach('params')
    return(nfeatures)
}

segmentPlate <- function (plate.raw, spot.seeds) {
    print('Segmenting Plate')
      
    tplate=plate.raw
    #thresh(plate.raw, 75,75,.01)
    if( quantile(plate.raw,.25)<50) {
        km = kmeans(log10(as.vector(plate.raw)+.001),c(1.5,2))
    } else{   km = kmeans(log10(as.vector(plate.raw)+.001),2)   }
    kmin = which.min(km$centers)
    kmax = which.max(km$centers)
    
    if(((km$size[kmax]/km$size[kmin])<.01)  | (max(km$centers)<1.8)  ) {km$cluster[km$cluster==kmax]=kmin }

    if(kmin==2) {km$cluster = abs(km$cluster-3)}
    
    tplate@.Data=matrix(km$cluster-1,nrow(tplate),ncol(tplate))

    seeds = matrix(0, nrow(plate.raw), ncol(plate.raw))
    seeds[is.numeric(seeds)]=0
    
    ycirc    <- makeBrush(9, 'disc', step=TRUE)

    for(i in 1:nrow(spot.seeds) ){
        x = spot.seeds[i,'X']
        y = spot.seeds[i, 'Y']
        
        # boundary cases 
        if( (x-4)<1) {x=5}
        if( (y-4)<1) {y=5}
        if(x>nrow(seeds)-4 ) {x = nrow(seeds)-5 }
        if(y>ncol(seeds)-4 ) {y = ncol(seeds)-5 }

        seeds[(x-4):(x+4), (y-4):(y+4) ] =ycirc*i
    }
    plate.mask = propagate( plate.raw, seeds, mask=tplate)
    return(plate.mask)
}

# Return size of objects in workspace in MB
plate.types  <- list('96'=c(8,12), '384'=c(16,24), '1536'=c(32,48) )

getObjS      <- function(y){  sort(sapply(y, function(x) { object.size(get(x))/1e6} ), decreasing=T)}

# Force garbage collection by calling it repeatedly , ghetto
gcp      <- function() { invisible(replicate(20, gc())) }

# store arrays as bit vectors to reduce memory usage 32 fold
# store original dimensions in 'dimBit' attribute
# if I cared I could force input to be a class I define
array2bitvec <- function (x) { 
    y <- as.bit(x)
    attr(y, 'dimBit') <- dim(x)
    return(y)
}

makeKey = function(params) {
    attach(params, warn.conflicts=FALSE)
    Key              <- read.delim(paste(keys.dir, key.file, sep=''), stringsAsFactors=FALSE, sep='\t', header=TRUE)
    names(Key)       <- toupper(names(Key))
    Key$FullFileName <- paste(images.dir, Key$FILENAME, sep='')
    detach('params')
    return(Key)
}

makeStrainLayouts = function(Key, params) {
    attach(params, warn.conflicts=FALSE)
    # Read in strain layout files 
    AllStrainLayouts        <- as.list(unique(Key$STRAINLAYOUT))
    names(AllStrainLayouts) <- unlist(AllStrainLayouts)
    AllStrainLayouts        <- lapply(AllStrainLayouts, function(x) { 
                            file.in <- paste(keys.dir, x, sep='')
                            lout = read.delim(file.in, sep='\t', header=FALSE, stringsAsFactors=FALSE) 
                            lout[lout=='']='BLANK'
                            return(lout)
                        })
    detach('params')                
    return(AllStrainLayouts)
}

makeColonyGrids<- function(params) {

    attach(params, warn.conflicts=FALSE)
    #read in colony grids
    AllColonyGrids <- read.delim(paste(keys.dir, coordinates.file, sep=''), stringsAsFactors=FALSE, sep='\t', header=TRUE)
    names(AllColonyGrids) <- toupper(names(AllColonyGrids))

    # inches to pixels conversion, if necessary
    if(inches==TRUE & hp==FALSE ){ AllColonyGrids$X=AllColonyGrids$X*(3200/8);     
                                   AllColonyGrids$Y=AllColonyGrids$Y*(4000/10) }
    if(inches==TRUE & hp==TRUE ) { AllColonyGrids$X=AllColonyGrids$X*(3149/7.87);  
                                   AllColonyGrids$Y=AllColonyGrids$Y*(4093/10.23) }
    
    # split by plate position and filename
    AllColonyGrids <- split(AllColonyGrids[,c('X','Y','PLATEPOSITION')], AllColonyGrids$FILENAME)
    AllColonyGrids <- lapply(AllColonyGrids, function(x) {split(x[,c('X','Y')], x$PLATEPOSITION) } )

    # Seed starting grid for spot detection
    AllColonyGrids <- makeGrids(AllColonyGrids, plate.type)

    detach('params')
    return(AllColonyGrids)
}


# convert bit vectors to arrays (undo array2bitvec)
bitvec2array <- function(x) { array(as.integer(x), attr(x,'dimBit') )  }

readScan <- function(filein, bit.depth=8, hp=FALSE) {   
    test =     imageData(readImage(filein)*(2^bit.depth))
    if(hp==TRUE) { test = flip(rotate(test,180)) }
    storage.mode(test)='raw'
    return(test)
}
makeGrids <- function(coords, plate.type) {
   grids = list()
   for (lay in names(coords) ) {
        # x and y confusion (X is really number of columns in Y space ... fix this) 
        avg.bounds   =   lapply(coords[[lay]], function(x){
                             x.c = sort(x$X)
                             y.c = sort(x$Y)
                             list( xmin=round(mean(x.c[c(1,2)])), xmax=round(mean(x.c[c(3,4)])),
                                   ymin=round(mean(y.c[c(1,2)])), ymax=round(mean(y.c[c(3,4)])))})
        for(i in 1:length(avg.bounds) ) {
              av=avg.bounds[[i]]
              xsz= av$xmax-av$xmin
              ysz= av$ymax-av$ymin
              if(xsz>ysz) {
               gridme = expand.grid( 
                      round(seq(av$xmin, av$xmax, length.out=max(plate.type)) ), 
                      round(seq(av$ymin, av$ymax, length.out=min(plate.type))))

              colnames(gridme)=c('X','Y') 
              } else {
                gridme = expand.grid( 
                      round(seq(av$xmin, av$xmax, length.out=min(plate.type))),
                        round(seq(av$ymin, av$ymax, length.out=max(plate.type)))) 

              colnames(gridme)=c('X','Y') 
              gridme = gridme[order(gridme[,1], gridme[,2]),]
              rownames(gridme)=seq_along(rownames(gridme))
            }
            #    plot(gridme, asp=T)
            # text(gridme, rownames(gridme), asp=TRUE)
             grids[[lay]][[names(avg.bounds)[i]]]=gridme
        }
  }
  return(grids)
}

#given a grid ... quadruplicate it like the rotor
rep4 = function(strnl) {
     strnl=  apply(strnl,c(1,2), as.character)
     strnl = gsub(" ", "", strnl)
     fillme = matrix(" ", nrow(strnl)*2, ncol(strnl)*2)
     xtoit = seq(1,nrow(fillme),2)
     ytoit = seq(1,ncol(fillme),2)
     for(x in 1:nrow(strnl)) {
           for(y in 1:ncol(strnl)) {
              xx=xtoit[x]
              yy=ytoit[y]
              repme = as.character(strnl[x,y])
              fillme[xx,yy]=repme
              fillme[xx+1,yy]=repme
              fillme[xx,yy+1]=repme
              fillme[xx+1,yy+1]=repme  }}

    return(fillme)
}

rotateME = function(i,img) {
       temp=rotate(img,i) 
       xmarg=apply(temp,1,sum)
       xmarg = xmarg[which(peaks(xmarg, 10)) ]
       xmarg = xmarg[5:(length(xmarg)-5)]
         
       ymarg=apply(temp,1,sum)
       ymarg = ymarg[which(peaks(ymarg, 10)) ]
       ymarg = ymarg[5:(length(xmarg)-5) ]
       return( sum(xmarg)+sum(ymarg))
}

alignGrid.X = function(i, img, ogrid.x) {
    xsum =    sum( apply(img[ogrid.x+i, ],1,sum))
    return(xsum)
}
 alignGrid.Y = function(i, img, ogrid.y) {
    ysum = sum( apply(img[,ogrid.y+i],2,sum))
    return(ysum)
}   


# return list of split reoriented images and corresponding reoriented grids
extractPlateImages <- function(image.raw, ColonyGrids, borderoffset=26, diagnose=FALSE, faceup=FALSE) {
    storage.mode(image.raw)='integer'
    fullscan = image.raw

    ranges = sapply(ColonyGrids, function(x){c(range(x$X), range(x$Y))} )
    ranges[c(1,3),] = ranges[c(1,3),]-borderoffset
    ranges[c(2,4),] = ranges[c(2,4),]+borderoffset
    rownames(ranges) = c('xmin', 'xmax', 'ymin', 'ymax')

    newColonyGrids = list()
    plates=list() 

    for (plate in names(ColonyGrids) )  {
        pranges = ranges[,plate]
        gridme=data.frame(
                ColonyGrids[[plate]]$X-pranges['xmin'],
                ColonyGrids[[plate]]$Y-pranges['ymin'])
        # if skinny make fat 
        if( (pranges['ymax']-pranges['ymin']) > (pranges['xmax']-pranges['xmin']) ) { 
            gridme = rev(gridme) 
            colnames(gridme)=c('X', 'Y')
            newColonyGrids[[plate]]=gridme
            plates[[plate]] = rotate( fullscan[ranges['xmin',plate]:ranges['xmax',plate], 
                                               ranges['ymin',plate]:ranges['ymax',plate]], 270)
         }
         else {
            colnames(gridme)=c('X', 'Y')
            newColonyGrids[[plate]]=gridme
            plates[[plate]] = fullscan[ranges['xmin',plate]:ranges['xmax',plate], 
                                               ranges['ymin',plate]:ranges['ymax',plate]]
         }
     }
     
     platesRotated=plates
     # gridAdjust=newColonyGrids

     spots = list()
     for (plate in names(ColonyGrids)){
         # print(plate)
         gridme = data.frame(
                 newColonyGrids[[plate]]$X,
                 #+as.numeric(gridAdjust['X', plate]),
                 newColonyGrids[[plate]]$Y)
                 #+as.numeric(gridAdjust['Y', plate]))
         colnames(gridme) = c('X', 'Y') 
         spots[[plate]]=gridme
        if(diagnose ==TRUE) {
            temp = platesRotated[[plate]]/256
            temp[spots[[plate]]$X, spots[[plate]]$Y]=1
            display(temp , title = paste(plate, 'rotated with fixed grid') ) 
        }

     }
    
     if(faceup==TRUE) { 
         for(plate in names(ColonyGrids) ) {
         #flop plate so that A1 is in upper right
             platesRotated[[plate]]=rotate(t(platesRotated[[plate]]), 90)
             
             tmpspot = spots[[plate]]
             tmpimg  = platesRotated[[plate]]
             tmpimg[is.numeric(tmpimg)]=0
             tmpimg[tmpspot$X, tmpspot$Y]=1
             tmpimg=rotate(t(tmpimg), 90)
             
             nspot=data.frame(which(tmpimg>0, arr.ind=TRUE))
             colnames(nspot)=c('X','Y')
             spots[[plate]]=nspot
         }
     }
     return(list('platesRotated'=platesRotated, 'spot.seeds'=spots) )
}

    
#segmentImages()--------------------------------------------------------
# input  = 1) array object of class 'Image' and 
#          2) 2d weighted smoothing diameter in pixels
#          3) vector of length 3 with width, height, and deviation from mean
#             for rolling local image thresholding
# output = array object of class 'Image' that contains only binary data



nameSpots = function(feats, layouts) {  rownames(feats)=as.vector(t(layouts)); return(feats) }

saveImages = function(plates.rotated, plate.masks, nfeatures, scann, images.dir) {
    dout = paste(images.dir, 'out/', sep='')
    dir.create(dout, showWarnings=FALSE)

    for(plate in 1:length(plates.rotated)) {
        print(plate)
        praw = plates.rotated[[plate]]
        mask = plate.masks[[plate]] 

        temp1 = paintObjects(mask, rgbImage(red=praw/256,blue=praw/256, green= praw/256 ))
        print('Done pO')
        temp2 = drawtext(temp1, cbind(nfeatures[[plate]][,3]-10, nfeatures[[plate]][,4]), 
                as.character(rownames(nfeatures[[plate]])))
        #, font=drawfont('sans', antialias=FALSE,size=7))
        print('Done dT')
        scann=gsub('tif','jpeg',scann)
        writeImage(temp2, paste(dout,'LAB_', names(plates.rotated)[plate], '_', scann, sep=''))
    }
 }


avgPlateControls = function(x) {
     tempvecs = rbind(sapply(x, as.vector))
     nmat = matrix(rowMeans(tempvecs, na.rm=T), dim(x[[1]]))
     attr(nmat, 'dimnames')=attr(x[[1]], 'dimnames')
     attr(nmat, "batch")= attr(x[[1]], 'batch')
     attr(nmat, "strainlayout")= attr(x[[1]], 'strainlayout')
     attr(nmat, "permutationgroup")= attr(x[[1]], 'permutationgroup')
     attr(nmat, "condition")= attr(x[[1]], 'condition')
     attr(nmat, "concentration")= attr(x[[1]], 'concentration')
     attr(nmat, "control")= attr(x[[1]], 'control')
     return(nmat)
}


processResults = function(Results) {
    nResults= foreach(n=1:length(Results)) %dopar%    { withinPlateNorm(Results[[n]]) }
    names(nResults)=names(Results)
    return(nResults)
}

#######################################################################################################
localPfit = function(vecin, x,y) {
   m = locfit(vecin~lp(x,y, scale=FALSE),lfproc=locfit.robust)
   predx = predict(m, data.frame(vecin, x,y))
   return(vecin-predx) }
#######################################################################################################

parsePlateName = function(n) {   unlist(strsplit(n, '::'))  }


#######################################################################################################
withinPlateNorm = function(pall) {

   pa=pall[,-c(19,20)]
   avg.int  = pa[,'m.int']/pa[,'m.pxs']
   pa = cbind(pa, avg.int)
   print(paste(attr(pall, 'scanfile'), attr(pall, 'plateposition'))) 
   
   normalized =  apply(pa[,c('g.effr', 'm.pxs')], 2, function(x){ 
                    tryCatch( { return(localPfit(x, pa[,'m.x'], pa[,'m.y'])+mean(x,na.rm=T))}, 
                       error = function(e){return(x)} ) })

   colnames(normalized)=paste(colnames(normalized),'norm', sep='.')
    
   r = cbind(pa, normalized)
   attr(r,'batch')                 =     attr(pall,'batch')
   attr(r,'strainlayout')           =    attr(pall,'strainlayout')
   attr(r,'permutationgroup')       =    attr(pall,'permutationgroup')
   attr(r,'scanfile')               =    attr(pall,'scanfile')
   attr(r,'plateposition')          =    attr(pall,'plateposition')
   attr(r,'condition')              =    attr(pall,'condition')
   attr(r,'concentration')          =    attr(pall,'concentration')
   attr(r,'control')                =    attr(pall,'control')
   return(r)
}
#######################################################################################################

extractAttributes = function(r) {
    batch            =attr(r,'batch')                
    strainlayout     =attr(r,'strainlayout')         
    permutationgroup =attr(r,'permutationgroup')     
    scanfile         =attr(r,'scanfile')             
    plateposition    =attr(r,'plateposition')        
    condition        =attr(r,'condition')            
    concentration    =attr(r,'concentration')        
    control          =attr(r,'control')
    return(c(batch,           
             strainlayout,    
             permutationgroup,
             scanfile,        
             plateposition,   
             condition,       
             concentration,   
             control))         
}
getAnnotation= function(r) {
    d=data.frame(do.call('rbind', lapply(r, extractAttributes)), stringsAsFactors=F)
    rownames(d)=NULL
    colnames(d)=c( 'batch', 'strainlayout', 'permutationgroup', 'scanfile', 'plateposition', 'condition', 'concentration', 'control')
    return(d)
}


attachAttributes = function(r,n) {
            attr(r,'batch')=n[1]
            attr(r,'strainlayout')=n[2]
            attr(r,'permutationgroup')=n[3]
            attr(r,'scanfile')=n[4]
            attr(r,'plateposition')=n[5]
            attr(r,'condition')=n[6]
            attr(r,'concentration')=n[7]
            attr(r,'control')=n[8]
            return(r)
}

#######################################################################################################
unlistify = function(rr, batch=1) {
    # Assuming input is an R binary file with variable called 'Results' from 
    # plate growth script
    plate.names =as.vector(unlist(sapply(rr, function(x){names(x)})))
    plate.names = gsub(' ', '', plate.names)
    plate.names = paste(batch, plate.names, sep='::') 
    R = list()
    cnt=1
    for(x in 1:length(rr)){
        for(y in 1:length(rr[[x]])){
            R[[plate.names[cnt]]]=rr[[x]][[y]]
            n=parsePlateName(plate.names[cnt])
            R[[plate.names[cnt]]]=attachAttributes(R[[plate.names[cnt]]], n)
            cnt=cnt+1
        }
    }
    R = R[!is.na(names(R))]
    return(R)   }
#######################################################################################################

unlistify2 = function(rr) {
    # Assuming input is an R binary file with variable called 'Results' from 
    # plate growth script
    R = list()
    plate.names=as.vector(unlist(sapply(rr, function(x){names(x)})))
    cnt=1
    for(x in 1:length(rr)){
        for(y in 1:length(rr[[x]])){
            R[[plate.names[cnt]]]=rr[[x]][[y]]
            n=parsePlateName(plate.names[cnt])
            R[[plate.names[cnt]]]=attachAttributes(R[[plate.names[cnt]]], n)
            cnt=cnt+1
        }
    }
    R = R[!is.na(names(R))]
    return(R)   }

reducePhenoData = function(phenoData, strains_used) {
    return(lapply(phenoData, function(x){
        mmeans=   names(x$segs.mean) %in% strains_used
        mall  =   names(x$segs.all) %in% strains_used
        x$segs.mean = x$segs.mean[mmeans]
        x$segs.all = x$segs.all[mall]
        return(x)
    }))
}


groupStrainPhenos = function(plist, pheno='jrad.locf'){
    # tech.rep.strains='YLK1993|YLK1950|YLK1879') {
    # plist.annot = data.frame(do.call('rbind',  strsplit(names(plist), '_') ) )
    # names(plist.annot) = c('SET', 'G384', 'FILE', 'SCANNERP', 'DRUG', 'COND')
    # names(plist)=paste( 1:nrow(plist.annot),  plist.annot[,2], sep='_')

    pall = do.call('rbind', plist)
    strains = factor(rownames(pall))
    #    p.384 = factor(rep(plist.annot$G384, each=384))
    #    p.96 = factor(as.vector(unlist(l384_96[plist.annot$G384])))

    strains.all  =  split(as.vector(pall[,pheno]), strains)

    segs.all = strains.all[grep('^A', names(strains.all))]
    segs.mean  =  sapply(segs.all, mean, na.rm=T)
    
    notsegs.all = strains.all[grep('^A', names(strains.all), invert=TRUE)]
    notsegs.mean  = sapply(notsegs.all, mean, na.rm=T)

    #wparents.sz = s[grep(tech.rep.strains, names(s))]
    #wparents.sz$segs = as.numeric(seg.median)
    #sdata = median per segregant
    #wparent.sz = parent data
    #s = all data (segregants and WI)
    #s.segs = segregants only
    return(list(segs.mean=segs.mean, notsegs.mean=notsegs.mean, segs.all=segs.all, notsegs.all = notsegs.all) )}


getTable = function(phenoData, value='segs.mean') {
    
    # without normalization for growth in YPD 10/04/11
    # segm = sapply(phenoData, function(x) { x$segs.median })
    
    # want conditions as colnames and segs as rownames
    segm = sapply(phenoData, function(x) {x[value]})
    
    segnames= sort(unique(as.vector(unlist(sapply(segm, names)))))
    
    segd= sapply(segm, function(x) {
        missingStrains = segnames[!segnames %in% names(x)]
        y=x
        if(length(missingStrains)>0) {
            missingD = rep(NA, length(missingStrains)) 
            names(missingD)=missingStrains
            y=c(x,missingD)
        }
        y=y[order(names(y))]
        return(y)  
    })  
    colnames(segd)=gsub(paste('.',value, sep=''),'',colnames(segd))
    return(segd)
}

getTable2 = function(pD, value='jrad.locf') {
    
    # want conditions as colnames and segs as rownames
    segm = lapply(pD, function(x) {x[,value]})
    segm2= lapply(segm, function(y) {
                #rmme=  grep('Par_YLK1950|Par_YLK1879|Par_YLK1993|Par_Blank',names(y))
             rmme=grep('^A',names(y))
             return(y[rmme])
         })

    segnames= sort(unique(as.vector(unlist(sapply(segm2, names)))))
    
    segd= sapply(segm2, function(x) {
        missingStrains = segnames[!segnames %in% names(x)]
        y=x
        if(length(missingStrains)>0) {
            missingD = rep(NA, length(missingStrains)) 
            names(missingD)=missingStrains
            y=c(x,missingD)
        }
        y=y[order(names(y))]
        return(y)  
    })  
    return(segd)
}

# Combine information from YPD normalization with existing feature information
makedtf = function(x,y) {
       lapply(x, function(z) { 
            print( paste(attr(z, 'condition'),attr(z, 'concentration'), attr(z, 'strainlayout')))
            dtf=z
            t1=lm(scale(z[,'g.effr.norm'])~scale(y[,'g.effr.norm']), na.action='na.exclude')
            #t1=lm(scale(z[,'g.effr.norm'])~scale(y[,'g.effr.norm']), na.action='na.exclude')
            do1=residuals(t1)
            #z[,'g.effr.norm']-(as.vector(predict(t1,y=y[,'g.effr.norm'])))
            dtf= cbind(dtf, do1)
            ncc = dim(dtf)[2]
            colnames(dtf)[ncc] ='g.effr.norm.ypdc'
            return(dtf)})
}
makedtf_include_control_data = function(x) { 
            r=cbind(x,x[,'g.effr.norm'])
            ncc = dim(r)[2]
            colnames(r)[ncc]='g.effr.norm.ypdc'
            attr(r,'batch')                 =     attr(x,'batch')
            attr(r,'strainlayout')           =    attr(x,'strainlayout')
            attr(r,'permutationgroup')       =    attr(x,'permutationgroup')
            attr(r,'scanfile')               =    attr(x,'scanfile')
            attr(r,'plateposition')          =    attr(x,'plateposition')
            attr(r,'condition')              =    attr(x,'condition')
            attr(r,'concentration')          =    attr(x,'concentration')
            attr(r,'control')                =    attr(x,'control')
            return(r)
}
      



writeTXTfiles = function(outDirectoryTXTfiles, Results.norm) {
    dir.create(outDirectoryTXTfiles)
    for(n in names(Results.norm)){
        n2 =n
        n2=gsub('/', '__', n2)
        n2=gsub(':', '_', n2)
        write.table(Results.norm[[n]], file=paste(outDirectoryTXTfiles,n2,'.txt',sep=''), 
                sep='\t', row.names=T,col.names=NA,quote=F) }
    }

# Diagnostics 
#all.size     = sapply(Results3, function(pall) {pall[,'m.pxs']})
#all.rad.diff = sapply(Results3, function(pall) { abs(log2((2*sqrt(pall[,'m.l1']))/(2*sqrt(pall[,'m.l2'])))  ) })
#all.avg.int  = sapply(Results3, function(pall) {pall[,'m.int']/pall[,'m.pxs']})
#all.jrad     = sapply(Results3, function(pall) {sqrt(pall[,'m.pxs']/pi)})
#all.avg.rad  = sapply(Results3, function(pall) {(2*sqrt(pall[,'m.l1'])+2*sqrt(pall[,'m.l2']))/2 })
#all.diff.exp = abs(all.avg.rad-all.jrad)
#all.ecc      = sapply(Results3, function(pall) {pall[,'m.ecc']})
#all.theta    = sapply(Results3, function(pall) {pall[,'m.theta']})
#all.int      = sapply(Results3, function(pall){pall[,'m.int']})
#max.rad      = sapply(Results3, function(pall){pall[,'m.l1']})
#
## investigate shape factor or acircularity measures as filters for errors 
#all.p     = sapply(Results3, function(pall) {pall[,'g.p']})
#all.expp  = all.
#
#all.sf     = sapply(Results3, function(pall) {pall[,'g.sf']})
#all.edge     = sapply(Results3, function(pall) {pall[,'g.edge']})
#all.acirc     = sapply(Results3, function(pall) {pall[,'g.acirc']})
#all.pdm     = sapply(Results3, function(pall) {pall[,'g.pdm']})
#all.pdsd     = sapply(Results3, function(pall) {pall[,'g.pdsd']})
#all.effr     = sapply(Results3, function(pall) {pall[,'g.effr']})
#
#all.expp=2*pi*all.avg.rad
#
#bs = all.edge>1 | (all.sf>1 & all.effr>20)

#Extracted object features are:
#        • ‘g.x ,g.y’ - coordinates of the geometric center.
#        • ‘g.s’ - size in pixels.
#        • ‘g.p’ - perimeter in pixels.
#        • ‘g.pdm’ - mean distance from the center to perimeter.
#        • ‘g.pdsd’ - standard deviation of the distance to perimeter.
#        • ‘g.effr’ - effective radius (is the radius of a circle with
#          the same area).
#        • ‘g.acirc’ - acircularity (fraction of pixels outside of the
#          circle with radius ‘g.effr’).
#        • ‘g.sf’ - shape factor, equals to (‘g.p/ ( 2*sqrt(i*g.s))’).
#        • ‘g.edge’ - number of pixels at the edge of the image.
#        • ‘g.theta’ - hull orientation angle, in radians. See above.
#        • ‘g.l1’ - largest eigeinvalue of the covariance matrix. See
#          above.
#        • ‘g.l2’ - lowest eigenvalue of the covariance matrix. See
#          above.
#        • ‘g.ecc’ - eccentricity, equals to sqrt(1-g.l2/g.l1). See
#          above.
#        • ‘g.I1, g.I2’ - first and second Hu's
#          translation/scale/rotation invariant moment. See above.

