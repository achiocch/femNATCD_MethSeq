require("kableExtra")
require("tidyverse")
require("compareGroups")
require("RCircos")
require("scales")

# preprocessing fucntions ####
yesnofac= function(x){factor(x, c(0,1), c("no", "yes"))}


# summary table function####

table_sumstat = function(DF, columns, groupfactor){
  tmpdf= DF %>% select(all_of(c(columns, groupfactor)))
  table <- compareGroups(formula(paste0(groupfactor, "~.")), data = tmpdf)
  pvals <- getResults(table, "p.overall")
  export_table <- createTable(table)
  return(export_table)
}


# scaling functions####
iqr = function (x){
  a = quantile(x,c(.25, .75))
  return(a[2]-a[1])}

robust_scaling = function(x){
  return((x-median(x))/iqr(x))}

minmax_scaling = function(x){
  a = range(x, na.rm = T)
  return((x-min(x,na.rm=T))/(a[2]-a[1]))}

# genomic functions####

getuniquegenes = function(vectorgenes){
  restmp = sapply(vectorgenes, strsplit, split=",")
  rescomb = unlist(restmp, recursive = TRUE)
  rescomb = unique(rescomb)
  return(rescomb)
}


chrtonum = function(vector) {
  x = gsub("M", 25,
           gsub("Y", 24,
                gsub("X", 23,
                     gsub("chr", "", vector))))
  x = as.numeric(x)
  return(x)}


chrtochar= function(vector) {
  x = paste0("chr", gsub(25, "M",
                         gsub(24, "Y",
                              gsub(23,"X",as.character(vector)))))
  return(x)}


data(UCSC.HG19.Human.CytoBandIdeogram)

add_legend=function(row_offset,ParamsRC, min.value, max.value, leglab){
  width = ParamsRC$plot.radius/5
  height = ParamsRC$plot.radius/25

  xright <- ParamsRC$plot.radius
  ytop<- ParamsRC$plot.radius - (row_offset*height+0.2*height)
  xleft <-   xright-width
  ybottom <- ytop-height

  colorscodes <- rev(RCircos.Get.Heatmap.Color.Scale(ParamsRC$heatmap.color))
  step=seq(xright, xleft, length.out = length(colorscodes)+1)
  for(s in 1:length(colorscodes)){
    rect(step[s],ybottom, step[s+1], ytop, col=colorscodes[s], border=NA)
  }
  text(xleft,(ybottom+ytop)/2, labels = leglab, pos=2, cex=0.5, offset=1)
  text(xleft,(ybottom+ytop)/2, labels = round(min.value,1), pos=2, cex=0.4, offset=0.2)
  text(xright,(ybottom+ytop)/2, labels = round(max.value,1), pos=4, cex=0.4, offset = 0.1)
}

plotcircos= function(plotdata, title = "CircosPlot",
                     targets = list(inner=c(), outer=c()),
                     labelsidx = c(),
                     labcol="gene",
                     pvalident = "log10_P",
                     pvallog=FALSE,
                     cutoffpval = 8,
                     filename=paste0(Home,"/output/circos_test.png")){

  targets_inner = targets[["inner"]]
  targets_outer = targets[["outer"]]

  labdat=plotdata[rownames(plotdata) %in% labelsidx, c(colnames(plotdata)[c(1:3)],labcol)]

  chr.exclude <- c("chrY")
  cyto.info <- UCSC.HG19.Human.CytoBandIdeogram
  tracks.inside <- length(targets_inner)
  tracks.outside <- length(targets_outer)
  if(length(labels)>0){tracks.outside = tracks.outside+2}

  RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)

  svg(filename, width = 7, height = 7)

  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()

  pvalindex=grepl(tolower(pvalident), tolower(colnames(plotdata)))
  pvalcols = colnames(plotdata)[pvalindex]

  if(!pvallog){
    tmp = plotdata %>% mutate_if(pvalindex, function(x){return(-log10(x))})
  } else {
    tmp = plotdata
  }

  min.value.pval = 0
  max.value.pval = cutoffpval

  if (length(targets_outer)>0){
    n=length(targets_outer)
    legendreverser = c(n:1)

    for (i in 1:length(targets_outer)){
      target1=targets_outer[i]
      tmp2=tmp[!is.infinite(plotdata[,target1]), ]
      if(target1 %in% pvalcols){
        min.value =  min.value.pval
        max.value= max.value.pval
        RC.param <- RCircos.Get.Plot.Parameters()
        RC.param['heatmap.color'] <- "BlackOnly"
        RCircos.Reset.Plot.Parameters(RC.param)

      } else {
        min.value = min(tmp2[,target1], na.rm=T)
        max.value = max(tmp2[,target1], na.rm=T)
        RC.param <- RCircos.Get.Plot.Parameters()
        RC.param['heatmap.color'] <- "YellowToRed"
        RCircos.Reset.Plot.Parameters(RC.param)
      }

      RCircos.Heatmap.Plot(tmp2,
                           data.col = which(colnames(tmp2)==target1),
                           min.value = min.value,
                           max.value = max.value,
                           genomic.columns=3,
                           track.num = i, side = "out")
      add_legend(row_offset = legendreverser[i]-1,
                 ParamsRC = RCircos.Get.Plot.Parameters(),
                 min.value = min.value,
                 max.value = max.value, leglab = target1)
    }
  }
  if(length(targets_inner)>0){

    for (i in 1:length(targets_inner)){
      target1=targets_inner[i]
      tmp2=tmp[!is.infinite(plotdata[,target1]), ]
      if(target1 %in% pvalcols){
        min.value = min.value.pval
        max.value=  max.value.pval
        RC.param <- RCircos.Get.Plot.Parameters()
        RC.param['heatmap.color'] <- "BlackOnly"
        RCircos.Reset.Plot.Parameters(RC.param)

      } else {
        min.value = min(tmp2[,target1], na.rm=T)
        max.value = max(tmp2[,target1], na.rm=T)
        RC.param <- RCircos.Get.Plot.Parameters()
        RC.param['heatmap.color'] <- "BlueWhiteRed"
        RCircos.Reset.Plot.Parameters(RC.param)
      }

      RCircos.Heatmap.Plot(tmp2,
                           data.col = which(colnames(tmp2)==target1),
                           min.value = min.value,
                           max.value = max.value,
                           genomic.columns=3,
                           track.num = i, side = "in")
      add_legend(row_offset = i+length(targets_outer),
                 ParamsRC = RCircos.Get.Plot.Parameters(),
                 min.value = min.value,
                 max.value = max.value, leglab = target1)
      }
  }
  title(main=title)

  if(length(labelsidx)>0){
    RCircos.Gene.Connector.Plot(labdat, genomic.columns=3,track.num=length(targets_outer)+1,side="out")
    RCircos.Gene.Name.Plot(labdat, name.col=4, genomic.columns=3,track.num=length(targets_outer)+2,side="out")
  }
  dev.off()
}


# convert  UCSC track data genomic ranges oject
convertUCSCtoGR = function(tabledat,
                           col.chr = "chrom",
                           col.start = "chromStart",
                           col.end = "chromEnd",
                           col.strand = NULL){

  GRobj = with(tabledat, GRanges(get(col.chr), IRanges(get(col.start), get(col.end))))

  if(length(col.strand)>0){
    strand(GRobj) = tabledat[,col.strand]
  }
  values(GRobj) = tabledat %>% select(-all_of(c(col.chr,col.start,col.end,col.strand)))
  return(GRobj)
}






#get tags matching a genelist ####
checkmatch = function(subject, query){
  return(any(strsplit(subject,",")[[1]] %in% query))
}

gettaglistforgenelist = function(genelist, ddsrangeobject, column="raw_gene"){
  Ranges = as.data.frame(rowRanges(ddsrangeobject))
  tmp=sapply(Ranges[[column]], checkmatch, query=genelist)
  return(which(tmp))
}


# OR plot horizontal ####

# to visualize pval #
convertpvaltostars=function(x){
  sapply(x, function(x){ifelse(x<=0.01, "**", ifelse(x<=0.05, "*", ""))})
}

multiORplot = function(datatoplot=FALSE, Pval = "Pval", Padj = "Padj", SE = "SE", beta="beta", pheno = "pheno"){
  starpval=convertpvaltostars(datatoplot[[Pval]])
  starpval[datatoplot[[Padj]]<0.05]="adj.p**"
  starpval[is.na(datatoplot[[beta]])]="n.a."
  CIUpper = datatoplot[[beta]] +1.96*datatoplot[[SE]]
  CILower = datatoplot[[beta]] -1.96*datatoplot[[SE]]
  xlim=range(c(CIUpper, CILower), na.rm=T)*1.2
  par(mar=c(5,12,5,3))
  betas = datatoplot[[beta]]

  plot(x=betas, y=1:length(betas),
       type="n", panel.first = grid(ny=NA),
       yaxt = "n", ylab="",
       xlim=xlim,
       xlab=expression(paste('log(OR)'%+-%95,"%CI")),
       main=paste(pheno))
  abline(v=0,col="black", lty=3)
  axis(2, at=1:length(betas),
       labels=base::rev(rownames(datatoplot)),
       las=1)
  arrows(x0=CILower, x1=CIUpper, y0=length(betas):1, y1=length(betas):1, col=rainbow(length(betas)), length=0, lwd=2,code = 3)
  points(y=length(betas):1, x=betas, pch=18, col="black")
  betas[is.na(betas)]=0
  text(y=(length(betas):1)+0.5, x=betas, labels=starpval, cex=0.7)
}


# display  functions ####
display_tab = function(df){
  df %>% kbl(digits = 3,align = "l") %>% kable_classic(full_width = F, position = "left") %>% scroll_box(width = "900px", height = paste0(40*nrow(df),"px"), fixed_thead = T)
}

display_tab_simple = function(df){
  df %>% kbl(digits = 3,align = "l") %>% kable_classic(full_width = F)
}



merge.with.order <- function(x,y, ..., sort = T, keep_order)
{
  # this function works just like merge, only that it adds the option to return the merged data.frame ordered by x (1) or by y (2)
  add.id.column.to.data <- function(DATA)
  {
    data.frame(DATA, id... = seq_len(nrow(DATA)))
  }
  # add.id.column.to.data(data.frame(x = rnorm(5), x2 = rnorm(5)))
  order.by.id...and.remove.it <- function(DATA)
  {
    # gets in a data.frame with the "id..." column.  Orders by it and returns it
    if(!any(colnames(DATA)=="id...")) stop("The function order.by.id...and.remove.it only works with data.frame objects which includes the 'id...' order column")

    ss_r <- order(DATA$id...)
    ss_c <- colnames(DATA) != "id..."
    DATA[ss_r, ss_c]
  }

  # tmp <- function(x) x==1; 1	# why we must check what to do if it is missing or not...
  # tmp()

  if(!missing(keep_order))
  {
    if(keep_order == 1) return(order.by.id...and.remove.it(merge(x=add.id.column.to.data(x),y=y,..., sort = FALSE)))
    if(keep_order == 2) return(order.by.id...and.remove.it(merge(x=x,y=add.id.column.to.data(y),..., sort = FALSE)))
    # if you didn't get "return" by now - issue a warning.
    warning("The function merge.with.order only accepts NULL/1/2 values for the keep_order variable")
  } else {return(merge(x=x,y=y,..., sort = sort))}
}



cerebroScale2<-
  function (x, clamp, divData, center_zero)
  {
    xmed <- median(x, na.rm = TRUE)
    xmad <- mad(x, constant = 1, na.rm = TRUE)
    xmin <- min(x, na.rm = TRUE)
    xmax <- max(x, na.rm = TRUE)
    avoidClamp <- max(abs(xmed - xmin), abs(xmed - xmax))/xmad
    fillMatrix <- x
    if (is.null(clamp)) {
      clamp <- avoidClamp + 1
    }
    outlrs <- clamp * xmad
    if (clamp <= 0)
      stop("clamp must be >0")
    pctOL <- round(length(which(x[!is.na(x)] <= (xmed - (outlrs)) |
                                  x[!is.na(x)] >= (xmed + (outlrs))))/length(x[!is.na(x)]) *
                     100, 2)
    if (pctOL > 0) {
      warning(paste("The clamp value of ", clamp, " will clamp ",
                    pctOL, "% of input values (outliers) to the min or max of the scaled range.",
                    sep = ""))
    }
    if (divData == TRUE & center_zero == FALSE) {
      abvMed <- x[x >= xmed & x <= (xmed + outlrs) & !is.na(x)]
      belMed <- x[x <= xmed & x >= (xmed - outlrs) & !is.na(x)]
      if (length(which(!is.na(x)))%%2 == 0) {
        rightsc <- rescale(c(xmed, abvMed), c(0.5, 1))[-1]
        fillMatrix[x >= xmed & x <= (xmed + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(c(xmed, belMed), c(-1, 0))[-1]
        fillMatrix[x <= xmed & x >= (xmed - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
      if ((length(which(!is.na(x))))%%2 == 1) {
        rightsc <- rescale(abvMed, c(-1, 1))
        fillMatrix[x >= xmed & x <= (xmed + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(belMed, c(-1, 0))
        fillMatrix[x <= xmed & x >= (xmed - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
    }
    if (divData == TRUE & center_zero == TRUE) {
      abvMed <- x[x >= -1 & x <= (-1 + outlrs) & !is.na(x)]
      belMed <- x[x <= -1 & x >= (-1 - outlrs) & !is.na(x)]
      if (length(which(!is.na(x)))%%2 == 0) {
        rightsc <- rescale(c(-1, abvMed), c(-1, 1))[-1]
        fillMatrix[x >= -1 & x <= (-1 + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(c(-1, belMed), c(-1, 0))[-1]
        fillMatrix[x <= -1 & x >= (-1 - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (-1 - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (-1 + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
      if ((length(which(!is.na(x))))%%2 == 1) {
        rightsc <- rescale(abvMed, c(-1, 1))
        fillMatrix[x >= -1 & x <= (-1 + outlrs) & !is.na(x)] <- rightsc
        leftsc <- rescale(belMed, c(-1, 1))
        fillMatrix[x <= -1 & x >= (-1 - outlrs) & !is.na(x)] <- leftsc
        fillMatrix[x < (-1 - outlrs) & !is.na(x)] <- -1
        fillMatrix[x > (-1 + outlrs) & !is.na(x)] <- 1
        xScaled <- fillMatrix
      }
    }
    if (divData == FALSE) {
      nonoutlrs <- x[x >= (xmed - outlrs) & x <= (xmed + outlrs) &
                       !is.na(x)]
      xsc <- rescale(nonoutlrs, c(-1, 1))
      fillMatrix[x >= (xmed - outlrs) & x <= (xmed + outlrs) &
                   !is.na(x)] <- xsc
      fillMatrix[x < (xmed - outlrs) & !is.na(x)] <- -1
      fillMatrix[x > (xmed + outlrs) & !is.na(x)] <- 1
      xScaled <- fillMatrix
    }
    return(xScaled)
  }
