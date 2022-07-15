library(dplyr)
library(readr)
library(RMySQL)
library(parallel)
library(GenomicRanges)
library(ggplot2)
options(useFancyQuotes = FALSE)
dbGroup <- 'intsites_miseq'

if(! file.exists('repeatMaskerTable.tsv.gz')) system('wget https://microb120.med.upenn.edu/data/export/everett/ViiV_HIV/repeatMaskerTable.tsv.gz')

# Read in repeatMasker result and convert into a GRange object.
repeats <- read_tsv('./repeatMaskerTable.tsv.gz')
repeats$strand <- ifelse(repeats$strand == 'C', '-', '+')
repeats <- makeGRangesFromDataFrame(repeats, keep.extra.columns = TRUE, seqnames.field = 'query_seq', start.field = 'query_start', end.field = 'query_end')

# Read in intSites.
# There are duplicate sites across samples: chr3+93470477, chr1-145287843, chr19+19343784, chr4+70003445
sites <- read_tsv('intSites_20220708.tsv')

groupings <- read.table('groupings.csv', header = FALSE, sep = ';')

sites <- subset(sites, GTSP %in% unlist(strsplit(groupings$V2, ',')))


# Expand intSites +/- 5 NTs for comparison to multi-hit positions since mult-hit
# positions are not standardized like the final intSite positions.
sites.expanded <- sites
sites.expanded$start <- sites.expanded$start - 5
sites.expanded$end <- sites.expanded$end + 5
sites.expanded <- makeGRangesFromDataFrame(sites.expanded, keep.extra.columns = TRUE)


# Find overlaps between intSites and repeatMaskMaker result.
sites <- makeGRangesFromDataFrame(sites, keep.extra.columns = TRUE)
o <- GenomicRanges::findOverlaps(sites, repeats, ignore.strand = TRUE)

sites.repeatHits <- dplyr::distinct(bind_rows(lapply(unique(queryHits(o)), function(x){
                      tibble(posid = sites[x]$posid, 
                             class = paste0(unique(repeats[subjectHits(o[which(queryHits(o) == x)])]$repeat_class), collapse = '; '),
                             seq = paste0(as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, 
                                                                        as.character(seqnames(repeats[subjectHits(o[which(queryHits(o) == x)])])), 
                                                                        start=start(repeats[subjectHits(o[which(queryHits(o) == x)])]),  
                                                                        end=end(repeats[subjectHits(o[which(queryHits(o) == x)])]))), collapse = '; '))
                    })))

# Retrieve multi-hits for the samples in the intSite data set.
sql <- paste0("select samples.sampleName, samples.refGenome, multihitpositions.multihitID, ",
              "multihitpositions.chr, multihitpositions.strand,  multihitpositions.position, ",
              "multihitlengths.length from multihitlengths left join multihitpositions on ",
              "multihitpositions.multihitID = multihitlengths.multihitID left join samples on ",
              "samples.sampleID = multihitpositions.sampleID where sampleName like ",
              paste0(sQuote(paste0(unique(sites$GTSP), '-%')), collapse = ' or sampleName like '))

dbConn  <- dbConnect(MySQL(), group = dbGroup)
sites.multi <- unique(dbGetQuery(dbConn, sql))
dbDisconnect(dbConn)


# Go through the replicate level multi hits and remove positions that overlap with 
# well established sites allowing for a +/- 5 NT difference. Cycle through remaining 
# positions and merge with a +/- 5 NT difference while tallying abundances.

cluster <- makeCluster(15)
clusterExport(cluster, c('sites.expanded'))

# A single abundance value is associated with a cluster.
multiHits  <- bind_rows(parLapply(cluster, split(sites.multi, sites.multi$multihitID), function(x){
  library(dplyr)
  library(GenomicRanges)
  
  # Find overlaps between positions in this multi-hit group with the established sites that
  # were expanded +/- 5 NTs. 
  q <- makeGRangesFromDataFrame(x, start.field = 'position', end.field = 'position', keep.extra.columns = TRUE)
  o <- GenomicRanges::findOverlaps(q, sites.expanded, ignore.strand = FALSE)
  
  # If a position in the multi-hit is found in the non-multihit data then skip the multi-hit.
  if(length(o) > 0) return(tibble())
  
  # Spread out the positions of the sites in this multi-hit.
  x$start <- x$position - 5
  x$end <- x$position + 5
  
  # Merge closely spaced sites.
  q  <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  q2 <- GenomicRanges::reduce(q, with.revmap = TRUE)
  
  # Tally the abundances of merged sites.
  bind_rows(lapply(1:length(q2), function(a){
         tibble(sample = x$sampleName[1], 
                multiHitID = x$multihitID[1],
                seqnames = unique(as.character(seqnames(q2[a]))), 
                strand = unique(as.character(strand(q2[a]))), 
                position = start(q2[a])+5, 
                estAbund = n_distinct(q[unlist(q2[a]$revmap)]$length))
  }))
}))


# Expand multi-hit meta data.
multiHits$GTSP  <- stringr::str_extract(multiHits$sample, 'GTSP\\d+')
multiHits$posid <- paste0(multiHits$seqnames, multiHits$strand, multiHits$position)


# Here we merge multi-hit clusters if they share one or more posids +/- 5 NTs
# to create sample-level multihits rather than techincal replicate-level multihits.

multiHits2 <- bind_rows(parLapply(cluster, split(multiHits, multiHits$GTSP), function(x){
  library(GenomicRanges)
  library(dplyr)
  
  x$start <- x$position - 5
  x$end <- x$position + 5
  
  a <- GenomicRanges::makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  r <- GenomicRanges::reduce(a, with.revmap = TRUE)
  
  # Tally the abundances of merged sites.
  # Most entries will be singletons (no merger).
  bind_rows(lapply(1:length(r), function(k){
    tibble(sample = x$GTSP[1], 
           multiHitID = paste0(unique(a[unlist(r[k]$revmap)]$multiHitID), collapse = ';'),
           seqnames = unique(as.character(seqnames(r[k]))), 
           strand = unique(as.character(strand(r[k]))), 
           position = start(r[k])+5, 
           estAbund = max(a[unlist(r[k]$revmap)]$estAbund))
  }))
}))



# Visualize members of multi-hits.
multiHits2$posid <- paste0(multiHits2$seqnames, multiHits2$strand, multiHits2$position)

d <- group_by(multiHits2, multiHitID) %>% 
     summarise(n = n_distinct(posid)) %>% 
     ungroup() %>% arrange(n) %>%
     mutate(multiHitID2 = factor(multiHitID, levels = unique(multiHitID)))

multiHitMemberPlot <- 
  ggplot(d, aes(multiHitID2, n)) + 
  geom_col() +
  labs(x = 'Multi-hit clusters', y = 'Positions within cluster') +
  ggtitle(paste0(n_distinct(d$multiHitID), ' Multi-hit clusters')) +
  theme(axis.text.x=element_blank())


classifyRepeat <- function(x){
  if(is.na(x) | is.null(x)){
    return('None')
  } else if(grepl('Alu', x, ignore.case = TRUE)){
    return('SINE/Alu')
  } else if(grepl('SINE', x, ignore.case = TRUE)){
    return('SINE/Other')
  } else if(grepl('LINE', x, ignore.case = TRUE)){
    return('LINE')
  } else if(grepl('LTR', x, ignore.case = TRUE)){
    return('LTR')
  } else if(grepl('Satellite', x, ignore.case = TRUE) | grepl('Satellite', x, ignore.case = TRUE)){
    return('Satellite/centr')
  } else if(grepl('complexity', x, ignore.case = TRUE)){
    return('Low complexity')
  } else if(grepl('DNA', x, ignore.case = TRUE)){
    return('DNA')
  } else {
    return('Other')
  }
}

# Determine the repeat classes of sites in multi-hit clusters.
# Not practical to use parLapply since we would have to export the repeats 
# object to all worker nodes.

t <- n_distinct(multiHits2$posid)
n <- 1

multiHits2 <- bind_rows(lapply(split(multiHits2, multiHits2$posid), function(x){
  message(n, '/', t); n <<- n + 1
  
  o <- GenomicRanges::findOverlaps(makeGRangesFromDataFrame(x[1,], 
                                                            start.field = 'position', 
                                                            end.field = 'position'),
                                   repeats, ignore.strand = TRUE)
  if(length(o) > 0){
    x$class <- classifyRepeat(paste0(repeats[subjectHits(o)]$repeat_class, collapse = ';'))
  } else {
    x$class <- 'None'
  }
  
  x
}))

multiHits2 <- bind_rows(lapply(split(multiHits2, multiHits2$multiHitID), function(x){
  t <- data.frame(table(x$class))
  t$p <- t$Freq / nrow(x)
  t <- arrange(t, desc(p))
  
  x$class2 <- 'Not clear'
  if(t[1,]$p >= 0.90) x$class2 <- t[1,]$Var1
  x
}))

save.image('image1.RData')


# Add repeat annotations to clearly called sites.
sites <- left_join(data.frame(sites), sites.repeatHits, by = 'posid')
sites$multiHitID <- 'None'

# Prep multi-hit clusters as pseudo-sites.
a <- select(multiHits2, sample, multiHitID, posid, class, class2)
a <- a[!duplicated(a$multiHitID), ]
a$posid <- 'multiHit'

# Prep clearly called sites.
b <- select(sites, GTSP, multiHitID, posid, class) %>% dplyr::rename(sample = GTSP)
b$class2 <- unlist(lapply(b$class, classifyRepeat))

tab <- bind_rows(a, b)

tab[which(tab$class2 == 'Not clear'),]$class2 <- 'Other'


tab$group <- 'None'
tab <- bind_rows(lapply(split(groupings, groupings$V1), function(x){
         o <- subset(tab, sample %in% unlist(strsplit(x$V2, ',')))
         o$group <- x$V1
         o
       }))

classOrder <- c('SINE/Alu', 'SINE/Other', 'LINE', 'LTR', 'Satellite/centr', 'Low complexity', 'DNA', 'Other', 'None')

o <- lapply(split(tab, tab$group), function(x){
       o <- data.frame(table(x$class2))
       m <- classOrder[! classOrder %in% o$Var1]
       if(length(m) > 0) o <- bind_rows(o, data.frame(Var1 = m, Freq = 0))
       o <- o[match(classOrder, o$Var1),]
       names(o) <- c('repeatClass', x$group[1])
       o
     })

countsTable <- bind_cols(lapply(o, '[', 2))
row.names(countsTable) <- classOrder
openxlsx::write.xlsx(countsTable, 'countTable.xlsx', rowNames = TRUE)



o <- lapply(split(tab, tab$group), function(x){
  o <- data.frame(table(x$class2))
  m <- classOrder[! classOrder %in% o$Var1]
  if(length(m) > 0) o <- bind_rows(o, data.frame(Var1 = m, Freq = 0))
  o <- o[match(classOrder, o$Var1),]
  o$Freq <- sprintf("%.1f%%", (o$Freq / nrow(x))*100)
  names(o) <- c('repeatClass', x$group[1])
  o
})

freqTable <- bind_cols(lapply(o, '[', 2))
row.names(freqTable) <- classOrder
openxlsx::write.xlsx(freqTable, 'freqTable.xlsx', rowNames = TRUE)

tab2 <- left_join(tab, select(sites.repeatHits, posid, seq), by = 'posid') 
tab2$seqLen <- nchar(tab2$seq)
tab2 <- select(tab2, sample, multiHitID, posid, class, class2, group, seqLen, seq)
openxlsx::write.xlsx(tab2[rev(1:nrow(tab2)),], file = 'Marquis_repeat_site_annotaions.xlsx')






