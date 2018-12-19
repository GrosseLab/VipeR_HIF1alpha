# conda env export > r35.yaml

library("GenomicRanges")
library("GenomicFeatures")
library("Rsamtools")
library("seqinr")
library("rtracklayer")
library("data.table")
library("BSgenome")


# GRCh38.82 ---------------------------------------------------------------
  anno.file <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38/Homo_sapiens.GRCh38.82.gtf'
  anno.format <- "gtf"
  fa.file <-  '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
  outputFile <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38/Homo_sapiens.GRCh38.82.EXON.fa'
  
  indexFa(fa.file)
  fa.seqs <- FaFile(fa.file)
  
  gr.db <- makeTxDbFromGFF(anno.file,format = anno.format)
  exon.list <- exonsBy(gr.db, "tx",use.names=T)
  
  # tmp <- purrr::map(names(exon.list),function(x) unlist(getSeq(fa.seqs,exon.list[[x]],as.character=TRUE)) )
  library(furrr)
  plan(multiprocess, workers = 15)
  tmp <- furrr::future_map(names(exon.list),function(x) unlist(getSeq(fa.seqs,exon.list[[x]],as.character=TRUE)),.progress = T)
  names(tmp) <- names(exon.list)
  tmp2 <- unlist(tmp, use.names=TRUE)
  writeXStringSet(DNAStringSet(tmp2), outputFile)

# GRCh38.94 ---------------------------------------------------------------
  
  anno.file <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38v94/Homo_sapiens.GRCh38.94.gtf'
  anno.file.out <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38v94/Homo_sapiens.GRCh38.94.fixed.gtf'
  anno.format <- "gtf"
  fa.file <-  '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38v94/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
  outputFile <- '/home/adsvy/GitHubRepo/SnakeWF_HIF/references/hg38v94/Homo_sapiens.GRCh38.94.fixed.EXON.fa'
  
  ### 1.
  ### problem : anno$gene_id, anno$gene_version can not be parsed correctly by STAR and salmon 
  ### solution: merge anno$gene_id and anno$gene_version
  anno <- import.gff(con=anno.file ,format=anno.format ) 
  
  anno$gene_IDV <-  paste0(anno$gene_id,'.',anno$gene_version)
  anno$transcript_IDV <-  paste0(anno$transcript_id,'.',anno$transcript_version)
  anno$protein_IDV <-  paste0(anno$protein_id,'.',anno$gene_version)
    
  anno2 <- anno
  elementMetadata(anno2) <- elementMetadata(anno2)[,setdiff(names(elementMetadata(anno2)),c('gene_id','gene_version','transcript_id','transcript_version','protein_id','gene_version') ) ]
  names(elementMetadata(anno2))[which(names(elementMetadata(anno2)) == 'gene_IDV' )] <- 'gene_id'
  names(elementMetadata(anno2))[which(names(elementMetadata(anno2)) == 'transcript_IDV' )] <- 'transcript_id'
  names(elementMetadata(anno2))[which(names(elementMetadata(anno2)) == 'protein_IDV' )] <- 'protein_id'
  export.gff(object=anno2,con=anno.file.out ,format=anno.format ) 
  
  ### 2.
  indexFa(fa.file)
  fa.seqs <- FaFile(fa.file)

  gr.db <- makeTxDbFromGFF(anno.file.out,format = anno.format)
  exon.list <- exonsBy(gr.db, "tx",use.names=T)
  
  # tmp <- purrr::map(names(exon.list),function(x) unlist(getSeq(fa.seqs,exon.list[[x]],as.character=TRUE)) )
  library(furrr)
  plan(multiprocess, workers = 15)
  tmp <- furrr::future_map(names(exon.list),function(x) unlist(getSeq(fa.seqs,exon.list[[x]],as.character=TRUE)),.progress = T)
  names(tmp) <- names(exon.list)
  tmp2 <- unlist(tmp, use.names=TRUE)
  writeXStringSet(DNAStringSet(tmp2), outputFile)


# ### still don't work ----------------------------------------------------

  ### still don't work
  # library("BSgenome.Hsapiens.UCSC.hg38") 
  # olaps<- as.data.table(anno2)
  # names(olaps)[1] <- 'names'
  # seqnames(BSgenome.Hsapiens.UCSC.hg38) <- as.character( sapply(names(BSgenome.Hsapiens.UCSC.hg38),function(x) stringr::str_remove(x,'chrUn_')))
  # seqnames(BSgenome.Hsapiens.UCSC.hg38) <- as.character( sapply(names(BSgenome.Hsapiens.UCSC.hg38),function(x) stringr::str_remove(x,'chr')))
  # seqnames(BSgenome.Hsapiens.UCSC.hg38) <- as.character( sapply(names(BSgenome.Hsapiens.UCSC.hg38),function(x) stringr::str_replace(x,'M','MT')))
  # seqnames(BSgenome.Hsapiens.UCSC.hg38) <- as.character( sapply(names(BSgenome.Hsapiens.UCSC.hg38),function(x) stringr::str_replace(x,'MTTT','MT')))
  # 
  # 
  # keepBSgenomeSequences <- function(genome, seqnames)
  # {
  #   stopifnot(all(seqnames %in% seqnames(genome)))
  #   genome@user_seqnames <- setNames(seqnames, seqnames)
  #   genome@seqinfo <- genome@seqinfo[seqnames]
  #   genome
  # }
  # genome <- keepBSgenomeSequences(BSgenome.Hsapiens.UCSC.hg38, 
  #                               seqlevels(BSgenome.Hsapiens.UCSC.hg38)[ (!sapply(names(BSgenome.Hsapiens.UCSC.hg38),function(x) stringr::str_detect(x,'alt') ) &
  #                                        !sapply(names(BSgenome.Hsapiens.UCSC.hg38),function(x) stringr::str_detect(x,'random') ) )
  #                                        ])
  # 
  # genome2 <- keepBSgenomeSequences(genome,  intersect(seqlevels(genome),levels(olaps$names))) 
  # olaps2 <- olaps[olaps$names %in%  intersect(seqlevels(genome),levels(olaps$names)) ,]
  # 
  # genomerngs = as(seqinfo(BSgenome.Hsapiens.UCSC.hg38), "GRanges")
  # genomerngs[genomerngs$] 
  # 
  # seq = BSgenome::getSeq(genome2, olaps2)
  
  
  
  
  
  
  ### still don't work
  # exon.seqs <- getSeq(fa.seqs, unlist(exon.list, use.names=TRUE),as.character=TRUE)
  # names(exon.seqs) <- names(unlist(exon.list, use.names=TRUE))
  # writeXStringSet(exon.seqs, outputFile)
  
  ##########################################
  testxxx <- function(){
  
  	# features <- binGenome(gr.db)
  	# trans.list <- transcriptsBy(gr.db, by=c("gene", "exon", "cds"), use.names=FALSE)
  	# trans.seqs <- getSeq(fa.seqs, unlist(trans.list, use.names=TRUE), as.character=TRUE)
  
  	if( packageDescription("GenomicRanges")$Version == "1.22.4"){
  	  exon.elt <- rep(names(trans.list), GenomicFeatures::elementLengths(trans.list)) ## elementLengths is the old version of elementNROWS()
  	}else{
  	  # trans.elt <- rep(names(trans.list), elementNROWS(trans.list))   
  	  exon.elt <- rep(names(exon.seqs), elementNROWS(exon.seqs))   
  	}
  
  	seq.list <- lapply(split(exon.seqs, exon.elt),paste,collapse="")
  
  	write.fasta(seq.list,file.out = outputFile,names = names(seq.list),nbchar = max(sapply(seq.list,nchar)),as.string = T)
  
  
  	#######
  
  	seq_gtf = function(gtf, seqs, feature='transcript', exononly=TRUE, 
  	    idfield='transcript_id', attrsep="; "){
  
  	    feature = match.arg(feature, c('transcript', 'exon'))
  
  	    gtfClasses = c("character", "character", "character", "integer", 
  	        "integer", "character", "character", "character", "character")
  	    if(is.character(gtf)){
  	        # read transcript structure from file:s
  	        gtf_dat = read.table(gtf, sep="\t", as.is=TRUE, quote="", header=FALSE, 
  	            comment.char="#", nrows= -1, colClasses=gtfClasses)
  	    } else if(is.data.frame(gtf)){
  	        # do what we can to check whether gtf really does represent a 
  	        # canonical GTF
  	        stopifnot(ncol(gtf) == 9)
  	        if(!all(unlist(lapply(gtf, class)) == gtfClasses)){
  	            stop("one or more columns of gtf have the wrong class")
  	        }
  	        gtf_dat = gtf
  	        rm(gtf)
  	    } else {
  	        stop("gtf must be a file path or a data frame")
  	    }
  
  	    colnames(gtf_dat) = c("seqname", "source", "feature", "start", "end",
  	        "score", "strand", "frame", "attributes")
  	    stopifnot(!any(is.na(gtf_dat$start)), !any(is.na(gtf_dat$end)))
  
  	    if(exononly){
  	        gtf_dat = gtf_dat[gtf_dat[,3]=="exon",]
  	    }
  
  	    # makes sure all chromosomes are present:
  	    chrs = unique(gtf_dat$seqname)
  	    if(is.character(seqs)){
  	        fafiles = list.files(seqs)
  	        lookingFor = paste0(chrs, '.fa')
  	    } else {
  	        fafiles = names(seqs)
  	        lookingFor = chrs
  	    }
  	    if(!(all(lookingFor %in% fafiles))){
  	        stop("all chromosomes in gtf must have corresponding sequences in seqs")
  	    }
  	    
  	    seqlist = lapply(chrs, function(chr){
  	        dftmp = gtf_dat[gtf_dat[,1] == chr,]
  	        if(is.character(seqs)){
  	            fullseq = readDNAStringSet(paste0(seqs, '/', chr, '.fa'))
  	        } else {
  	            fullseq = seqs[which(names(seqs) == chr)]
  	        }
  	        if(feature == 'exon'){
  	            dftmp = dftmp[!duplicated(dftmp[,c(1,4,5,7)]),] #unique exons
  	        }
  	        these_seqs = subseq(rep(fullseq, times=nrow(dftmp)), 
  	            start=dftmp$start, end=dftmp$end)
  	        if(feature == 'transcript'){
  	            names(these_seqs) = getAttributeField(dftmp$attributes, idfield, 
  	                attrsep=attrsep)
  	            if(substr(names(these_seqs)[1],1,1) == '"'){
  	                x = names(these_seqs)
  	                names(these_seqs) = substr(x, 2, nchar(x)-1)
  	            }
  	        }else{
  	            names(these_seqs) = paste0(dftmp[,1], ':', dftmp[,4], '-', 
  	                dftmp[,5], '(', dftmp[,7], ')')
  	        }
  	        revstrand = which(dftmp$strand == '-')
  	        these_seqs[revstrand] = reverseComplement(these_seqs[revstrand])
  	        these_seqs
  	    })
  
  	    full_list = do.call(c, seqlist)
  
  	    if(feature == 'exon'){
  	        return(full_list)
  	    }else{
  	        split_list = split(full_list, names(full_list))
  	        return(DNAStringSet(lapply(split_list, unlist)))
  	    }
  	}
  	 library(Biostrings)
  	bioSeq <- seq_gtf(anno.file, fa.file, feature = "transcript", exononly = TRUE,idfield = "transcript_id", attrsep = "; ")
  
  }
