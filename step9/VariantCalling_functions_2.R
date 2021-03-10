
Parse_VCF_File <- function(path) {

	require(tidyverse)

	vcf <- read.table(path,sep="\t")
	vcf$V1 <- as.character(vcf$V1)
	vcf$V2 <- as.numeric(vcf$V2)
	vcf$V3 <- as.character(vcf$V3)
	vcf$V4 <- as.character(vcf$V4)
	vcf$V5 <- as.character(vcf$V5)
	vcf$V6 <- as.character(vcf$V6)
	vcf$V7 <- as.character(vcf$V7)
	vcf$V8 <- as.character(vcf$V8)
	vcf$V9 <- as.character(vcf$V9)
	vcf$V10 <- as.character(vcf$V10)
	vcf$V11 <- as.character(vcf$V11)
	
	#mutations_info <- as.character(mutations$V8)
	
	return(vcf)
}

Reformat_Annotated_Aggregated_VCF <- function(listcl_mutations) {
  
  # ONLY Returns SNPS!!!
  # Removes INDELS
  
  mutations <- listcl_mutations %>% 
    mutate(ref_len = nchar(V4), alt_len = nchar(V5)) %>%
    filter(ref_len == 1 & alt_len == 1) %>% 
    separate(8,c('CONTQ','DP','ECNT','GERMQ','MBQ','MFRL','MMQ','MPOS','NALOD','NLOD','POPAF','SAAF','SAP','TLOD','ANN'),sep=';') %>%
    separate(22, as.character(seq(1,15,1)), sep='\\|')
  colnames(mutations) <- c('Chr','POS','ID','REF','ALT','QUAL','FILTER',
                           'CONTQ','DP','ECNT','GERMQ','MBQ','MFRL','MMQ','MPOS','NALOD','NLOD','POPAF','SAAF','SAP','TLOD','ANN',
                           'VAR_TYPE','MOD','GENE','ENSEMBL_GENE_ID','READ_TYPE','ENST_ID','FUNCTION','ratio','NUC_CHANGE','AA_CHANGE',
                           'ratio1','ratio2','ratio3','number','stat1','stat2','stat3','donor','bc','ref_len','alt_len')
  mutations <- mutations %>% mutate(CONTQ = parse_number(CONTQ), DP = parse_number(DP), ECNT = parse_number(ECNT), GERMQ = sub(".*=", "", GERMQ),
                                         MBQ = sub(".*=", "", MBQ), MFRL = sub(".*=", "", MFRL), MMQ = sub(".*=", "", MMQ), MPOS = parse_number(MPOS), 
                                         NALOD = parse_number(NALOD), NLOD = parse_number(NLOD), POPAF = parse_number(POPAF), SAAF = sub(".*=", "", SAAF), 
                                         SAP = sub(".*=", "", SAP), TLOD = parse_number(TLOD), ANN = sub(".*=", "", ANN))
  
  
  return(mutations)
  
  # mutations_info <- as.character(mutations$V8)
  # mutations_split <- as_tibble(mutations_info) #
  # mut_split <- separate(mutations_split, value, into = as.character(seq(1,20,1)), sep = ';')
  # 
  # split_info <- strsplit(mutations_info, ';', fixed = T)
  # data.frame(do.call(rbind, split_info))
  # split_in <- as.data.frame(split_info)
  # for (i in c(1:length(split_info))
  # tail(split_info[[1]], n =1)
  # 
}

score.mutations <- function(mutations, point, read) {
  
  require(tidyverse) 
  sam = unique(mutations$donor)
  
  # Count number of each mutation seen in Sample:
    count.mutations <- mutations %>%
      group_by(Chr,POS,REF,ALT) %>%
      summarise(n())
    
    # Generate Running Filter Dataframe:
    
    mutations.filtered <- mutations
    
    # Filter mutations for candidates for correction:
    filter.mutect <- mutations %>%
      mutate(mutect_filter = ifelse(TLOD >= 5.3 & DP > 10 & ECNT > 1, 'pass', 'fail')) %>%
      filter(mutect_filter == 'pass') %>%
      dplyr::select(bc,Chr,POS,ALT,mutect_filter)
    mutations.filtered <- mutations.filtered %>% left_join(filter.mutect, by = c('bc','Chr','POS','ALT'))
   
   print('passed mutect filter')
   
    # For Candidates for Correction:
    # Transform Point Mutations Read Data:
    # Join Point Mutations and Read Information:
    colnames(point) <- c('read','zero','flag','Chr','p','ALT','qual','POS','i','j')
    point <- point %>% filter(ALT != '.')
    
    read <- read %>% separate(X1, into = c('read','bc','umi'), sep = '___') %>% 
      mutate(bc = substr(bc,6,24), umi = substr(umi,6,16)) %>%
      distinct() %>% filter(str_length(umi) == 10) %>% dplyr::select(read, bc, umi)
  
    point.reads <- point %>%
      inner_join(read, by = 'read') %>%
      mutate(donor = sam)
      
    print('mutations reformatted')
    
    # FILTER:
    
    # Filter Point Mutations using the candidate list:
    # (Just making sure the candidates are the same as original data)
    candidate_positions <- mutations.filtered %>% filter(mutect_filter == 'pass') %>% dplyr::select(bc, Chr, POS)
    point.reads.filter <- inner_join(point.reads, candidate_positions, by = c('bc','Chr','POS'))
  
    # Remove any point mutations with less than 2 reads in a UMI:
    # (These can't be corrected)
    filter.reads.in.umi <- point.reads.filter %>% group_by(bc,Chr,POS,umi) %>%
      summarise(reads_in_umi = n()) %>% mutate(reads_in_umi_filter = ifelse(reads_in_umi > 1, 'pass','fail'), verbose = F) %>% #inner_join(point.reads, by = c('bc','umi'))
      filter(reads_in_umi_filter == 'pass') 
    joining_reads <- filter.reads.in.umi %>%
      dplyr::select(bc,Chr,POS,reads_in_umi_filter) %>% 
      distinct()
    mutations.filtered <- mutations.filtered %>% left_join(joining_reads, by = c('bc','Chr','POS'))
    
    print('read_number filter complete')
    
    # Remove point mutations with fewer than 50% of supporting reads in a UMI:
    # (These are likely errors)
    read.umi.fraction.filter <- filter.reads.in.umi %>% 
      ungroup() %>% distinct() %>% dplyr::select(bc,umi) %>% inner_join(point.reads, by = c('bc','umi')) %>%
      group_by(bc,Chr,POS,umi,ALT) %>%
      summarise(num = n(), verbose = F) %>%
      group_by(bc,Chr,POS,umi) %>%
      mutate(umi_fraction = num / sum(num)) %>% # Fraction of each base in UMI
      mutate(umi_fraction_filter = ifelse(umi_fraction > 0.5, 'pass', 'fail')) %>% # Keep bases with more than 50% reads in UMI
      filter(umi_fraction_filter == 'pass') %>%
      ungroup() 
    joining.umi.fraction.filter <- read.umi.fraction.filter %>%
      dplyr::select(bc,Chr,POS,umi_fraction_filter) %>%
      distinct()
    mutations.filtered <- mutations.filtered %>% left_join(joining.umi.fraction.filter, by = c('bc','Chr','POS'))
    
    print('umi_fraction filter complete')
    
    # Remove point mutations with only 1 UMI:
    # (These are possible errors - This filter likely does not improve signal/noise)
    umi.number.filter <- read.umi.fraction.filter %>% 
      dplyr::select(bc,Chr,POS,umi) %>% 
      distinct() %>%
      group_by(bc,Chr,POS) %>% 
      summarize(umi_count = n()) %>%
      mutate(umi_count_filter = ifelse(umi_count > 1, 'pass', 'fail')) %>%
      filter(umi_count_filter == 'pass') %>% 
      ungroup()
    joining.umi.number.filter <- umi.number.filter %>%
      dplyr::select(bc,Chr,POS,umi_count_filter) %>%
      distinct()
    mutations.filtered <- mutations.filtered %>% left_join(joining.umi.number.filter, by = c('bc','Chr','POS'))
    
    print('umi_number filter complete')
    
    # Remove point mutations with conflicting UMI's:
    conflict.umi.fraction.filter <- umi.number.filter %>% 
      dplyr::select(bc,Chr,POS) %>% distinct() %>% inner_join(read.umi.fraction.filter, by = c('bc','Chr','POS')) %>%
      dplyr::select(bc,Chr,POS,umi,ALT) %>% distinct() %>% 
      group_by(bc,Chr,POS,ALT) %>% summarize(number_umi_ALT = n()) %>%
      group_by(bc,Chr,POS) %>% mutate(conflict_umi_fraction = number_umi_ALT / sum(number_umi_ALT)) %>%
      mutate(conflict_umi_filter = ifelse(conflict_umi_fraction > 0.2, 'pass', 'fail')) %>%
      filter(conflict_umi_filter == 'pass') %>%
      ungroup()
    joining.conflict.umi.fraction.filter <- conflict.umi.fraction.filter %>%
      dplyr::select(bc,Chr,POS,conflict_umi_filter) %>%
      distinct()
    mutations.filtered <- mutations.filtered %>% left_join(joining.conflict.umi.fraction.filter, by = c('bc','Chr','POS'))
    
    print('umi_conflict filter')
    
    #
    return(mutations.filtered)
  }


replace_low_qual_bases <- function(path, Q_threshold, n_cores) {

	require(doParallel)
	require(parallel)
	require(foreach)
  	require(tidyverse)
  	require(Rsamtools)
  	library(stringr)

  # Sample Values:
  #path <- path
  #Q_threshold <- 25
  #n_cores <- 29
	
  cl <- makeCluster(n_cores)
	registerDoParallel(cl, cores = n_cores)

	bam <- BamFile(path)
	yieldSize(bam) <- 100
	#ScanBamParam(tag=c("NM", "MD:Z", "MC:Z", "AS:i", "XS:i", "SA:Z", ), what="flag")
	open(bam)
	reads <- scanBam(bam)

	low_scores <- phred[phred$q_score < Q_threshold,]$symbol
	vals <- paste(low_scores, collapse = '')
	regex <- paste('[',vals,']',sep= '')

	positions <- list()
	clusterExport(cl = cl, c("reads","positions"), envir = environment())

	for (j in c(1:length(reads))) {
	  foreach (i=1:length(reads[[j]]$seq), .packages = c('tidyverse', 'Rsamtools', 'stringr')) %dopar% {
	  #for (i in c(1:length(reads[[j]]$seq))) {
	  	
	    
	    locations <- str_locate_all(reads[[j]]$qual[[i]], regex)[[1]][,1]
	  	print(locations)
	  	len <- length(locations)
	  	positions[[as.character(i)]] <- locations
	  
		#if (len > 0) {
			#reads[[j]]$seq[[i]] <- replaceLetterAt(reads[[j]]$seq[[i]], locations, rep('N',len))
		#}
		}
	}

#stopCluster(cl)
final <- list(reads, positions)
return(final)
}
