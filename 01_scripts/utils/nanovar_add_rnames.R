#setwd("/mnt/ibis/lbernatchez/users/lalec31/RDC_Romaine/02_long_reads/SVs_long_reads")

argv <- commandArgs(T)
VCF <- argv[1]
READS_SUP <- argv[2]
OUT <- argv[3]


# 1. Import ---------------------------------------------------------------
vcf <- data.table::as.data.table(read.table(VCF, comment.char = "#", stringsAsFactors = FALSE,
                                            col.names = c('CHROM', 'POS', 'ID', 'END')))

reads <- data.table::as.data.table(read.table(READS_SUP, comment.char = "#", stringsAsFactors = FALSE, header = TRUE,
                                              col.names = c('SV_ID', 'READS'), sep = "\t"))


# 2. Format read names to remove index ------------------------------------
reads$reads_reformat <- 
  sapply(X = reads$READS, FUN = function(x) {
    gsub(pattern = "~[0-9A-Za-z]+,", replacement = ",", x)
    
  }
  )

reads$reads_reformat <- 
  sapply(X = reads$reads_reformat, FUN = function(x) {
    gsub(pattern = "~[0-9A-Za-z]+$", replacement = "", x)
    
  }
  )



# 3. Merge each VCF entry with corresponding SV name in reads list --------
merged_ids <- merge(x = vcf, y = reads, by.x = 'ID', by.y = 'SV_ID', all.x = TRUE)

# Reorder correctly
merged <- as.data.frame(merged_ids[match(vcf$ID, merged_ids$ID), ])


# 4. Export ---------------------------------------------------------------
write.table(x = merged[, c('CHROM', 'POS', 'ID', 'END', 'reads_reformat')], file = OUT,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
