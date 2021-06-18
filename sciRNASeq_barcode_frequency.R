library(readr)
library(dplyr)

calculate_barcode_frequency = function(barcode_file, reaper_whitelist) {
  RT_barcodes <- read_csv(barcode_file, col_names = FALSE)
  RT_grouped = RT_barcodes %>% group_by(X1) %>% summarise(count = n())
  
  reaper_barcodes <- read_delim(reaper_whitelist, "\t", escape_double = FALSE, trim_ws = TRUE)
  reaper_barcodes$present = T
  reaper_barcodes = reaper_barcodes %>% select(barcode, present)

  barcode_seq = merge(RT_grouped, reaper_barcodes, by.x = "X1", by.y = "barcode", all.x = T, all.y = F)
  barcode_seq = barcode_seq[order(barcode_seq$count, decreasing = T),]
  barcode_seq$present[is.na(barcode_seq$present)] <- F
  rownames(barcode_seq) = NULL
  barcode_seq
}


barcode_seq = calculate_barcode_frequency("~/sci_rnaseq_round_02/reaper_out/lint_barcodes.txt", "~/sci_rnaseq_round_02/reaper_barcodes_4.txt")


write.table(barcode_seq, "lint_barcode_seq.txt", sep="\t", quote = F, row.names = T)
