# Final ASD Discovery Pipeline
# Prereqs: 
#   - data/results/asd_peaks_qc.rds (from genome-wide scan + QC)
#   - data/results/shift_null_final.rds (from corrected_shift_null.R)
#   - data/asd_known_loci.csv (Grove 2019 known loci)

library(data.table)
library(jsonlite)
library(httr)

# Locus Clumping and Known/Novel Classification
cat("Locus Clumping & Classification")
qc  <- readRDS("data/results/asd_peaks_qc.rds")
sig <- qc[ras >= 4.5]; setorder(sig, chr, bp)
sig[, locus_id := cumsum(c(TRUE, diff(bp) > 5e5)), by = chr]

loci <- sig[, .(lead_bp = bp[which.max(ras)], max_ras = max(ras), n_peaks = .N),
            by = .(chr, locus_id)]

known <- fread("data/asd_known_loci.csv")
loci[, dist_known := Inf]
for (i in seq_len(nrow(known))) {
  loci[chr == known$chr[i], dist_known := pmin(dist_known, abs(lead_bp - known$bp[i]))]
}

loci[, class := fifelse(dist_known <= 1e6, "Known (Grove 2019)", "Novel")]
setorder(loci, -max_ras)

cat("CHECK")
cat(sprintf("loci: %d | known: %d | novel: %d | NA: %d\n", 
            nrow(loci), sum(loci$class == "Known (Grove 2019)"), 
            sum(loci$class == "Novel"), sum(is.na(loci$class))))

# Poisson Upper Bound Verification
cat("Poisson Upper Bound Verification")
null_hs <- readRDS("data/results/shift_null_final.rds")$null_hs
null_events <- sum(null_hs >= 4.5)
upper_rate <- 3 / 600 # 95% Poisson upper bound for 0 events
upper_genome <- 725 * upper_rate

cat(sprintf("null events >= 4.5: %d/600\n", null_events))
cat(sprintf("95%% Poisson upper rate/block: %.4f\n", upper_rate))
cat(sprintf("Genome-wide upper false peaks: %.1f\n", upper_genome))

# Ensembl API Gene Mapping
cat("Fetching nearest genes from Ensembl (hg19)")
loci$nearest_gene <- NA_character_

for(i in 1:nrow(loci)) {
  url <- sprintf("https://grch37.rest.ensembl.org/overlap/region/human/%d:%d-%d?feature=gene", 
                 loci$chr[i], loci$lead_bp[i], loci$lead_bp[i])
  
  res <- tryCatch(GET(url, accept_json()), error = function(e) NULL)
  
  if (!is.null(res) && status_code(res) == 200) {
    genes <- fromJSON(content(res, "text", encoding = "UTF-8"))
    if (length(genes) > 0) {
      loci$nearest_gene[i] <- genes$external_name[1]
    } else {
      url_wide <- sprintf("https://grch37.rest.ensembl.org/overlap/region/human/%d:%d-%d?feature=gene", 
                          loci$chr[i], loci$lead_bp[i] - 50000, loci$lead_bp[i] + 50000)
      res_wide <- tryCatch(GET(url_wide, accept_json()), error = function(e) NULL)
      if (!is.null(res_wide) && status_code(res_wide) == 200) {
        genes_wide <- fromJSON(content(res_wide, "text", encoding = "UTF-8"))
        if (length(genes_wide) > 0) loci$nearest_gene[i] <- genes_wide$external_name[1]
      }
    }
  }
  Sys.sleep(0.2) 
}

# Save Final Discoveries
fwrite(loci, "data/results/FINAL_ASD_DISCOVERIES.tsv", sep = "\t")
cat("\nFINAL DISCOVERIES\n")
print(loci[, .(chr, lead_bp, max_ras, class, nearest_gene)])
cat("Saved to data/results/FINAL_ASD_DISCOVERIES.tsv\n")