# Get amino acid sequence coordinate ranges for classical HLA loci
# Uses the HLAtools package to build protein alignments from IPD-IMGT/HLA

library(HLAtools)

data(HLAgazetteer)

classical_loci <- HLAgazetteer$classical
cat("Classical HLA loci:", paste(classical_loci, collapse = ", "), "\n")
cat("IPD-IMGT/HLA Database version:", HLAgazetteer$version, "\n\n")

# Build protein alignments for all classical loci
alignments <- alignmentFull(loci = classical_loci, alignType = "prot")

# Metadata columns to exclude when extracting position numbers
meta_cols <- c("locus", "allele", "trimmed_allele", "allele_name")

# Extract position ranges for each locus
results <- data.frame(
  locus = character(),
  min_position = integer(),
  max_position = integer(),
  num_positions = integer(),
  stringsAsFactors = FALSE
)

for (loc in classical_loci) {
  if (!is.null(alignments$prot[[loc]])) {
    aln <- alignments$prot[[loc]]
    pos_cols <- setdiff(colnames(aln), meta_cols)
    positions <- as.integer(pos_cols)
    results <- rbind(results, data.frame(
      locus = loc,
      min_position = min(positions),
      max_position = max(positions),
      num_positions = length(positions),
      stringsAsFactors = FALSE
    ))
  } else {
    cat("WARNING: No protein alignment available for", loc, "\n")
  }
}

cat("Amino Acid Sequence Position Ranges for Classical HLA Loci\n")
cat("===========================================================\n")
print(results, row.names = FALSE)
