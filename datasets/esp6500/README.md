1. Download CADD ESP6500 annotated file.
2. R

r <- read.table("ESP6500SI_inclAnno.tsv", header=T, skip=1, sep="\t", stringsAsFactors=F, comment.char="")
r.maf.5_10 <- r[which(r$ESP_AF >= 0.05 & r$ESP_AF <= 0.1),]
chrom.pos <- paste(r.maf.5_10$X.Chrom, r.maf.5_10$Pos, sep=":")
r.maf.5_10$cp <- chrom.pos
rsorted <- r.maf.5_10[sort.list(r.maf.5_10$cp),]
nondupIdx <- which(!duplicated(rsorted$cp))
nondups <- rsorted[nondupIdx,]
nonbenign.idx <- which(nondups$SIFTcat == "deleterious" | nondups$PolyPhenCat == "possibly_damaging" | nondups$PolyPhenCat == "probably_damaging")
nrow(nondups[-nonbenign.idx,])
esp.benign <- nondups[-nonbenign.idx,]
