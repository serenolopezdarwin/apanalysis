library(ggseqlogo)

pdf("figs/Fig1D.logo.all.original.100K.pdf",width=20,height=8)
# par(mar=c(14, 14, 10, 10), mgp = c(5, 1, 0))
a = data <- read.table(file("stdin"), header=F, quote="")
ggseqlogo( as.character(a$V2) )
dev.off()

(mat = Biostrings::consensusMatrix( toupper(as.character(a$V2)), as.prob=T))

pdf("figs/Fig1D.PFM.all.original.100K.pdf",width=20,height=8)
mat=mat[rownames(mat) %in% c("A","T","C","G"),]
# par(mar=c(14, 14, 10, 10), mgp = c(5, 1, 0))
matplot(t(mat), type='l', lty=1, lwd=2)
legend("topright", rownames(mat), col=seq_len(4), fill=seq_len(4))
dev.off()
