# Load libraries
library('GenomicFeatures')
library('tximport')
library('tidyverse')
library('txdbmaker')

# Build a TxDb from triticale gff3
txdb <- makeTxDbFromGFF('./data/Triticale_Svevov1_RyeLo72018v1p1p1.Ensembl59.gff3', format = 'gff3')

# Extract transcript to gene mapping
k <- keys(txdb, keytype = 'TXNAME')
tx2gene <- select(txdb, keys = k, keytype = 'TXNAME', columns = 'GENEID')

tx2gene %>% View()

write.table(file = './result/annotation.txt', tx2gene, sep = ' ', row.names = FALSE)











