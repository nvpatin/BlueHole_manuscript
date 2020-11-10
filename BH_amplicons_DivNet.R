#load
library(phyloseq)
library(breakaway)
library(DivNet)
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(data.table)
library(textshape)
library(magrittr) # to use %>% as a pipe
library(vegan)
library(doParallel)
library(dplyr)

# import ASV table (eg "table_IDs_sterivex") and ID-to-taxonomy table

microbe <- column_to_rownames(table_hole_2019_water_phyloseq_filt, 'ID')
tax_table <- column_to_rownames(taxonomy_phyloseq, 'ID')
metadata <- column_to_rownames(BH_hole_water_metadata_phyloseq, 'Sample')

# convert all to matrices
microbe <- as.matrix(microbe)
tax_table <- as.matrix(tax_table)
metadata <- as.data.frame(metadata)

# combine into phyloseq object
ASV = otu_table(microbe, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
MET = sample_data(metadata)
physeq = phyloseq(ASV, TAX, MET)

### ALPHA DIVERSITY ###

# Aggregate taxa at the order or family level
water <- physeq %>% tax_glom("Order")
water <- physeq %>% tax_glom("Family")

ba <- breakaway(water)

plot(ba, water, color = "DepthGroup")

# Run DivNet at various taxonomic levels: phylum, genus, and ASV
# For Blue Hole manuscript I ran it at the ASV level using the ASV ID provided
divnet_phylum <- divnet(tax_glom(physeq, taxrank="Phylum"), X="Location", ncores=6)
divnet_genus <- divnet(tax_glom(physeq, taxrank="Genus"), X="Depth", ncores=4)
divnet_asv <- divnet(physeq, base="580c721cd115033cebadc0d909e4cd83", ncores=4)

# Alpha diversity
plot(divnet_asv$shannon, 
     physeq, 
     col = "DepthGroup")

divnet_depth <- physeq %>%
  divnet(X = "DepthGroup", ncores = 2)

plot(divnet_depth$shannon, 
     physeq, 
     col = "DepthGroup")

plot(divnet_depth$simpson, 
     physeq, 
     col = "DepthGroup")

testDiversity(divnet_depth, "shannon")

# Make table of alpha diversity metrics
estimates_shan <- divnet_asv$shannon %>% summary %$% estimate
estimates_simp <- divnet_asv$simpson %>% summary %$% estimate
alphadiv_table <- subset(dftr_meta, select = V1)

alphadiv_table[["Shannon"]] <- estimates_shan
alphadiv_table[["Simpson"]] <- estimates_simp

alphadiv_table <- subset(alphadiv_table, select=c("Shannon","Simpson"))
alphadiv_table[["DepthGroup"]] <- BH_hole_water_metadata_phyloseq_v2[["DepthGroup"]]
alpha <- tibble::rownames_to_column(alphadiv_table, "Sample")

amelt <- melt(alpha, id=c("Sample", "DepthGroup"))

p1 <- ggplot(amelt, aes(variable, value, colour=DepthGroup))
p1 + geom_boxplot() + theme_classic() + facet_wrap(~variable, scales="free")

################ BETA DIVERSITY #################

# Pull out Bray-Curtis distances as a square distance matrix
bc <- divnet_asv$'bray-curtis'

# Base R PCA function from distance matrix
PCA <- prcomp(x=bc)
PCAi <- data.frame(PCA$x)

# Base R plotting for quick check
autoplot(PCA)
p <- ggplot(PCAi, aes(x=PC1, y=PC2)) + geom_point(size=3)

# Pretty plot in ggplot2
# copy data frame and add any metadata columns
bc_meta <- bc
# make sure the metadata is in the correct sample order
metadata <- metadata[ order(row.names(metadata)), ]

bc_meta$Depth = metadata$DepthGroup
bc_meta$Month = as.factor(metadata$Month)

PCAi <- data.frame(PCA$x, group=bc_meta$Depth, shape=bc_meta$Month)

# Make a custom color pallette
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

p <- ggplot(PCAi, aes(x=PC1, y=PC2)) + 
    geom_point(alpha=0.8, size=4, aes(fill=group, shape=shape)) +
    scale_shape_manual(values=c(21, 24)) + 
    scale_fill_manual(values=cbp1) + labs(fill="Depth Group", shape="Month") +
    guides(fill=guide_legend(override.aes=list(shape=21)))
    
p + theme_classic() + labs(x = "PC1", y = "PC2") 