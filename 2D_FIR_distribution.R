#!/usr/bin/env Rscript


# 2D FIR distribution

# global variables
title             <- "RXLRs"
all_gene_file     <- "example_all_genes.gff3"
subset_gene_file  <- "example_gene_subset.gff3"
genome_sizes_file <- "example_genome.dat"
expression_file   <- "example_expression.dat"  # containing percentile data
threshold         <- 10000   # maximum size of intergenic distance for plotting
iterations        <- 1e4     # number of iterations in Monte Carlo simulations


# packages
require(hexbin)
require(RColorBrewer)
require(GISTools)
require(grid)


# functions
# calculate probabilities
calc_prop <- function(vector, observed_value){
	
	n_greater <- length(which(vector >= observed_value))
	n_less    <- length(which(vector <= observed_value))
	
	message(paste('Greater than', observed_value, ': p =', n_greater / length(vector)))
	message(paste('Less than', observed_value, ': p =', n_less / length(vector)))
}

# monte carlo simulator
monte_carlo <- function(dataframe, sample.size, iterations, observed.median){
	
	# make distribution vector
	medians <- vector(length=iterations)
	
	for (i in 1:iterations) {
		
		sim <- dataframe[sample(nrow(dataframe), sample.size), ]
		medians[i] <- median(c(sim$prime5, sim$prime3), na.rm=T)
	}
	
	# calculate how many simulated datasets are better than the observed number
	calc_prop(medians, observed.median)
}


# extract genes only from the GFF3 file (all genes)
command <- paste0("grep \"gene\" ", all_gene_file, " > all_genes.gff3")
cat(command,"\n")
try(system(command))

# get 5' and 3' closest genes (all genes)
command <- paste0("bedtools closest -id -io -D \"ref\" -g ", genome_sizes_file, " -t first -a all_genes.gff3 -b all_genes.gff3 | cut -f 14 > all_genes_5prime.dat")
cat(command,"\n")
try(system(command))

command <- paste0("bedtools closest -iu -io -D \"ref\" -g ", genome_sizes_file, " -t first -a all_genes.gff3 -b all_genes.gff3 | cut -f 13 > all_genes_3prime.dat")
cat(command,"\n")
try(system(command))

command <- "paste all_genes.gff3 all_genes_5prime.dat all_genes_3prime.dat > all_genes.dat"
cat(command,"\n")
try(system(command))


# extract genes only from the GFF3 file (subset genes)
command <- paste0("grep \"gene\" ", subset_gene_file, " > subset_genes.gff3")
cat(command,"\n")
try(system(command))

# get 5' and 3' closest genes (subset genes)
command <- paste0("bedtools closest -id -io -D \"ref\" -g ", genome_sizes_file, " -t first -a subset_genes.gff3 -b all_genes.gff3 | cut -f 14 > subset_genes_5prime.dat")
cat(command,"\n")
try(system(command))

command <- paste0("bedtools closest -iu -io -D \"ref\" -g ", genome_sizes_file, " -t first -a subset_genes.gff3 -b all_genes.gff3 | cut -f 13 > subset_genes_3prime.dat")
cat(command,"\n")
try(system(command))

command <- "paste subset_genes.gff3 subset_genes_5prime.dat subset_genes_3prime.dat > subset_genes.dat"
cat(command,"\n")
try(system(command))


# read data file (all genes)
d <- read.table("all_genes.dat", header=F)
d[d == -1] <- NA

# calculate distances to nearest genes (5' and 3') (all genes)
d$prime5 <- d[,4] - d[,10]
d$prime3 <- d[,11] - d[,5]


# read data file (subset genes)
s <- read.table("subset_genes.dat", header=F)
s[s == -1] <- NA

# calculate distances to nearest any gene (5' and 3') (subset genes)
s$prime5 <- s[,4] - s[,10]
s$prime3 <- s[,11] - s[,5]

# calculate median (subset genes)
d.med <- median(c(d$prime5, d$prime3), na.rm=T)
d.med

s.med <- median(c(s$prime5, s$prime3), na.rm=T)
s.med

# delete outliers to define plotting boundaries (all genes)
threshold
d$prime5[d$prime5 >  threshold] <- NA
d$prime5[d$prime5 < -threshold] <- NA
d$prime3[d$prime3 >  threshold] <- NA
d$prime3[d$prime3 < -threshold] <- NA


# read expression data
e <- read.table(expression_file, header=T)

# clean up gene names in dataframe from subset gene GFF3
s[,9] <- gsub(";+[A-Za-z0-9=]*$", "", s[,9])
s[,9] <- gsub("ID=", "", s[,9])

# make new percentile column
s$perc <- rep(0, length(s[,1]))

# add percentile numbers to subset gene dataframe
for( i in 1:length(s[,1])){
	
	s[i,]$perc <- e[which(e$gene == s[i,9]),]$percentile
}

# set color range
colors <- rep(0, 100)
colfunc <- colorRampPalette(c("#FFF5F0", "#67000D")) # Reds
#colfunc <- colorRampPalette(c("#FCFBFD", "#3F007D")) # Purples
#colfunc <- colorRampPalette(c("#FFFFCC", "#800026")) # YlOrRd
colors[0:100] <- colfunc(100)

# make new color column (with zero as white)
s$color <- rep(0, length(s[,1]))
for( i in 1:length(s[,1])){
	
	if(s[i,]$perc == 0){
		s[i,]$color <- "#FFFFFF"
	} else {
		s[i,]$color <- colors[s[i,]$perc]
	}
}

# color legend code
# plot(1:100, rep(1,100), pch=15, cex=3, col=colors)

# hexbin plot settings
bin <- hexbin(d$prime3/1000, d$prime5/1000, xbins=40)
my_colors <- colorRampPalette(rev(add.alpha(brewer.pal(11,'Spectral'),0.5)), alpha=T)

# hexbin plot
p <- plot(bin, main=title, colramp=my_colors, legend=F, xlab="3ʹ Flanking Intergenic Distance (Kb)", ylab="5ʹ Flanking Intergenic Distance (Kb)") 
pushHexport(p$plot.vp)
grid.points(s$prime3/1000, s$prime5/1000, pch=21, size=unit(6, "pt"), gp=gpar(fill=s$color))
upViewport()


# can change this to a log scale, but it can look quite misleading
# https://stackoverflow.com/questions/20651920/how-do-i-change-how-bins-are-assigned-in-hexbin-plot
# for log scales, axis tick marks may need to be adjusted manually

# hexbin plot settings
bin <- hexbin(log10(d$prime3/1000), log10(d$prime5/1000), xbins=40)
my_colors <- colorRampPalette(rev(add.alpha(brewer.pal(11,'Spectral'),0.5)), alpha=T)

# hexbin plot
p <- plot(bin, main=title, colramp=my_colors, legend=F, xlab="3ʹ Flanking Intergenic Distance (Kb)", ylab="5ʹ Flanking Intergenic Distance (Kb)") 
pushHexport(p$plot.vp)
grid.points(log10(s$prime3/1000), log10(s$prime5/1000), pch=21, size=unit(6, "pt"), gp=gpar(fill=s$color))
upViewport()


# monte carlo simulation (bootstrap without replacement)
# distance to the nearest gene (any), *not* the distance to the nearest subset gene
# otherwise subsetting obviously leads to larger distances

# run monte carlo
monte_carlo(d, sample.size=length(s$prime5), iterations, observed.median=s.med)

