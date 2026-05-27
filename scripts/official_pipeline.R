library(openPrimeR)
# Specify a FASTA file containing the templates:
fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                          "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
# Load the template sequences from 'fasta.file'
seq.df.simple <- read_templates(fasta.file)
seq.df.simple$Header[1]
hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
seq.df <- read_templates(fasta.file, hdr.structure, delim = "|", id.column = "GROUP")
# Show the loaded metadata for the first template
c(seq.df$Accession[1], seq.df$Group[1], seq.df$Species[1], seq.df$Function[1])
seq.df$ID[1]
seq.df$Sequence[1]
seq.df$Allowed_fw[1]
seq.df$Allowed_rev[1]

l.fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                            "Homo_sapiens_IGH_functional_leader.fasta", package = "openPrimeR")
template.df <- assign_binding_regions(seq.df, fw = l.fasta.file, rev = NULL)


list.files(system.file("extdata", "settings", package = "openPrimeR"), pattern = "*\\.xml")
settings.xml <- system.file("extdata", "settings", 
                            "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
settings <- read_settings(settings.xml)


# Define the FASTA primer file to load
primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                               "Ippolito2012.fasta", package = "openPrimeR")
# Load the primers
primer.df <- read_primers(primer.location, fw.id = "_fw")
# Allow off-target coverage
conOptions(settings)$allowed_other_binding_ratio <- c("max" = 1.0)
# Evaluate all constraints found in 'settings'
constraint.df <- check_constraints(primer.df, template.df, 
                                   settings, active.constraints = names(constraints(settings)))
constraint.df$primer_coverage
template.df <- update_template_cvg(template.df, constraint.df)
template.df$primer_coverage[1:5]
as.numeric(get_cvg_ratio(constraint.df, template.df))
cvg.stats <- get_cvg_stats(constraint.df, template.df, for.viewing = TRUE)
cvg.stats


# Optimal primer subsets
primer.subsets <- subset_primer_set(constraint.df, template.df)
plot_primer_subsets(primer.subsets, template.df)
plot_constraint_fulfillment(constraint.df, settings)