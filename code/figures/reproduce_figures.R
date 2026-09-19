# Reproduce the five analysis figures from the frozen source tables.
# Required packages: tidyverse, patchwork, jsonlite, grid, gridExtra
# Run from the project root or set DATA and OUT below.

required <- c("tidyverse", "patchwork", "jsonlite", "gridExtra")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Installing missing R packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
suppressPackageStartupMessages({
  library(tidyverse); library(patchwork); library(jsonlite); library(gridExtra)
})

DATA <- "outputs/manuscript_figures_only_2026-09-06/source_data"
OUT  <- "outputs/manuscript_figures_only_2026-09-06/R_exports"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
blue <- "#286FA6"; dark <- "#173D5B"; gold <- "#C37A1F"
theme_pub <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(face = "bold", colour = dark),
        axis.title = element_text(colour = "#162D3D"),
        plot.title = element_text(face = "bold", colour = dark),
        legend.position = "none")

rank <- read_csv(file.path(DATA, "sequence_rankings.csv"), show_col_types = FALSE)
ev <- read_csv(file.path(DATA, "candidate_evidence.csv"), show_col_types = FALSE)
tcr <- read_tsv(file.path(DATA, "tcr_metrics.tsv"), show_col_types = FALSE)
gate <- fromJSON(file.path(DATA, "control_gate.json"), simplifyVector = FALSE)

# Figure 1: the workflow is intentionally schematic; counts come from the frozen tables.
workflow <- ggplot() + coord_cartesian(xlim=c(0,10), ylim=c(0,7), expand=FALSE) + theme_void() +
  annotate("rect", xmin=1.9, xmax=8.1, ymin=5.95, ymax=6.9, fill="#E6F0F7", colour=blue) +
  annotate("text", x=2.1, y=6.55, label="Peptide library", hjust=0, fontface="bold", colour=blue, size=4) +
  annotate("text", x=2.1, y=6.2, label="IEDB records + canonical human protein sequences", hjust=0, size=3.2) +
  annotate("rect", xmin=.1, xmax=4.45, ymin=3.95, ymax=5.25, fill="#E6F0F7", colour=blue) +
  annotate("text", x=.3, y=4.85, label="6,400 sequence comparisons", hjust=0, fontface="bold", colour=blue, size=4) +
  annotate("text", x=.3, y=4.35, label="1,600 per allele\nDRB1*03:01 · *08:01 · *13:03 · *15:01", hjust=0, size=3.1) +
  annotate("rect", xmin=5.55, xmax=9.9, ymin=3.95, ymax=5.25, fill="#E6F0F7", colour=blue) +
  annotate("text", x=5.75, y=4.85, label="Structural analysis", hjust=0, fontface="bold", colour=blue, size=4) +
  annotate("text", x=5.75, y=4.35, label="Predicted peptide–HLA complexes\nSeparate V3 structural shortlist", hjust=0, size=3.1) +
  annotate("rect", xmin=.1, xmax=4.45, ymin=1.95, ymax=3.2, fill="#E6F0F7", colour=blue) +
  annotate("text", x=.3, y=2.8, label="8 candidates reviewed", hjust=0, fontface="bold", colour=blue, size=4) +
  annotate("text", x=.3, y=2.3, label="Binding · registers · provenance\nAdditional evidence and gaps", hjust=0, size=3.1) +
  annotate("rect", xmin=5.55, xmax=9.9, ymin=1.95, ymax=3.2, fill="#E6F0F7", colour=blue) +
  annotate("text", x=5.75, y=2.8, label="Control benchmarks", hjust=0, fontface="bold", colour=blue, size=4) +
  annotate("text", x=5.75, y=2.3, label="Surface electrostatics\nTCR docking", hjust=0, size=3.1) +
  annotate("rect", xmin=.1, xmax=2.25, ymin=.03, ymax=1.13, fill="#E6F0F7", colour=blue) +
  annotate("text", x=.25, y=.75, label="2 medium priority\nBinding / register studies", hjust=0, size=3.2, fontface="bold", colour=blue) +
  annotate("rect", xmin=2.55, xmax=4.45, ymin=.03, ymax=1.13, fill="#FAEBD6", colour=gold) +
  annotate("text", x=2.7, y=.75, label="6 on hold\nEvidence gaps", hjust=0, size=3.2, fontface="bold", colour=gold) +
  annotate("rect", xmin=5.55, xmax=9.9, ymin=.03, ymax=1.13, fill="#FAEBD6", colour=gold) +
  annotate("text", x=5.75, y=.75, label="Control limitations\nNo additional candidate support", hjust=0, size=3.2, fontface="bold", colour=gold)
ggsave(file.path(OUT,"Figure_1_workflow.pdf"), workflow, width=10, height=7.3)

# Figure 2: identical bins across alleles; markers are the two reviewed BALF5–TALDO1 records.
leads <- rank %>% filter(ebv_protein == "BALF5", self_protein == "TALDO1")
lead_scores <- leads %>% group_by(allele) %>% slice_min(hla_rank, n=1) %>% ungroup()
p2 <- rank %>% mutate(allele = factor(allele, levels=c("HLA-DRB1*03:01","HLA-DRB1*08:01","HLA-DRB1*13:03","HLA-DRB1*15:01"))) %>%
  ggplot(aes(primary_score)) + geom_histogram(bins=30, fill=blue, colour="white") +
  facet_wrap(~allele, ncol=2, scales="free_y") + theme_pub +
  labs(x="BLOSUM62-based similarity score", y="Pair comparisons") +
  geom_vline(data=lead_scores, aes(xintercept=primary_score), colour=gold, linetype="dashed")
ggsave(file.path(OUT,"Figure_2_sequence_rankings.pdf"), p2, width=10, height=8)

# Figure 3: evidence matrix for the eight rows.
long <- ev %>% mutate(label=paste0(ebv_protein,"–",self_protein,"\n",sub("HLA-DRB1\\*","*",allele))) %>%
  transmute(label, `EBV binding`=if_else(ebv_binding_consensus,"Yes","No"), `Self binding`=if_else(self_binding_consensus,"Yes","No"),
            `EBV register`=if_else(ebv_register_consensus_matches_declared,"Yes","No"), `Self register`=if_else(self_register_consensus_matches_declared,"Yes","No"),
            `IEDB arms`=as.character(iedb_exact_hla_positive_arm_count), `Ligand arms`=as.character(immunopeptidome_exact_hla_compatible_arm_count),
            `Stage 1`=case_when(stage1_status=="stage1_medium_priority"~"Medium priority", TRUE~"Hold")) %>%
  pivot_longer(-label, names_to="metric", values_to="value")
p3 <- ggplot(long, aes(metric, fct_rev(label), fill=value)) + geom_tile(colour="white") + geom_text(aes(label=value), size=2.5) +
  scale_fill_manual(values=c("Yes"=blue,"No"="#EDF0F2","0"="#EDF0F2","1"=blue,"Hold"="#FAEBD6","Medium priority"="#E6F0F7"), na.value="#EDF0F2") +
  theme_pub + theme(axis.text.x=element_text(angle=35,hjust=1), axis.title=element_blank())
ggsave(file.path(OUT,"Figure_3_candidate_evidence.pdf"), p3, width=12, height=7.5)

# Figure S1: reconstruct status from the archived gate, rather than the legacy table.
fail <- map_dfr(gate$failures, as_tibble) %>% mutate(status="FAIL")
missing <- map_dfr(gate$missing, as_tibble) %>% mutate(status="NE")
gate_rows <- bind_rows(fail, missing) %>% select(endpoint, layer, pair_id, rank, status)
pass <- tibble(endpoint=c("peptide","composite","shape"), layer="af_314159", pair_id="PASS_CONTROL", rank=1:3, status="PASS")
s1 <- bind_rows(gate_rows, pass) %>% mutate(status=factor(status, levels=c("PASS","FAIL","NE")))
pS1 <- ggplot(s1, aes(endpoint, fct_rev(pair_id), fill=status)) + geom_tile(colour="white") + geom_text(aes(label=if_else(status=="NE","Not evaluable",paste0(if_else(status=="PASS","Rank ","Rank "),rank," · ",status))), size=2.4) +
  scale_fill_manual(values=c(PASS=blue,FAIL="#FAEBD6",NE="#EDF0F2")) + theme_pub + theme(axis.title=element_blank(), axis.text.y=element_text(size=6))
ggsave(file.path(OUT,"Figure_S1_electrostatics_controls.pdf"), pS1, width=10, height=7)

# Figure S2: archived calibration RMSDs.
tlong <- tcr %>% select(reference, pMHC_CA_RMSD_A, TCR_placement_CA_RMSD_A) %>% pivot_longer(-reference, names_to="metric", values_to="rmsd") %>% mutate(metric=recode(metric,pMHC_CA_RMSD_A="pMHC geometry",TCR_placement_CA_RMSD_A="TCR placement"))
pS2 <- ggplot(tlong, aes(rmsd, fct_rev(reference), colour=metric)) + geom_segment(aes(x=0,xend=rmsd,yend=reference), linewidth=1.2) + geom_point(size=3) + geom_text(aes(label=sprintf("%.2f",rmsd)), nudge_y=.15, show.legend=FALSE) + facet_wrap(~metric, nrow=1) + scale_colour_manual(values=c("pMHC geometry"=blue,"TCR placement"=gold)) + theme_pub + labs(x="Cα RMSD (Å)", y=NULL)
ggsave(file.path(OUT,"Figure_S2_TCR_docking_controls.pdf"), pS2, width=10, height=6.6)

message("Wrote five PDFs to ", normalizePath(OUT))
