library(ggplot2)
out <- '/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/manuscript_figures_v2'
blue <- '#2F6B9A'; gold <- '#B87918'; ink <- '#233142'; grey <- '#637381'; pale <- '#E9EFF3'
theme_paper <- function(){theme_classic(base_size=13, base_family='Helvetica')+theme(plot.title=element_blank(),plot.subtitle=element_blank(),legend.title=element_blank(),legend.position='none',axis.title=element_text(colour=ink, size=13),axis.text=element_text(colour=ink, size=11),axis.line=element_line(colour=ink),plot.margin=margin(12,36,12,12))}
# Figure 2: within-allele score distributions
ranked <- read.csv('processed/literature_grounded_hla2_rankings_v3_2026-08-27/v3_all_hla_ranked_pairs.csv',check.names=FALSE)
ranked$allele_short <- sub('HLA-','',ranked$allele)
ranked$highlight <- ifelse((ranked$allele=='HLA-DRB1*13:03' & ranked$ebv_protein=='BALF5' & ranked$self_protein=='TALDO1') | (ranked$allele=='HLA-DRB1*15:01' & ranked$ebv_protein=='BALF5' & ranked$self_protein=='TALDO1'),'BALF5-TALDO1','Other comparisons')
hits <- subset(ranked, highlight=='BALF5-TALDO1' & primary_rank %in% c(13,14))
p2 <- ggplot(ranked,aes(allele_short,primary_score))+
 geom_violin(fill=pale,colour='#93A6B5',linewidth=.35,trim=FALSE)+
 geom_jitter(data=subset(ranked, highlight=='Other comparisons'),width=.12,height=0,alpha=.12,size=.45,colour='#7890A0')+
 geom_point(data=hits,aes(colour=allele_short),size=3.8)+
 geom_label(data=hits,aes(label=paste0('rank ',primary_rank,'/1,600'),colour=allele_short),fill='white',label.size=0,size=4,nudge_y=.02,show.legend=FALSE)+
 scale_colour_manual(values=c('DRB1*03:01'='#7A8B99','DRB1*08:01'='#7A8B99','DRB1*13:03'=blue,'DRB1*15:01'=gold))+
 labs(x=NULL,y='TCR-facing BLOSUM62 similarity')+theme_paper()+theme(axis.text.x=element_text(angle=20,hjust=1))
ggsave(file.path(out,'figure_2_score_distributions_restyled.pdf'),p2,width=7.4,height=4.6,device='pdf',useDingbats=FALSE)
# Figure 4: 4-metric slope chart with true reverse log scale and collision-free labels
s <- data.frame(
  method_num = rep(1:4, 2),
  rank = c(13, 173, 2, 15,   # DRB1*13:03
           14, 3, 4, 33),     # DRB1*15:01
  context = rep(c('DRB1*13:03', 'DRB1*15:01'), each=4)
)

# Explicit non-colliding positions for rank labels
lbls <- data.frame(
  x = c(0.85, 1.15, 2.00, 2.00, 2.85, 3.15, 3.88, 3.88),
  y = c(12.0, 15.0, 280,  1.9,  1.7,  4.8,  12.0, 42.0),
  label = c('13', '14', '173', '3', '2', '4', '15', '33'),
  context = c('DRB1*13:03', 'DRB1*15:01', 'DRB1*13:03', 'DRB1*15:01',
              'DRB1*13:03', 'DRB1*15:01', 'DRB1*13:03', 'DRB1*15:01'),
  hjust = c(1, 0, 0.5, 0.5, 1, 0, 1, 1),
  vjust = c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
)

p4 <- ggplot(s, aes(x=method_num, y=rank, group=context, colour=context)) +
  geom_hline(yintercept=c(1, 3, 10, 30, 100, 300, 1000, 1600), colour='#EBF0F3', linewidth=0.35) +
  geom_hline(yintercept=16, linetype='dashed', colour='#8B9DAE', linewidth=0.55) +
  geom_hline(yintercept=80, linetype='dotted', colour='#8B9DAE', linewidth=0.65) +
  geom_line(linewidth=1.35) +
  geom_point(size=4.2) +
  geom_text(data=lbls, aes(x=x, y=y, label=label, colour=context, hjust=hjust, vjust=vjust),
            fontface='bold', size=4.2, inherit.aes=FALSE) +
  annotate('segment', x=0.75, xend=1.05, y=950, yend=950,
           linetype='dashed', colour='#8B9DAE', linewidth=0.65) +
  annotate('text', x=1.10, y=950, label='Top 1% threshold (rank <= 16)',
           colour='#506070', size=3.4, hjust=0) +
  annotate('segment', x=2.35, xend=2.65, y=950, yend=950,
           linetype='dotted', colour='#8B9DAE', linewidth=0.75) +
  annotate('text', x=2.70, y=950, label='Top 5% threshold (rank <= 80)',
           colour='#506070', size=3.4, hjust=0) +
  annotate('text', x=4.14, y=15, label='DRB1*13:03\nBALF5-TALDO1 (108-122)',
           hjust=0, vjust=0.5, colour=blue, size=3.9, fontface='bold', lineheight=0.95) +
  annotate('text', x=4.14, y=33, label='DRB1*15:01\nBALF5-TALDO1 (216-230)',
           hjust=0, vjust=0.5, colour=gold, size=3.9, fontface='bold', lineheight=0.95) +
  scale_colour_manual(values=c('DRB1*13:03'=blue, 'DRB1*15:01'=gold)) +
  scale_x_continuous(
    breaks=1:4,
    labels=c('TCR-facing\nBLOSUM62', 'Full-core\nBLOSUM62', 'TCR-facing\nidentity', 'Grantham\ndistance')
  ) +
  scale_y_continuous(
    trans=c('log10', 'reverse'),
    limits=c(1600, 1),
    breaks=c(1, 3, 10, 30, 100, 300, 1000, 1600),
    labels=c('1', '3', '10', '30', '100', '300', '1,000', '1,600')
  ) +
  labs(x=NULL, y='Within-allele rank (1 = highest)') +
  theme_paper() +
  theme(
    axis.text.x = element_text(size=11, colour=ink, face='bold'),
    plot.margin = margin(14, 155, 12, 14)
  ) +
  coord_cartesian(xlim=c(0.65, 4.05), ylim=c(1600, 1), clip='off')
ggsave(file.path(out,'figure_4_rank_sensitivity_restyled.pdf'),p4,width=7.8,height=4.8,device='pdf',useDingbats=FALSE)
ggsave(file.path(out,'figure_4_rank_sensitivity.pdf'),p4,width=7.8,height=4.8,device='pdf',useDingbats=FALSE)

# Figure 5: structural RMSD range with Å unit and explicit median callouts
d5 <- data.frame(
  label=factor(c('DRB1*13:03\nBALF5','DRB1*13:03\nTALDO1 (108-122)','DRB1*15:01\nBALF5','DRB1*15:01\nTALDO1 (216-230)'),
               levels=c('DRB1*13:03\nBALF5','DRB1*13:03\nTALDO1 (108-122)','DRB1*15:01\nBALF5','DRB1*15:01\nTALDO1 (216-230)')),
  low=c(.316,.518,.154,.311),
  median=c(.445,.900,1.712,.888),
  high=c(.788,20.176,4.274,1.515),
  allele=rep(c('DRB1*13:03','DRB1*15:01'),each=2),
  med_str=c('0.45 \u00c5', '0.90 \u00c5', '1.71 \u00c5', '0.89 \u00c5')
)
p5 <- ggplot(d5,aes(label,median,colour=allele))+
  geom_hline(yintercept=c(0.1, 0.3, 1, 3, 10, 30), colour='#EBF0F3', linewidth=0.35) +
  geom_linerange(aes(ymin=low,ymax=high),linewidth=1.5)+
  geom_point(shape=21,fill='white',size=5.2,stroke=1.8)+
  geom_text(aes(label=med_str), hjust=-0.35, vjust=0.45, size=4.1, fontface='bold', show.legend=FALSE) +
  scale_colour_manual(values=c('DRB1*13:03'=blue,'DRB1*15:01'=gold))+
  scale_y_log10(breaks=c(.1,.3,1,3,10,30),labels=c('0.1','0.3','1','3','10','30'),limits=c(0.1,30))+
  labs(x=NULL,y='Parent-versus-nested peptide RMSD (\u00c5)')+
  theme_paper()+
  theme(
    axis.text.x=element_text(size=11, colour=ink, face='bold'),
    plot.margin=margin(14, 38, 12, 14)
  )+
  coord_cartesian(clip='off')
ggsave(file.path(out,'figure_5_structural_input_sensitivity_restyled.pdf'),p5,width=7.4,height=4.6,device='pdf',useDingbats=FALSE)
ggsave(file.path(out,'figure_5_structural_input_sensitivity.pdf'),p5,width=7.4,height=4.6,device='pdf',useDingbats=FALSE)
# Figure 6: binding-percentile quadrant scatter
b <- read.csv('processed/high_yield_candidate_evidence_2026-08-28/predictor_register_comparison.csv',check.names=FALSE)
b$group <- ifelse(b$target_id %in% c('HY13_SEQ_02','HY15_SEQ_02'),'BALF5-TALDO1 arms','Other audit arms')
b$label <- ''
b$label[b$arm_id=='HY13_SEQ_02__ebv'] <- 'DR13 viral'
b$label[b$arm_id=='HY15_SEQ_02__ebv'] <- 'DR15 viral'
b$label[b$arm_id=='HY15_SEQ_02__self'] <- 'DR15 self'
b$label_x <- b$netmhciipan_el_percentile
b$label_y <- b$mixmhc2pred_context_percentile
b$label_x[b$arm_id=='HY13_SEQ_02__ebv'] <- .33
b$label_y[b$arm_id=='HY13_SEQ_02__ebv'] <- .16
b$label_x[b$arm_id=='HY15_SEQ_02__ebv'] <- 6.6
b$label_y[b$arm_id=='HY15_SEQ_02__ebv'] <- 3.7
b$label_x[b$arm_id=='HY15_SEQ_02__self'] <- 14.8
b$label_y[b$arm_id=='HY15_SEQ_02__self'] <- 1.45
p6 <- ggplot(b,aes(netmhciipan_el_percentile,mixmhc2pred_context_percentile))+
 annotate('rect',xmin=.003,xmax=5,ymin=.003,ymax=5,fill='#DDEDE7',alpha=.45)+
 annotate('rect',xmin=5,xmax=20,ymin=5,ymax=20,fill='#F7EEDB',alpha=.5)+
 geom_vline(xintercept=c(5,20),linetype=c('dashed','dotted'),colour=c(gold,grey),linewidth=.45)+geom_hline(yintercept=c(5,20),linetype=c('dashed','dotted'),colour=c(gold,grey),linewidth=.45)+
 geom_point(aes(colour=group),size=3.8,alpha=.9)+
 geom_text(data=subset(b,label!=''),aes(x=label_x,y=label_y,label=label,colour=group),hjust=0,size=3.8,show.legend=FALSE)+
 scale_colour_manual(values=c('BALF5-TALDO1 arms'=blue,'Other audit arms'='#93A6B5'))+
 scale_x_log10(limits=c(.003,100),breaks=c(.01,.1,1,5,20,100),labels=c('0.01','0.1','1','5','20','100'))+
 scale_y_log10(limits=c(.003,100),breaks=c(.01,.1,1,5,20,100),labels=c('0.01','0.1','1','5','20','100'))+
 labs(x='NetMHCIIpan-4.3 EL percentile rank',y='MixMHC2pred context percentile rank')+theme_paper()+theme(legend.position=c(.75,.12),legend.background=element_rect(fill='white',colour=NA),legend.text=element_text(size=10))
ggsave(file.path(out,'figure_6_initial_audit_binding_predictions_restyled.pdf'),p6,width=6.8,height=5.3,device='pdf',useDingbats=FALSE)
