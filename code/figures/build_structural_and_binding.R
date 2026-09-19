out_dir <- '/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/manuscript_figures_v2'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
blue <- '#2F6B9A'; gold <- '#B87918'; ink <- '#233142'; gray <- '#637381'; light_gray <- '#D9E0E5'; red <- '#AE4A45'
open_device <- function(stem, ext, w, h) {
  p <- file.path(out_dir, paste0(stem, '.', ext)); if (ext=='pdf') pdf(p,width=w,height=h,useDingbats=FALSE) else png(p,width=w,height=h,units='in',res=400)
}
for (ext in c('pdf')) {
  open_device('figure_5_structural_input_sensitivity',ext,7.4,4.5)
  par(mar=c(5.1,4.8,.6,.4),xaxs='i',yaxs='i')
  d <- data.frame(
    label=c('DRB1*13:03\nBALF5', 'DRB1*13:03\nTALDO1', 'DRB1*15:01\nBALF5', 'DRB1*15:01\nTALDO1'),
    low=c(.316,.518,.154,.311), median=c(.445,.900,1.712,.888), high=c(.788,20.176,4.274,1.515),
    col=c(blue,blue,gold,gold))
  x<-1:4
  plot(x,d$median,log='y',ylim=c(.1,30),xlim=c(.55,4.45),type='n',xaxt='n',yaxt='n',xlab='',ylab='Parent-versus-nested peptide RMSD (A)')
  axis(1,at=x,labels=d$label,tick=FALSE,cex.axis=.86); axis(2,at=c(.1,.3,1,3,10,30),labels=c('0.1','0.3','1','3','10','30'),las=1)
  abline(h=c(.1,.3,1,3,10,30),col=light_gray,lwd=.7)
  for(i in x){segments(i,d$low[i],i,d$high[i],lwd=3,col=d$col[i]); points(i,d$median[i],pch=21,bg='white',col=d$col[i],cex=1.8,lwd=2)}
  dev.off()
}
for (ext in c('pdf')) {
  open_device('figure_6_initial_audit_binding_predictions',ext,6.6,5.3)
  d <- read.csv('processed/high_yield_candidate_evidence_2026-08-28/predictor_register_comparison.csv',check.names=FALSE)
  x <- d$netmhciipan_el_percentile; y <- d$mixmhc2pred_context_percentile
  special <- d$target_id %in% c('HY13_SEQ_02','HY15_SEQ_02')
  par(mar=c(4.7,5,1.0,.7),xaxs='i',yaxs='i')
  plot(x,y,log='xy',xlim=c(.003,100),ylim=c(.003,100),type='n',xlab='NetMHCIIpan-4.3 EL percentile rank',ylab='MixMHC2pred context percentile rank',xaxt='n',yaxt='n')
  ticks<-c(.01,.1,1,5,20,100); labs<-c('0.01','0.1','1','5','20','100')
  axis(1,at=ticks,labels=labs);axis(2,at=ticks,labels=labs,las=1)
  abline(v=c(5,20),h=c(5,20),col=c(gold,gray),lty=c(2,3),lwd=1)
  points(x[!special],y[!special],pch=21,bg='#AAB7C2',col='white',cex=1.35)
  points(x[special],y[special],pch=21,bg=blue,col='white',cex=1.55)
  idx13<-which(d$arm_id=='HY13_SEQ_02__ebv'); idx15v<-which(d$arm_id=='HY15_SEQ_02__ebv'); idx15s<-which(d$arm_id=='HY15_SEQ_02__self')
  text(x[idx13],y[idx13], 'DR13 viral',pos=4,offset=.45,cex=.72,col=blue)
  text(x[idx15v],y[idx15v], 'DR15 viral',pos=4,offset=.45,cex=.72,col=blue)
  text(12.5,3.15, 'DR15 self',cex=.72,col=blue)
  legend('bottomright',legend=c('other arms in the eight-pair audit','BALF5-TALDO1 arms'),pch=21,pt.bg=c('#AAB7C2',blue),pt.cex=1.3,bty='n',cex=.72)
  dev.off()
}
