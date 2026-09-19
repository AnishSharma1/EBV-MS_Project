out_dir <- '/Users/anishsharma/Library/Mobile Documents/com~apple~CloudDocs/ebv_ms_publication/manuscript_figures_v2'
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

ink <- '#233142'; blue <- '#2F6B9A'; light_blue <- '#E6F0F7'; gold <- '#B87918'; light_gold <- '#F7EEDB'; gray <- '#637381'; red <- '#AE4A45'
open_device <- function(stem, ext, w, h) {
  p <- file.path(out_dir, paste0(stem, '.', ext))
  if (ext == 'png') png(p, width=w, height=h, units='in', res=400, type='cairo')
  if (ext == 'pdf') pdf(p, width=w, height=h, useDingbats=FALSE)
  if (ext == 'svg') svg(p, width=w, height=h)
}
close_device <- function() dev.off()
for (ext in c('png','pdf','svg')) {
  open_device('figure_1_branched_workflow', ext, 8.6, 4.8)
  par(mar=c(0,0,0,0), xaxs='i', yaxs='i', family='Helvetica')
  plot.new(); plot.window(xlim=c(0,8.6), ylim=c(0,3.9))
  box <- function(x,y,w,h,txt,fill,border=ink,cex=.82,bold=FALSE) {
    rect(x,y,x+w,y+h,col=fill,border=border,lwd=1.4)
    text(x+w/2,y+h/2,txt,cex=cex,col=ink,font=2)
  }
  arrow <- function(x1,y1,x2,y2,col) arrows(x1,y1,x2,y2,length=.09,lwd=1.6,col=col)
  box(.20,1.18,1.25,1.05,'6,400\nallele specific\ncomparisons',light_blue,cex=.86,bold=TRUE)
  box(1.70,2.20,1.45,.92,'8 selected for\nevidence audit',light_blue,blue)
  box(3.70,2.20,1.65,.92,'2 medium priority\n6 holds\n0 high priority',light_blue,blue)
  box(1.70,.78,1.45,.92,'49 pair\nexpanded cohort',light_gold,gold)
  box(3.70,.78,1.65,.92,'44 with register\npredictor agreement',light_gold,gold)
  box(5.75,.64,1.45,1.20,'2 pairs\nwith complete\ncomputational\ncriteria',light_gold,gold)
  box(7.30,.78,1.15,.92,'7 provisional\ntiers\n42 untiered',light_gold,gold,cex=.72)
  arrow(1.45,1.70,1.70,2.66,blue); arrow(1.45,1.70,1.70,1.24,gold)
  arrow(3.15,2.66,3.70,2.66,blue); arrow(3.15,1.24,3.70,1.24,gold); arrow(5.35,1.24,5.75,1.24,gold); arrow(7.20,1.24,7.30,1.24,gold)
  rect(2.22,3.54,2.32,3.64,col=light_blue,border=blue,lwd=1.1)
  text(2.40,3.59,'Blue = initial eight pair audit',adj=c(0,0.5),cex=.64,col=ink,font=2)
  rect(4.98,3.54,5.08,3.64,col=light_gold,border=gold,lwd=1.1)
  text(5.16,3.59,'Orange = later 49 pair expansion',adj=c(0,0.5),cex=.64,col=ink,font=2)
  close_device()
}

for (ext in c('png','pdf','svg')) {
  open_device('figure_4_rank_sensitivity', ext, 7.2, 4.2)
  par(mar=c(4.4,4.6,0.5,.5), xaxs='i', yaxs='i')
  ranks <- rbind(c(13,173,2,15),c(14,3,4,33))
  cols <- c(blue,gold); xx <- 1:4
  plot(xx, ranks[1,], log='y', ylim=c(1600,1), xlim=c(.8,4.2), type='n', xaxt='n', yaxt='n', xlab='', ylab='Within-allele rank (1 = highest)')
  axis(1,at=xx,labels=c('TCR-facing\nBLOSUM62','Full-core\nBLOSUM62','TCR-facing\nidentity','Grantham\ndistance'),tick=FALSE,cex.axis=.80)
  axis(2,at=c(1,3,10,30,100,300,1000,1600),labels=c('1','3','10','30','100','300','1,000','1,600'),las=1)
  abline(h=c(1,3,10,30,100,300,1000,1600),col='#D9E0E5',lwd=.7)
  for (i in 1:2) {
    lines(xx,ranks[i,],type='b',pch=16,lwd=2,col=cols[i])
    for (j in 1:4) text(xx[j] + ifelse(i == 1, -.06, .06), ranks[i,j], labels=ranks[i,j], pos=ifelse(ranks[i,j]<30,1,3), offset=.65, cex=.82, col=cols[i], font=2)
  }
  legend('bottomleft',legend=c('DRB1*13:03','DRB1*15:01'),col=cols,lty=1,pch=16,lwd=2,bty='n',cex=.78)
  close_device()
}
