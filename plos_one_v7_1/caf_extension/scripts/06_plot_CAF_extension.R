args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1],winslash="/",mustWork=TRUE)
extra_library <- Sys.getenv("FAP_R_LIBRARY", "")
if (nzchar(extra_library)) .libPaths(c(extra_library, .libPaths()))
library(ggplot2);library(patchwork)
o <- file.path(root,"figures");dir.create(o,showWarnings=FALSE)
rdir <- file.path(root,"results/adapted_v02")
s <- read.csv(file.path(rdir,"primary_single_cell_results.csv"))
n <- read.csv(file.path(rdir,"matched_null_summary.csv"))
b <- read.csv(file.path(rdir,"TCGA_CAF_associations.csv"))
keys <- c("Telocyte_CAF","Troph_MMPminus","Troph_MMPplus")
labels <- c("Telocyte-like","Trophocyte-like MMP-","Trophocyte-like MMP+")
for (df in c("s","n","b")) {
 d <- get(df);d$program <- factor(d$program,levels=keys,labels=labels);assign(df,d)
}
pal <- c(GSE132465="#26759E",GSE166555="#C77C3B")
theme_set(theme_classic(base_size=8,base_family="Arial")+
 theme(axis.line=element_line(linewidth=.3),axis.ticks=element_line(linewidth=.3),
 axis.text=element_text(colour="#303030",size=7),axis.title=element_text(size=8),
 plot.title=element_text(face="bold",size=9,margin=margin(b=5)),
 plot.subtitle=element_text(size=7.3,margin=margin(b=5)),
 strip.background=element_blank(),strip.text=element_text(face="bold",size=8),
 legend.title=element_blank(),legend.text=element_text(size=7),
 plot.margin=margin(5,6,5,5),plot.tag=element_text(face="bold",size=11)))
s$q_label <- ifelse(s$q<.001,"q < 0.001",sprintf("q = %.3f",s$q))
s$cohort_label <- factor(s$cohort,levels=names(pal),labels=c("GSE132465 | 15 patients","GSE166555 | 10 patients"))
a <- ggplot(s,aes(r,program,colour=cohort))+
 geom_vline(xintercept=0,colour="grey65",linetype=2,linewidth=.3)+
 geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=.13,linewidth=.55)+
 geom_point(size=2)+geom_text(aes(x=.71,label=q_label),hjust=1,size=2.45,show.legend=FALSE)+
 facet_grid(.~cohort_label)+scale_colour_manual(values=pal)+
 scale_x_continuous(limits=c(-.23,.73),breaks=c(-.2,0,.2,.4,.6))+
 labs(title="Within-patient CAF program-SenMayo associations",
 subtitle="Source-detectability-adapted programs; shared genes removed; depth and MKI67 adjusted",
 x="Mean within-patient partial correlation (bootstrap 95% CI)",y=NULL)+
 theme(legend.position="none")
n$row <- as.numeric(n$program)*2+ifelse(n$cohort=="GSE132465",.24,-.24)
n$cohort <- factor(n$cohort,levels=names(pal))
bp <- ggplot(n,aes(colour=cohort))+
 geom_vline(xintercept=0,colour="grey75",linetype=2,linewidth=.3)+
 geom_segment(aes(x=null_low,xend=null_high,y=row,yend=row),linewidth=1.8,alpha=.3)+
 geom_point(aes(x=null_median,y=row),shape=1,size=2,stroke=.6)+
 geom_point(aes(x=observed_r,y=row),shape=16,size=2)+
 scale_colour_manual(values=pal)+scale_y_continuous(breaks=c(2,4,6),labels=labels)+
 scale_x_continuous(limits=c(-.21,.53),breaks=c(-.2,0,.2,.4))+
 labs(title="Comparison with matched random programs",
 subtitle="Filled: observed | Open: null median\nBand: central 95% of 1,000 random programs",
 x="Adjusted within-patient correlation",y=NULL)+
 theme(legend.position="bottom",legend.key.width=unit(8,"mm"))
b <- b[b$target=="SenMayo",]
b$model <- factor(b$model,levels=c("marginal","fib5_adjusted"),labels=c("Marginal","fib5-adjusted"))
cp <- ggplot(b,aes(rho,program,colour=model,shape=model))+
 geom_vline(xintercept=0,colour="grey70",linetype=2,linewidth=.3)+
 geom_errorbar(aes(xmin=ci_low,xmax=ci_high),orientation="y",width=.12,linewidth=.5,position=position_dodge(width=.5))+
 geom_point(size=2,position=position_dodge(width=.5))+
 scale_colour_manual(values=c("Marginal"="#737373","fib5-adjusted"="#26759E"))+
 scale_shape_manual(values=c(1,16))+
 scale_x_continuous(limits=c(-.12,.95),breaks=c(0,.3,.6,.9))+
 labs(title="Bulk association and composition adjustment",subtitle="TCGA-COAD/READ | 380 primary-tumour patients",
 x="CAF program-SenMayo Spearman correlation",y=NULL)+
 theme(legend.position="bottom",legend.key.width=unit(8,"mm"))
# Three full-width rows keep long biological labels and uncertainty bars readable.
fig <- a / bp / cp + plot_layout(heights=c(1.1,1,1))+
 plot_annotation(tag_levels="a")
w <- 183/25.4;h <- 210/25.4
svglite::svglite(file.path(o,"CAF_extension_exploratory_v01.svg"),width=w,height=h);print(fig);dev.off()
grDevices::cairo_pdf(file.path(o,"CAF_extension_exploratory_v01.pdf"),width=w,height=h,family="Arial");print(fig);dev.off()
ragg::agg_png(file.path(o,"CAF_extension_exploratory_v01_preview.png"),width=w,height=h,units="in",res=200,background="white");print(fig);dev.off()
ragg::agg_tiff(file.path(o,"CAF_extension_exploratory_v01_600dpi.tiff"),width=w,height=h,units="in",res=600,compression="lzw",background="white");print(fig);dev.off()
write.csv(s,file.path(o,"panel_a_source_data.csv"),row.names=FALSE)
write.csv(n,file.path(o,"panel_b_source_data.csv"),row.names=FALSE)
write.csv(b,file.path(o,"panel_c_source_data.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,"qa/R_figure_sessionInfo.txt"))
cat("FIGURE EXPORTS COMPLETE\n")
