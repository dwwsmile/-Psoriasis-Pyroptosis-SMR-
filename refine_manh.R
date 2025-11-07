
# R script to draw Manhattan plot for SMR analysis
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggrepel)

Manhattan_refine <- function(Background_data, highlight_names = "",highlight_size = 2,
                             pval_names = "p value",thresholds = 0.05,thresholds_line_color = "red",
                             points_color = "skyblue", highlight_color = "orange",
                             pointsize = 2,pointshape=18,label_size=4.5,force=1,
                             legend_labels = "SMR: cis-eQTL Plasma",legend_position=c(0.9,0.8),
                             legend_text_size = 13,Below_the_line = FALSE){
  
  highlight_names <- highlight_names
  
  Background_data <- Background_data
  colnames(Background_data) <- c("Gene","probeID","CHR","POS","P")
  Background_data[,c(3:5)] <- mutate_all(Background_data[,c(3:5)], as.numeric)
  
  do_manh <- Background_data %>% 
    group_by(CHR) %>% 
    summarise(chr_length = max(POS)) %>%
    mutate(total=cumsum(chr_length)-chr_length) %>%
    dplyr::select(-chr_length) %>%
    left_join(Background_data, ., by=c("CHR"="CHR")) %>%
    arrange(CHR, POS) %>%
    mutate(BPcum=POS+total) %>%
    # Add highlight and annotation information
    mutate(is_highlight=ifelse(probeID %in% highlight_names, "yes", "no"))
  
    # Prepare X axis准备X轴
    x_axis <- do_manh %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum))/2 )
    
    # 准备目标文本标签位置
    label_position <- do_manh %>% 
      group_by(CHR) %>% 
      summarise(chr_pos= min(BPcum))
    
    label_position <- label_position[c(1,4,8,12),]
    high_light <- do_manh[do_manh$is_highlight=="yes",]
  
    high_light$label_position <- 0
    high_light$label_position[high_light$CHR <= 4] <- label_position$chr_pos[1] - high_light$BPcum[high_light$CHR <= 4]
    high_light$label_position[high_light$CHR > 4 & high_light$CHR <= 8] <- label_position$chr_pos[2] - high_light$BPcum[high_light$CHR > 4 & high_light$CHR <= 8] 
    high_light$label_position[high_light$CHR > 8 & high_light$CHR <= 12] <- label_position$chr_pos[3] - high_light$BPcum[high_light$CHR > 8 & high_light$CHR <= 12] 
    high_light$label_position[high_light$CHR > 12] <- label_position$chr_pos[4] - high_light$BPcum[high_light$CHR > 12] 
   
   high_light$lables <- paste0(high_light$probeID," (",high_light$Gene,")")
   
   if(length(grep("pQTL",legend_labels)) > 0){
     high_light$lables <- high_light$probeID
   }
   
   ### 需要调整
   # text_y_pos <- ifelse(max(-log10(do_manh$P)) > 150,-log10(0.05)+(log10(max(-log10(do_manh$P)))/250*(max(-log10(do_manh$P)))),-log10(0.05)+(log10(max(-log10(do_manh$P)))/10*(max(-log10(do_manh$P)))))
     
   position_tag <- 0.001*max(-log10(do_manh$P))*log(max(-log10(do_manh$P)))
   
   if(max(-log10(do_manh$P)) <= 5){
     position_tag <- 0.04*max(-log10(do_manh$P))*log(max(-log10(do_manh$P)))
   }
   if(max(-log10(do_manh$P)) >= 5 & max(-log10(do_manh$P)) < 50 ){
     position_tag <- 0.02*max(-log10(do_manh$P))*log(max(-log10(do_manh$P)))
   }  
   
   position_tag <- ifelse(Below_the_line, -position_tag,position_tag)
   
   text_y_pos <- -log10(0.05) + position_tag
   
   #ggplot2绘图
  ggplot(do_manh, aes(x=BPcum, y=-log10(P))) +
    # Show all points展示点位置
  geom_point(color=points_color, 
             alpha=0.7, 
             size=pointsize,
             show.legend = F,
             shape = pointshape) +
    #scale_color_manual(values = rep(points_color, 44)) +
    # custom X axis:设置x轴
    scale_x_continuous(label = c("CHR",x_axis$CHR), 
                       breaks= c(min(do_manh$BPcum),
                                 x_axis$center)) +
    scale_y_continuous(expand = c(0, 0)) +     # remove space between plot area and x axis
    # Add highlighted points, 新开一组点，展示Legend
    geom_point(data=subset(do_manh, is_highlight=="yes"), 
               aes(color = as.factor(is_highlight)),
               size=highlight_size,
               show.legend = T) +
    scale_color_manual(values = highlight_color,
                       labels = legend_labels)+
    # Add label using ggrepel to avoid overlapping
    geom_label_repel(data=high_light,
                     aes(BPcum,label = lables),
                     nudge_x = high_light$label_position,
                     alpha = 0.9, point.padding = 0.8, box.padding = 1, 
                     direction = "y",max.overlaps = 30,segment.color = 'grey50', ## 指向线色
                     ylim = c(-log10(0.05),0.7*(-log(min(do_manh$P)))),
                     size = label_size,force = force)+
    
    geom_hline(yintercept = -log10(0.05),
               lty = 2,
               color = thresholds_line_color,
               size=0.7) + 
    annotate("text",label = paste0(pval_names," = ",thresholds),
             x=max(do_manh$BPcum)*0.95,
             y=text_y_pos,
             cex = 5,  
             fontface = "italic")+ 
    xlab("")+ 
    ylab(expression(-log[10]~"(p)"))+ 
    #scale_color_discrete(labels = "SMR between cis-eQTL & GWAS")+
    guides(color = guide_legend(title = ""))+ 
    # Custom the theme:#设置主题
    theme_classic() +
    theme(legend.position=legend_position,
          legend.text = element_text(size = legend_text_size,face = "bold"),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.title.y = element_text(size = 15),
          axis.text.x = element_text(face = "bold",size = 10)) 

}

