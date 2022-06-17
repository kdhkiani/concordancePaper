theme_set(theme_classic())

gonzalezScatterPlots <- function(deseqtib_path, peaksanno_path, peaksdat_path, joinedTib, outputLoc){
  
  deseqtib    <- read_tsv(deseqtib_path) 
  peaksanno   <- read_tsv(peaksanno_path)
  peaksdat    <- read_tsv(peaksdat_path)
  
  peaktib <- left_join(peaksanno, peaksdat)
  
  peaktib <- peaktib %>% 
    dplyr::rename(CD34_dhs = CD34, CD14_dhs= CD14) %>% 
    dplyr::select(- H1hesc, - CD56, - CD3, -CD19)
  
  deseqtib <- deseqtib %>% 
    dplyr::rename(CD34_rna = CD34, CD14_rna= CD14) %>% 
    dplyr::select(- H1hesc, - CD56, - CD3, -CD19)
  
  bigtib <- left_join(peaktib, deseqtib, by = "symbol")
  
  bigtib <- bigtib %>% 
    dplyr::mutate(peakFilter = grepl(pattern = "CD[1-3]4", x = accessPattern)) %>% 
    dplyr::filter(peakFilter == TRUE)
  
  if (joinedTib) {
    bigtib <- read_tsv(here('extractedData', 'Gonzalez2015', '100KBGonzalez_joinedTablePeaksNearGenes.tsv')) 
    outputLoc <- here('plots','scatter','100kb_')
  }
  
  ##massage and prepare data from bigtib
  
  #add 1 to avoid infinite values when doing log2 fold changes
  bigtib <- bigtib %>% 
    dplyr::mutate(CD34_rna = CD34_rna + 1) %>% 
    dplyr::mutate(CD34_dhs = CD34_dhs + 1) %>% 
    dplyr::mutate(CD14_rna = CD14_rna + 1) %>% 
    dplyr::mutate(CD14_dhs = CD14_dhs + 1)
  
  #annotate max peak value between CD34, CD14 
  bigtib <- bigtib %>% 
    rowwise() %>% 
    dplyr::mutate(maxValDHS = max(CD14_dhs, CD34_dhs)) %>% 
    dplyr::mutate(maxValRNA = max(CD14_rna, CD34_rna))
  
  bigtib <- bigtib %>% 
    dplyr::filter(maxValRNA > 5)
  
  ##massage and prepare data from bigtib
  
  #calculate log2FC for tpm values between CD34,CD14 and create a filtering column
  bigtib <- bigtib %>% 
    dplyr::mutate(log2FC = log2(CD14_rna/CD34_rna)) 
  #calculate log2FC for accessibility and create a filtering column
  bigtib <- bigtib %>% 
    dplyr::mutate(log2FCdhs = log2(CD14_dhs/CD34_dhs))
  
  #calculate median change in DHS 
  plotTib <- bigtib %>% 
    group_by(symbol) %>% 
    dplyr::mutate(medianDHSlog2FC = median(log2FCdhs))
  
  plotTib <- plotTib %>% 
    dplyr::select(log2FC, medianDHSlog2FC) %>% 
    unique()
  
   p1 <- ggscatter(plotTib, x = "medianDHSlog2FC", y = "log2FC", color = "#ee8833", alpha = 0.3, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
    stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
    ggtitle("Median fold change Acessibility vs. fold change Expression CD14 vs CD34") +
    xlab("log2 FC Accessibility") +
    ylab("log2 FC Expression")
  
   ggsave(filename = paste0(outputLoc, "gonzalezMedianScatter.svg"), plot = p1)
   
   #calculate bootstrapped R values
   bootTib <- plotTib %>% 
     dplyr::select("medianDHSlog2FC", "log2FC") %>% 
     na.omit()
   
   boot <- boot(bootTib, statistic = function(bootTib, i){
     cor(bootTib[i, "medianDHSlog2FC"], bootTib[i, "log2FC"], method = "pearson")
   },
   R=10000)
   
   boot.ci(boot, type = 'basic') 
   
   ##look at max value change in DHS
   plotTib <- bigtib %>% 
     group_by(symbol) %>% 
     dplyr::mutate(maxDHSlog2FC = max(log2FCdhs))
   
   plotTib <- plotTib %>% 
     dplyr::select(log2FC, maxDHSlog2FC) %>% 
     unique()
   
   p2 <- ggscatter(plotTib, x = "maxDHSlog2FC", y = "log2FC", color = "#ee8833", alpha = 0.3, shape = 16, size = 2, add="reg.line", cor.method = "pearson") +
     stat_cor(label.x = 0, label.y = 11) + stat_regline_equation(label.x = 0,label.y = 12) + 
     ggtitle("Maximum fold change Acessibility vs. fold change Expression CD14 vs CD34") +
     xlab("log2 FC Accessibility") +
     ylab("log2 FC Expression")
   
   
   ggsave(filename = paste0(outputLoc, "gonzalezMaxScatter.svg"), plot = p1)
   
   #calculate bootstrapped R values
   bootTib <- plotTib %>% 
     dplyr::select("maxDHSlog2FC", "log2FC") %>% 
     na.omit()
   
   boot <- boot(bootTib, statistic = function(bootTib, i){
     cor(bootTib[i, "maxDHSlog2FC"], bootTib[i, "log2FC"], method = "pearson")
   },
   R=10000)
   
   boot.ci(boot, type = 'basic') 
 
}
 