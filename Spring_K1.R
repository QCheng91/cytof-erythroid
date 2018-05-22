install.packages('DeLorean')
library("cytofkit") 
library(flowCore)
library(tidyverse)
library(stringr)

get_pro_data <- function() {
  
  pro0 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c17_K4+5_03_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE) 
  pro1 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c09_K4+5_03_0_0_FLI1Day4_Ungated_Ungated.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c10_K4+5_03_0_0_FLI1Day6_Ungated_Ungated.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c11_K4+5_03_0_0_FLI1Day8_Ungated_Ungated.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c12_K4+5_03_0_0_FLI1Day10_Ungated_Ungated.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c13_K4+5_03_0_0_FLI1Day11_Ungated_Ungated.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c14_K4+5_03_0_0_FLI1Day12_Ungated_Ungated.fcs", transformation=FALSE)
  pro7 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c16_K4+5_03_0_0_FLI1Day14_Ungated_Ungated.fcs", transformation=FALSE)
  
  #pro7 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_High/c11_K4+5_03_0_0_FLI1Day8_Flag_hi.fcs", transformation=FALSE) 
  
  # pro0 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c01_K4+5_03_0_0_SPi1Day4_Ungated_Ungated.fcs", transformation=FALSE) 
  # pro1 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c01_K4+5_03_0_0_SPi1Day4_Ungated_Ungated.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c02_K4+5_03_0_0_SPi1Day6_Ungated_Ungated.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c03_K4+5_03_0_0_SPi1Day8_Ungated_Ungated.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c05_K4+5_03_0_0_SPi1Day10_Ungated_Ungated.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c06_K4+5_03_0_0_SPi1Day11_Ungated_Ungated.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c07_K4+5_03_0_0_SPi1Day12_Ungated_Ungated.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c08_K4+5_03_0_0_Spi1Day14_Ungated_Ungated.fcs", transformation=FALSE)
  # 
  # pro8 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_Low/c17_K4+5_03_0_0_Day0_Flag_low.fcs", transformation=FALSE)
  # pro9 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_Low/c09_K4+5_03_0_0_FLI1Day4_Flag_low.fcs", transformation=FALSE)
  # pro10 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_Low/c10_K4+5_03_0_0_FLI1Day6_Flag_low.fcs", transformation=FALSE)
  # pro11 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_Low/c11_K4+5_03_0_0_FLI1Day8_Flag_low.fcs", transformation=FALSE)
  # pro12 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_Low/c12_K4+5_03_0_0_FLI1Day10_Flag_low.fcs", transformation=FALSE)
  # pro13 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_Low/c13_K4+5_03_0_0_FLI1Day11_Flag_low.fcs", transformation=FALSE)
  # pro14 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_Low/c14_K4+5_03_0_0_FLI1Day12_Flag_low.fcs", transformation=FALSE)
  # pro15 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Flag_Low/c16_K4+5_03_0_0_FLI1Day14_Flag_low.fcs", transformation=FALSE)
  # 
  pro1
  
  pro0_df = exprs(pro0)
  pro1_df = exprs(pro1)
  pro2_df = exprs(pro2)
  pro3_df = exprs(pro3)
  pro4_df = exprs(pro4)
  pro5_df = exprs(pro5)
  pro6_df = exprs(pro6)
  pro7_df = exprs(pro7)
  # pro8_df = exprs(pro8)
  # pro9_df = exprs(pro9)
  # pro10_df = exprs(pro10)
  # pro11_df = exprs(pro11)
  # pro12_df = exprs(pro12)
  # pro13_df = exprs(pro13)
  # pro14_df = exprs(pro14)
  # pro15_df = exprs(pro15)
  
  g0 = data.frame(pro0_df) 
  g1 = data.frame(pro1_df) 
  g2 = data.frame(pro2_df) 
  g3 = data.frame(pro3_df) 
  g4 = data.frame(pro4_df) 
  g5 = data.frame(pro5_df) 
  g6 = data.frame(pro6_df) 
  g7 = data.frame(pro7_df) 
  # g8 = data.frame(pro8_df) 
  # g9 = data.frame(pro9_df) 
  # g10 = data.frame(pro10_df) 
  # g11 = data.frame(pro11_df) 
  # g12 = data.frame(pro12_df) 
  # g13 = data.frame(pro13_df) 
  # g14 = data.frame(pro14_df) 
  # g15 = data.frame(pro15_df) 
  ###############################################################################
  sur_marker = c("CD34", "CD71", "CD45RA", "CD123", "CD49f", 
                 "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4")
  
  other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "TAL1", "GATA2", 
                  "RUNX1", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "Flag_Tag")
  
  cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")
  
  isotops = c("Sm149Di", "Lu175Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di", "Sm147Di", "Yb171Di",
              
              "Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di", "Nd142Di",
              
              "Nd148Di", "Nd150Di", "Er168Di", "I127Di",
              
              "Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di")
  
  pro_marker = c(sur_marker, other_marker, tran_factor, cell_cycle, DNA)
  
  seq=c(1:37)
  
  day_levels = paste0("Day", c("0", "4", "6", "8", "10", "11", "12", "14"))
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    protein1 = bind_rows( g0[marker] %>% mutate(obstime = "Day0", col_id = 1:nrow(g0)),
                          g1[marker] %>% mutate(obstime = "Day4", col_id = 1:nrow(g1)), 
                          g2[marker] %>% mutate(obstime = "Day6", col_id = 1:nrow(g2)),
                          g3[marker] %>% mutate(obstime = "Day8", col_id = 1:nrow(g3)),
                          g4[marker] %>% mutate(obstime = "Day10", col_id = 1:nrow(g4)), 
                          g5[marker] %>% mutate(obstime = "Day11", col_id = 1:nrow(g5)),
                          g6[marker] %>% mutate(obstime = "Day12", col_id = 1:nrow(g6)),
                          g7[marker] %>% mutate(obstime = "Day14", col_id = 1:nrow(g7))) %>% mutate(obstime = factor(obstime, levels = day_levels))
    
    
    protein1 = protein1 %>% mutate(prot_name = pro_marker[i])
    colnames(protein1) = c("Intensity", "obstime", "cell_id", "cell.type")
    prot_list = c(prot_list, list(protein1))
  }   
  return(prot_list)
}

all_data = as_tibble(bind_rows(get_pro_data()))

#all_data = all_data %>% mutate(obstime = as.numeric(obstime))
all_data = all_data %>% mutate(cell = factor(paste(cell_id, obstime)))
all_data = all_data %>% mutate(capture = obstime)
all_data = all_data %>% mutate(Intensity = asinh(Intensity))
all_data$obstime

# check_data = filter(all_data, all_data$cell.type == "CD34") #all_data$obstime == 20, 
# check_data = filter(all_data, all_data$obstime == "low_0", all_data$cell.type == "CD34")
# check_data

gene_data0 = all_data %>% spread(cell.type, Intensity) 

####################################################################
mean_irid = function(ir1, ir2){
  reg = lm(ir1~ir2)
  # reg = lm(data = gated_data, Ir191~Ir193)
  (predict(reg) + ir1)/2
}

gene_data = gene_data0 %>% filter(Ce140 <= 5.0) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt <= 5.0 | HBA >= 6.5 | CD235ab >= 6.5) %>% 
  filter(mean_Ir >= 5.0 | HBA >= 6.5 | CD235ab >= 6.5,
         mean_Ir <= 8.0 | HBA >= 6.5 | CD235ab >= 6.5)

flag_hi = gene_data %>% filter(Flag_Tag >= 4.5) %>% group_by(obstime) %>% 
  mutate(cell_id = 1:n()) %>% ungroup() %>% mutate(flag = if_else(obstime == "Day0", 0, 1)) %>% gather(cell.type, Intensity, ATRX:TAL1) %>% filter(cell_id<=800)

flag_low = gene_data %>% filter(Flag_Tag < 4.5) %>% group_by(obstime) %>% 
  mutate(cell_id = 1:n()) %>% ungroup() %>% mutate(flag = 0) %>% gather(cell.type, Intensity, ATRX:TAL1) %>% filter(cell_id<=800)

# xx = all_data %>% filter(obstime == "Day14" & cell.type == "Flag_Tag")
# range(xx$Intensity)

all_data %>% spread(cell.type, Intensity) %>% ggplot(aes(x=Flag_Tag, y=PU1)) +
  geom_point(alpha = 0.1) +
  geom_density_2d(aes(colour = "red"))+
  facet_wrap(~obstime,  ncol = 3)
  
gene_data %>% ggplot(aes(x=Flag_Tag)) + 
  geom_histogram(bins=100) +
  facet_wrap(~obstime)
  
gene_data %>% filter(obstime == "Day0") %>% select(Flag_Tag) %>% range(.)  
  ggplot(aes(x = obstime, y = asinh(Intensity), fill = obstime)) +
  geom_violin(size=0.2)  +
  #scale_y_log10() +  
  #stat_summary(fun.y=mean, geom="point", size=2, colour="red") +
  #geom_boxplot(width=0.1) +
  facet_wrap(~cell.type, ncol = 6, scales = "free") 

#expression data
express_data0 = rbind(flag_low, flag_hi) %>% spread(cell.type, Intensity)
express_data = express_data0 %>% select(-Ce140, -cyclin_B1, -H3, -HBB, -IdU, -Ir191, -Ir193, 
                                        -p_HH3, -p_Rb, -Pt194, -Pt195, -Flag_Tag, -mean_Ir, -mean_Pt, -flag)

# express_data %>% ggplot() +
#   geom_jitter(mapping = aes(x=obstime, y = CD41), shape = 16, alpha = 0.5) +
#   geom_smooth(mapping = aes(x=obstime, y = CD41))
# 
# express_data %>% ggplot() +
#   geom_jitter(mapping = aes(x=obstime, y = ratio2), shape = 16, alpha = 0.5)

x = express_data$ratio1
express_data$obstime <- NULL
express_data$cell_id <- NULL
express_data$cell <- NULL
express_data$capture <- NULL

#x1=t(express_data)
x1=express_data
write.table(x1, file='SPRING/datasets/flag_fli1//FLI1_Flag.csv',quote=FALSE, sep='\t',col.names = FALSE, row.names = FALSE)

#cell groupings according to days
day_data = express_data0 #filter_data %>% spread(cell.type, Intensity)
#day_data = day_data %>% mutate(flag = Flag_Tag>= 4.5)
y <- day_data$obstime
z <- as.numeric(day_data$flag)

x2 = data.frame(obstime = y, flag_tag = as.numeric(z))
x2 = t(x2)
rownames(x2) <- c('Day', 'Flag_tag')
#write.csv(x, "K3.csv", col.names = FALSE)
write.table(x2, file='SPRING/datasets/flag_fli1/FLI1_Flag_day.csv',quote=FALSE, sep=',',col.names = FALSE)

all_prot = names(express_data)
#######################################################
#KL Divergence Analysis
install.packages("entropy")
library("entropy")

express_data = all_data %>% spread(cell.type, Intensity)

gated_data= express_data %>% filter(Ce140 <= 5.0) %>% gather(cell.type, Intensity, ATRX:TAL1) 

entropy_disc = function(x, range_vals = c(0, 12), bins = 100)
{
  y1 = discretize(x, numBins=bins, r=range_vals)
  entropy(y1)
}

entro_data = gated_data %>% group_by(obstime, cell.type) %>% summarise(entro = entropy_disc(Intensity)) %>% ungroup() %>% 
  
  ggplot(aes(x = obstime, y = entro)) +
  geom_point(size=2)+
  facet_wrap(~cell.type, ncol = 4, scales = "free") 

gated_data %>% ggplot(aes(x = obstime, y = Intensity, fill=obstime)) +
  geom_violin(size=0.1, scale = "width", trim = FALSE, bw = 0.1)  +
  facet_wrap(~cell.type, ncol = 4, scales = "free") 

figure = gated_data %>% mutate(mean_Pt = (Pt194 + Pt195)/2, mean_Ir = mean_irid(Ir191, Ir193)) %>%
  ggplot(aes(mean_Ir, mean_Pt, colour = HBA)) +
  geom_point(alpha = 0.1)+
  #geom_density_2d(aes(colour = "red"))+
  facet_wrap(~obstime,  ncol = 4)

figure = gated_data %>% mutate(mean_Pt = (Pt194 + Pt195)/2, mean_Ir = mean_irid(Ir191, Ir193)) %>%
  ggplot(aes(p_Rb, CD235ab, colour = mean_Ir)) +
  geom_point(alpha = 0.1)+
  #geom_density_2d(aes(colour = "red"))+
  facet_wrap(~obstime,  ncol = 4)

png(filename = "Desktop/Oct24/Overexpression/test.png", width = 1500, height = 1200, units = "px", pointsize = 20, bg = "white")
plot(figure)
dev.off()

mean_irid = function(ir1, ir2){
  reg = lm(ir1~ir2)
  # reg = lm(data = gated_data, Ir191~Ir193)
  (predict(reg) + ir1)/2
}

