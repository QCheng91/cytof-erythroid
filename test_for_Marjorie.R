#####Plot for the comparison between datasets#####

library(tidyverse)
time = c(4, 6, 8, 10, 11, 12, 14)

Flag_high = data.frame(Flag_High_total = c(3699, 6411, 6367, 1321, 641, 158, 288), 
                       Flag_High_cd41 = c(82, 465, 595, 300, 213, 56, 145))

Flag_low = data.frame(Flag_Low_total = c(88568, 52694, 134740, 169970, 202883, 115542, 264530),
                       Flag_Low_cd41 = c(1345, 1199, 2185, 2843, 3222, 1578, 3604))

Flag_control = data.frame(Flag_control_total = c(49578, 52622,
                                                 40917,
                                                 51722,
                                                 71893,
                                                 28932,
                                                 43465),
                                                 
                          Flag_control_cd41 = c(844,
                                                864,
                                                876, 
                                                1207,
                                                1264,
                                                882,
                                                974))

Flag_total = cbind(Flag_high, Flag_low) %>% mutate(fli1_total = Flag_High_total+Flag_Low_total, fli1_cd41= Flag_High_cd41+Flag_Low_cd41) %>% 
  mutate(ratio =fli1_cd41/fli1_total, Flag_FLI1 = "Total") %>% select(ratio, Flag_FLI1) %>% mutate(time = time) 

Flag_high = Flag_high %>% mutate(ratio =Flag_High_cd41/Flag_High_total, Flag_FLI1 = "Flag+") %>% select(ratio, Flag_FLI1) %>% mutate(time = time) 
Flag_low  = Flag_low %>% mutate(ratio =Flag_Low_cd41/Flag_Low_total, Flag_FLI1 = "Flag-") %>% select(ratio, Flag_FLI1) %>% mutate(time = time) 


Flag_control  = Flag_control %>% mutate(ratio =Flag_control_cd41/Flag_control_total, flag = "pTRIP") %>% select(ratio, flag) %>% mutate(time = time) 

library(ggforce)

rbind(Flag_high, Flag_low) %>% as.tibble() %>% mutate(type = factor(Flag_FLI1)) %>% 
  ggplot(aes(x = time, y = ratio*100))+
  geom_point(aes(colour = Flag_FLI1), size =4) +
  ylab("Percentage of cells in Mega branch (%)") +
  xlab("Day") +
  theme_grey(base_size = 15)
  facet_zoom(y = ratio >=0 & ratio <= 0.04)

data_x = gated_data %>%
  as.tibble() %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% 
  left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 

data_x1 = data_x %>% select(-Day, -day_index) %>% as.matrix() #Ir191, Ir193, IdU, p_Rb

day_indx = data_x %>%  pull(day_index)

all_cells = predict_cluster(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx) # 

data_x %>% mutate(cluster = all_cells) %>% group_by(cluster) %>% summarise(count = n())

rbind(Flag_high, Flag_low, Flag_total, Flag_control) %>% as.tibble() %>% mutate(type = factor(flag)) %>% 
  ggplot(aes(x = time, y = ratio*100))+
  geom_point(aes(colour = flag), size =3.5) +
  ylab("Percentage of cells in Mega branch") +
  facet_zoom(y = ratio >=0 & ratio <= 0.04)

coordinates = read.table(file = "Results_by_date/Dec2018/May/K3/coordinates.txt", header = FALSE, sep = ",", dec = ".")

express_data = read.table(file='Software/SPRING_new_version_master/datasets/May_K3/prot_data.csv', header = FALSE, sep='\t')
class = read.table(file='Software/SPRING_new_version_master/datasets/May_K3/classification.csv',  header = FALSE, sep=',') %>% as.tibble() %>% 
  select(-V1) %>% as.matrix() %>% 
  t() %>% as.tibble()

all_prot = read.table(file='Software/SPRING_new_version_master/datasets/May_K3/prot_name.csv', sep='\t') %>% t()

colnames(express_data) <- all_prot
classify = class$V2 %>% as.numeric()

coord_x = coordinates$V2
coord_y = coordinates$V3
  
express_data %>% mutate(flag = as.factor(classify), x = coord_x, y = coord_y) %>% 
  filter(!flag %in% c(7, 11, 13, 4, 19, 2)) %>% 
  #pull(flag) %>% unique()
  group_by(flag) %>% 
 # summarise(count = n())
  sample_n(100) %>% 
   ungroup() %>% 
  ggplot(aes(x = y , y = -x, fill = flag))+
  geom_point(shape = 21, size = 8.5, stroke = 0.8) +
  #facet_wrap(~flag) +
  theme_void() +
  #scale_fill_brewer(palette="Spectral")
  scale_fill_manual(values=c("yellow1", "gold1", "burlywood", "tan2", "salmon2",
                             "plum2", "blue2", "navyblue", "blueviolet", "red",
                             "darkmagenta", "red2", "red3"))


  scale_color_gradient2(mid="black",
                        high="green", space ="Lab" )
  
##############################################################################
#install.packages('DeLorean')
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

  # pro0 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c17_K4+5_03_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  # pro1 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c01_K4+5_03_0_0_SPi1Day4_Ungated_Ungated.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c02_K4+5_03_0_0_SPi1Day6_Ungated_Ungated.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c03_K4+5_03_0_0_SPi1Day8_Ungated_Ungated.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c05_K4+5_03_0_0_SPi1Day10_Ungated_Ungated.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c06_K4+5_03_0_0_SPi1Day11_Ungated_Ungated.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c07_K4+5_03_0_0_SPi1Day12_Ungated_Ungated.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-SPI1/Ungated/c08_K4+5_03_0_0_Spi1Day14_Ungated_Ungated.fcs", transformation=FALSE)
  
  # pro0 <- read.FCS("Research/RBC/Data Analysis/overexpression-Kinetics_10_17/Flag-FLI1/Ungated/c17_K4+5_03_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  # pro1 <- read.FCS("Backup/cytof-erythroid/Overexpression-Results/6-Exp-Overexpression-pTRIPcontrol-Kinetic3-14-9-2017/Ungated/c11_K3_06_0_pTRIPDay4_Ungated_Ungated.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Backup/cytof-erythroid/Overexpression-Results/6-Exp-Overexpression-pTRIPcontrol-Kinetic3-14-9-2017/Ungated/c12_K3_06_0_pTRIPDay6_Ungated_Ungated.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Backup/cytof-erythroid/Overexpression-Results/6-Exp-Overexpression-pTRIPcontrol-Kinetic3-14-9-2017/Ungated/c13_K3_06_0_pTRIPDay8_Ungated_Ungated.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Backup/cytof-erythroid/Overexpression-Results/6-Exp-Overexpression-pTRIPcontrol-Kinetic3-14-9-2017/Ungated/c14_K3_06_0_pTRIPDay10_Ungated_Ungated.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Backup/cytof-erythroid/Overexpression-Results/6-Exp-Overexpression-pTRIPcontrol-Kinetic3-14-9-2017/Ungated/c15_K3_06_0_pTRIPDay11_Ungated_Ungated.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Backup/cytof-erythroid/Overexpression-Results/6-Exp-Overexpression-pTRIPcontrol-Kinetic3-14-9-2017/Ungated/c16_K3_06_0_pTRIPDay12_Ungated_Ungated.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Backup/cytof-erythroid/Overexpression-Results/6-Exp-Overexpression-pTRIPcontrol-Kinetic3-14-9-2017/Ungated/c20_K3_06_0_pTRIPDay14_Ungated_Ungated.fcs", transformation=FALSE)

  pro0_df = flowCore::exprs(pro0)
  pro1_df = flowCore::exprs(pro1)
  pro2_df = flowCore::exprs(pro2)
  pro3_df = flowCore::exprs(pro3)
  pro4_df = flowCore::exprs(pro4)
  pro5_df = flowCore::exprs(pro5)
  pro6_df = flowCore::exprs(pro6)
  pro7_df = flowCore::exprs(pro7)
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
                 "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4",
                 "CD36", "CD38")
  
  other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "TAL1", "GATA2",  #"PU1"
                  "RUNX1", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "Flag_Tag",
                  "FLI1", "NFE2p45")
  
  cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")
  
  isotops = c("Sm149Di", "Lu175Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di", "Sm147Di", "Yb171Di",
              "Gd155Di", "Yb172Di",
              
              "Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di", "Nd142Di",
              "Tm169Di", "Sm154Di",
              
              "Nd148Di", "Nd150Di", "Er168Di", "I127Di",
              
              "Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di")
  
  pro_marker = c(sur_marker, other_marker, tran_factor, cell_cycle, DNA)
  
  seq=c(1:41)
  
  day_levels = c("0", "4", "6", "8", "10", "11", "12", "14")
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    protein1 = bind_rows( g0[marker] %>% mutate(obstime = "0", col_id = 1:nrow(g0)),
                          g1[marker] %>% mutate(obstime = "4", col_id = 1:nrow(g1)), 
                          g2[marker] %>% mutate(obstime = "6", col_id = 1:nrow(g2)),
                          g3[marker] %>% mutate(obstime = "8", col_id = 1:nrow(g3)),
                          g4[marker] %>% mutate(obstime = "10", col_id = 1:nrow(g4)), 
                          g5[marker] %>% mutate(obstime = "11", col_id = 1:nrow(g5)),
                          g6[marker] %>% mutate(obstime = "12", col_id = 1:nrow(g6)),
                          g7[marker] %>% mutate(obstime = "14", col_id = 1:nrow(g7))) %>% mutate(obstime = factor(obstime, levels = day_levels))
    
    
    protein1 = protein1 %>% mutate(prot_name = pro_marker[i])
    colnames(protein1) = c("Intensity", "obstime", "cell_id", "cell.type")
    prot_list = c(prot_list, list(protein1))
  }   
  return(prot_list)
}

get_pro_data <- function() {
  
  #K1
  pro1 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c01_K1_02_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c02_K1_02_0_0_Day2_Ungated_Ungated.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c03_K1_02_0_0_Day4_Ungated_Ungated.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c04_K1_02_0_0_Day6_Ungated_Ungated.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c05_K1_02_0_0_Day8_Ungated_Ungated.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c06_K1_02_0_0_Day10_Ungated_Ungated.fcs", transformation=FALSE)
  pro7 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c07_K1_02_0_0_Day11_Ungated_Ungated.fcs", transformation=FALSE)
  pro8 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c08_K1_02_0_0_Day12_Ungated_Ungated.fcs", transformation=FALSE)
  pro9 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c09_K1_02_0_0_Day14_Ungated_Ungated.fcs", transformation=FALSE)
  pro10 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c10_K1_02_0_0_Day16_Ungated_Ungated.fcs", transformation=FALSE)
  pro11 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c15_K1_02_0_0_Day18_Ungated_Ungated.fcs", transformation=FALSE)
  pro12 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c16_K1_02_0_0_Day20_Ungated_Ungated.fcs", transformation=FALSE)
  pro13 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Controls/c19_K1_02_0_0_MNCs_Ungated_Ungated.fcs", transformation=FALSE)
  pro14 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Controls/c17_K1_02_0_0_Jurkat_Ungated_Ungated.fcs", transformation=FALSE)

  # # #K2
  # pro1 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c01_K2_01_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c02_K2_01_0_0_Day2_Ungated_Ungated.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c03_K2_01_0_0_Day4_Ungated_Ungated.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c04_K2_01_0_0_Day6_Ungated_Ungated.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c05_K2_01_0_0_Day8_Ungated_Ungated.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c06_K2_01_0_0_Day10_Ungated_Ungated.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c07_K2_01_0_0_Day11_Ungated_Ungated.fcs", transformation=FALSE)
  # pro8 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c08_K2_01_0_0_Day12_Ungated_Ungated.fcs", transformation=FALSE)
  # pro9 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c09_K2_01_0_0_Day14_Ungated_Ungated.fcs", transformation=FALSE)
  # pro10 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c10_K2_01_0_0_Day16_Ungated_Ungated.fcs", transformation=FALSE)
  # pro11 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c15_K2_01_0_0_Day18_Ungated_Ungated.fcs", transformation=FALSE)
  # pro12 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c20_K2_01_0_0_Day20_Ungated_Ungated.fcs", transformation=FALSE)

  #tibble(desc = as.character(pro13@parameters@data$desc), name = names(pro13@parameters@data$desc)) %>% separate(desc, c("Isotope", "Protein"), "_") %>% View()
  
  pro1_df = exprs(pro1)
  pro2_df = exprs(pro2)
  pro3_df = exprs(pro3)
  pro4_df = exprs(pro4)
  pro5_df = exprs(pro5)
  pro6_df = exprs(pro6)
  pro7_df = exprs(pro7)
  pro8_df = exprs(pro8)
  pro9_df = exprs(pro9)
  pro10_df = exprs(pro10)
  pro11_df = exprs(pro11)
  pro12_df = exprs(pro12)
  # pro13_df = exprs(pro13)
  # pro14_df = exprs(pro14)
  
  g1 = data.frame(pro1_df)
  g2 = data.frame(pro2_df)
  g3 = data.frame(pro3_df)
  g4 = data.frame(pro4_df)
  g5 = data.frame(pro5_df)
  g6 = data.frame(pro6_df)
  g7 = data.frame(pro7_df)
  g8 = data.frame(pro8_df)
  g9 = data.frame(pro9_df)
  g10 = data.frame(pro10_df)
  g11 = data.frame(pro11_df)
  g12 = data.frame(pro12_df)
  # g13 = data.frame(pro13_df)
  # g14 = data.frame(pro14_df)
  ###############################################################################
  sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f", 
                 "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4")
  
  other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                  "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "MEF2C")
  
  cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")
  
  
  isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di", "Sm147Di", "Yb171Di",
              
              "Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Sm154Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di", "Nd146Di",
              
              "Nd148Di", "Nd150Di", "Er168Di", "I127Di",
              
              "Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di") 
  
  pro_marker = c(sur_marker, other_marker, tran_factor, cell_cycle, DNA)
  
  seq=c(1:41)
  
  #day_levels = paste0("", c("MNCs", "Jurkats")) #paste0("", c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16, 18 ,20, "MNCs", "Jurkats"))
  day_levels = paste0("", c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16, 18 ,20))
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    protein1 = bind_rows( g1[marker] %>% mutate(Day = "0", col_id = 1:nrow(g1)), 
                          g2[marker]%>% mutate(Day = "2", col_id = 1:nrow(g2)),
                          g3[marker]%>% mutate(Day = "4", col_id = 1:nrow(g3)),
                          g4[marker]%>% mutate(Day = "6", col_id = 1:nrow(g4)),
                          g5[marker]%>% mutate(Day = "8", col_id = 1:nrow(g5)),
                          g6[marker]%>% mutate(Day = "10", col_id = 1:nrow(g6)),
                          g7[marker]%>% mutate(Day = "11", col_id = 1:nrow(g7)),
                          g8[marker]%>% mutate(Day = "12", col_id = 1:nrow(g8)),
                          g9[marker]%>% mutate(Day = "14", col_id = 1:nrow(g9)),
                          g10[marker]%>% mutate(Day = "16", col_id = 1:nrow(g10)), 
                          g11[marker]%>% mutate(Day = "18", col_id = 1:nrow(g11)),
                          g12[marker]%>% mutate(Day = "20", col_id = 1:nrow(g12))) %>% 
      # g13[marker]%>% mutate(Day = "MNCs", col_id = 1:nrow(g13)),
      # g14[marker]%>% mutate(Day = "Jurkats", col_id = 1:nrow(g14))) %>%
      mutate(Day = factor(Day, levels = day_levels))
    
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

gene_data0 = all_data %>% spread(cell.type, Intensity) 
####################################################################
mean_irid = function(ir1, ir2){
  reg = lm(ir1~ir2)
  # reg = lm(data = gated_data, Ir191~Ir193)
  (predict(reg) + ir1)/2
}

gene_data = gene_data0 %>% filter(Ce140 <= sinh(5.0)) %>% 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) %>% 
  filter(mean_Ir >= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5),
         mean_Ir <= sinh(8.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

gene_data %>% group_by(obstime) %>% summarise(count = n())

gene_data %>% ggplot(aes(x = asinh(CD41)))+
  geom_histogram() + 
  facet_wrap(~obstime)

gene_data %>% as.tibble() %>% filter(obstime != 0) %>% 
  #filter(Flag_Tag < sinh(5.0)) %>%
  filter(CD41 > sinh(4.5)) %>% 
  group_by(obstime) %>% 
  summarise(count=n())

flag_day_0 = gene_data %>% as.tibble() %>% filter(obstime == 0) %>% 
  filter(CD34 >= asinh(5)) %>% sample_n(1000) %>% 
  ungroup() %>% 
  mutate(flag = 0)

flag_hi_1 = gene_data %>% as.tibble() %>% filter(obstime != 0, obstime != 11, obstime != 12, obstime != 14) %>% 
  filter(Flag_Tag >= sinh(5.0)) %>% group_by(obstime) %>% 
  sample_n(1000) %>% 
  ungroup() %>% 
  mutate(flag = 1)

flag_hi_2 = gene_data %>% as.tibble() %>% filter(obstime == 11 | obstime == 12 | obstime == 14) %>% 
  filter(Flag_Tag >= sinh(5.0)) %>% 
  mutate(flag = 1)

flag_low = gene_data %>% as.tibble() %>% filter(obstime != 0) %>% 
  filter(Flag_Tag < sinh(5.0)) %>% group_by(obstime) %>% 
  sample_n(1000) %>% 
  ungroup() %>% 
  mutate(flag = 0)

express_data0 = rbind(flag_day_0, flag_low, flag_hi_1, flag_hi_2)
########################
oct_k1 = gene_data %>% select(markers) # select(-Day)

markers = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f", 
            "CD90", "CD44", "CD41", "CD235ab",
            "GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
            "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B")

data_x = gene_data %>%
  as.tibble() %>% 
  mutate(Day = as.numeric(as.character(obstime))) %>%
  select(Day, cell, markers) %>% 
  gather(prot_name, Intensity, markers) %>% 
  spread(prot_name, Intensity)
  left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 

data_x1 = data_x %>% select(ATRX:TAL1) %>% as.matrix() #Ir191, Ir193, IdU, p_Rb

day_indx = data_x %>%  pull(day_index)

all_cells = predict_cluster(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx) #

#print(table(all_cells)) %>% sum()
cluster_prob = all_cells %>% as.tibble()

full_clusters = all_cells %>% as.tibble() %>% unique() %>% pull(V1)

###using our cluster methods to find the cells in pTRIP group

all_data %>% group_by(Day) %>% sample_n(5000) %>% ungroup() %>%
  filter(prot_name %in% c("DNA-1", "Oct4", "Klf4", "c-Myc", "Sox2", "CyclinB1", "pH3", "IdU", "Ki67", "Nanog",
                          "CD54", "CD73", "Lin28", "CD24", "CD140a")) %>% 
  mutate(Intensity = asinh(Intensity)) %>%
  ggplot(aes(x = Day, y = Intensity, fill = Day))+
  geom_violin(size=0.2, scale = "width", trim = FALSE, bw = 0.5)  +
  facet_wrap(~prot_name, ncol = 5) 

express_data = all_data %>% mutate(Intensity = if_else(Intensity < 0, 0, Intensity)) %>%
  spread(prot_name, Intensity)

express_data %>% ggplot(aes(x = asinh(IdU)))+
  geom_histogram(bins = 100)+
  facet_wrap(~Day) #, scales = "free")

express_data %>% group_by(Day) %>% sample_n(5000) %>% ungroup() %>% 
  ggplot(aes(x = asinh(Ki67), y = asinh(IdU)))+
  geom_point(alpha = 0.25) #+
facet_wrap(~Day)

# cell_data = express_data %>% filter(Ce140 <= sinh(5.0)) %>%
#   mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193))
# 
# gated_data1 = cell_data %>% filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

all_data_1 = express_data %>% select(`DNA-1`, IdU, Ki67, Cell_id, Day)
#####################################################################
sur_marker = c("CD34", "CD71", "CD45RA", "CD123", "CD49f", 
               "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4",
               "CD36", "CD38")

other_marker = c("HBB", "HBA", "H3")

tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "TAL1", "GATA2",
                "RUNX1", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "Flag_Tag",
                "FLI1", "NFE2p45")

cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")

DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")

gated_data = express_data0 %>%
  mutate(Day = obstime) %>% 
  select(-cell_id, -cell, -capture, -obstime) %>% 
  select(-cell_cycle, -DNA, -other_marker, -mean_Pt, -mean_Ir) %>% 
  mutate_if(is.numeric, ~ floor(.x)) %>% 
  filter((!Day %in% c("MNCs", "Jurkats")))

prot_name = gated_data %>% select(-Day, -flag) %>% colnames() %>% as.character()

prot_list = data.frame(prot_name = prot_name, prot_index = 1:length(prot_name)) %>% as.tibble() %>% 
  mutate(prot_name = as.character(prot_name))

gated_data %>% group_by(Day) %>% summarise(count=n())

cells_per_day = 5000
abv_max = gated_data %>% group_by(Day) %>% filter(n()<=5000)
sub_data = gated_data %>% group_by(Day) %>% filter(n()>5000) %>% sample_n(cells_per_day) %>% ungroup() %>%
  bind_rows(abv_max) %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% 
  arrange(Day) 

sub_data = sub_data %>% group_by(Day) %>% sample_n(cells_per_day) %>% ungroup()

x_data   = sub_data %>% select(-Day) %>% as.matrix()
vec_d    = sub_data %>% mutate(new_indx = as.numeric(factor(as.character(Day), levels = unique(as.numeric(Day))))) %>% pull(new_indx)
indx_tp  = sub_data %>% group_by(Day) %>% summarise(num_cells = n()) %>% mutate(num_cells  = cumsum(num_cells)) %>% 
  mutate(from = c(0, num_cells[-n()]) + 1) %>% select(from, num_cells) %>% as.matrix()

x_indx = x_data %>% as_tibble() %>% mutate_all(~as.numeric(as.factor(.x)))
x_unique_aux = x_data %>% as_tibble() %>% map2(1:ncol(.), ~ as.factor(.x) %>% levels() %>% 
                                                 tibble(unique_count = ., prot =.y)) %>% bind_rows()
unique_prot_x = x_unique_aux %>% mutate(unique_count = as.numeric(unique_count)) %>% pull(unique_count)
indx_prot  = x_unique_aux %>% group_by(prot) %>% summarise(num_vals = n()) %>% mutate(num_vals  = cumsum(num_vals)) %>% 
  mutate(from = c(0, num_vals[-n()]) + 1) %>% select(from, num_vals) %>% as.matrix()

data_cyto = list(c = 20, k = ncol(x_data), N = nrow(x_data), x = x_data, d = max(vec_d), 
                 day_index = vec_d, indx_tp = indx_tp,
                 num_unique_counts = length(unique_prot_x),
                 unique_prot_x = unique_prot_x,
                 indx_prot = indx_prot,
                 x_indx = x_indx)

vec_d %>% unique()

library(rstan)

options(mc.cores = parallel::detectCores())
doMC::registerDoMC(cores = 14)
rstan_options(auto_write = TRUE)
# map_cyto = optimizing(model2, data = data_cyto)
model   = stan_model("~/mix_nb_model_Ed_vec_V4.stan") #V4
model   = stan_model("~/Software/cluster_model/mix_nb_model_Ed_vec_V4.stan") 

#vb a function in stan model
runs = c(1:4)

library(doParallel)

foreach(number = runs)%dopar%{
  start_time <- Sys.time()
  map_est4   = optimizing(model, data = data_cyto, iter = 3000, as_vector = F)
  #map_est4   = vb(model, data = data_cyto, algorithm = "fullrank")
  #print(start_time - Sys.time())
  
  direct = paste("Results_by_date/Nov2018/Oct_data/", "Flag_pTRIP_", number, ".rda", sep="") #Desktop/Sep2018/Controls_Oct/Oct_K1/", "K1_mncs_nn_"
  save(map_est4, file = direct)
  
  return(start_time - Sys.time())
}

save(sub_data, file = "Desktop/Sep2018/pheno_cluster/K1_500_in_5000.rda")

###############################################################################
load(file = "Research/Data_Analysis/Over-expressed/cluster_model/FLI1_4.rda")
load(file = "Results_by_date/Nov2018/Oct_data/Flag_Fli1_balance_1.rda")
#########analyse the data#########################
c = 20

day_vals = sub_data$Day %>% unique()
day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)

run_fit = map_est4

cluster_plot <- function(run_fit){
  
  vals = run_fit$par 
  
  time = tibble(Day = day_index_0$day_index,  real_time= day_index_0$Day) 
  
  cluster = vals$theta %>% as.tibble() %>% #t()
    pull(value) %>% as.tibble() %>%
    rename(vals = value) %>%
    mutate(cluster = 1:c)
  
  # cluster = vals$theta %>% t() %>% as.tibble() %>%
  # rename(vals = V1) %>%
  # mutate(cluster = 1:c)
  
  phi = vals$phi_s_cp %>% as.tibble() %>% 
    mutate(cluster = 1:c) %>% 
    gather(protein, std, -cluster) %>%
    mutate(protein = as.double(str_remove_all(protein, "V")))
  
  signal = vals$mu_s_cp %>% as.tibble() %>% 
    mutate(cluster = 1:c) %>% 
    gather(protein, vals, -cluster) %>%
    mutate(protein = as.numeric(str_remove_all(protein, "V"))) %>% 
    left_join(prot_list, by = c("protein" = "prot_index")) %>% 
    left_join(phi) 
  
  signal = signal %>% 
    group_by(prot_name) %>% 
    mutate(max_vals = max(vals), min_vals = min(vals), mid_vals = mean(vals), sd_vals = sd(vals)) %>% 
    ungroup()
  
  day_vals = gated_data$Day %>% unique()
  day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)
  
  mu_c_fit = vals$mu_s_cp %>% t()
  
  theta0_c_fit = vals$theta0_s_cp %>% t()
  
  phi_c_fit = vals$phi_s_cp %>% t()
  
  mu_bg_fit = vals$mu_bg %>% as.vector()
  
  st_fit = vals$scaling_t
  st_fit = c(1, st_fit)
  
  theta_fit = rep(cluster$vals, length(day_vals)) %>% 
    matrix(nrow = c)
  
  data_x = gated_data %>%
    as.tibble() %>% 
    mutate(Day = as.numeric(as.character(Day))) %>% 
    left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 
  
  data_x1 = data_x %>% select(-Day, -day_index, -flag) %>% as.matrix() #Ir191, Ir193, IdU, p_Rb
  
  day_indx = data_x %>%  pull(day_index)
  
  all_cells = predict_cluster(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx) # 
  
  #print(table(all_cells)) %>% sum()
  cluster_prob = all_cells %>% as.tibble()
  
  full_clusters = all_cells %>% as.tibble() %>% unique() %>% pull(V1)
  
  # p1 = signal %>% mutate(valss = (vals-min_vals)/(max_vals-min_vals)) %>% 
  #   filter(cluster %in% full_clusters) %>% 
  #   ggplot(aes(x = prot_name, y= vals)) + 
  #   geom_point(aes(colour = prot_name)) +
  #   scale_y_log10()+
  #   geom_text(aes(label = prot_name), size=3.0)+
  #   facet_wrap(~cluster, nrow = 1) +
  #   theme(legend.position = "none")
  key_proteins = c("DNA-1", "Oct4", "Klf4", "Sox2", "Ki67", "Nanog", "CD54", "CD73", "CD24", "CD140a")
  
  p1 = signal %>% mutate(valss = (vals-mid_vals)/(sd_vals)) %>% 
    #filter(cluster %in% full_clusters) %>% 
    #filter(prot_name %in% key_proteins) %>% 
    ggplot(aes(x = prot_name, y= vals)) + 
    geom_point(aes(colour = prot_name)) +
    scale_y_log10()+
    geom_text(aes(label = prot_name), size=3.0)+
    facet_wrap(~cluster, nrow = 2) +
    theme(legend.position = "none")
  
  p2 = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
    mutate(prob = num/sum(num)) %>% ungroup() %>% 
    select(-num) %>% 
    spread(cluster, prob, fill = 0) %>% 
    gather(cluster, prob, -Day) %>% 
    mutate(cluster = as.numeric(cluster)) %>% 
    ggplot(aes(x = Day, y = prob))+
    geom_point()+
    facet_wrap(~cluster, nrow=2)
  
  gridExtra::grid.arrange(p1, p2, ncol=1)
}

cluster_plot(map_est4)

##############SPRING plot data generation###################
express_data = data_x %>% mutate(cluster = all_cells)

express_data$Day <- NULL
express_data$day_index <- NULL
express_data$cluster <- NULL
express_data$Flag_Tag <- NULL
express_data$flag <- NULL
express_data$FLI1 <- NULL

x1=express_data %>% asinh()
write.table(x1, file='Software/SPRING_new_version_master/datasets/Flag_FLI1/prot_data.csv',quote=FALSE, sep='\t',col.names = FALSE, row.names = FALSE)

#cell groupings according to days
day_data = data_x %>% mutate(cluster = all_cells)
y <- day_data$Day
z <- day_data$cluster
m <- day_data$flag

x2 = data.frame(obstime = y, cluster = z, flag = m)
x2 = t(x2)
rownames(x2) <- c('Day', 'Cluster', 'Flag') 
write.table(x2, file='Software/SPRING_new_version_master/datasets/Flag_FLI1/classification.csv',quote=FALSE, sep=',',col.names = FALSE)

all_prot = names(express_data)
write(all_prot, file='Software/SPRING_new_version_master/datasets/Flag_FLI1/prot_name.csv', ncolumns = length(all_prot), sep='\t')

flag_data = data_x %>% select(Flag_Tag, CD41, CD235ab, flag) %>% asinh()
write.table(flag_data, file='Software/SPRING_new_version_master/datasets/Flag_FLI1/flag.csv', sep='\t', col.names = TRUE, row.names = FALSE)

############New plot for May Data K3#############################
library(flowCore)
library(tidyverse)
library(stringr)

get_pro_data <- function() {
  
  #K1
  pro1 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c01_sample-1_01_0_Day2_live.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c02_sample-1_01_0_Day4_live.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c03_sample-1_01_0_Day6_live.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c04_sample-1_01_0_Day8_live.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c05_sample-1_01_0_Day10_live.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c06_sample-1_01_0_Day11_live.fcs", transformation=FALSE)
  pro7 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c16_sample-1_01_0_Day12_live.fcs", transformation=FALSE)
  pro8 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c17_sample-1_01_0_Day14_live.fcs", transformation=FALSE)
  pro9 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c18_sample-1_01_0_Day16_live.fcs", transformation=FALSE)
  pro10 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c19_sample-1_01_0_Day18_live.fcs", transformation=FALSE)
  pro11 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c20_sample-1_01_0_Day20_live.fcs", transformation=FALSE)
  
  #K2
  # pro1 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c01_sample_01_0_Day0_live.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c02_sample_01_0_Day2_live.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c03_sample_01_0_Day4_live.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c04_sample_01_0_Day6_live.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c05_sample_01_0_Day8_live.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c07_sample_01_0_Day10_live.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c08_sample_01_0_Day11_live.fcs", transformation=FALSE)
  # pro8 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c09_sample_01_0_Day12_live.fcs", transformation=FALSE)
  # pro9 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c10_sample_01_0_Day14_live.fcs", transformation=FALSE)
  # pro10 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c11_sample_01_0_Day16_live.fcs", transformation=FALSE)
  # pro11 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c12_sample_01_0_Day18_live.fcs", transformation=FALSE)
  # pro12 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c13_sample_01_0_Day20_live.fcs", transformation=FALSE)
  # pro13 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K3/export_c14_sample_01_0_Day22_live.fcs", transformation=FALSE)
  
  pro1_df = flowCore::exprs(pro1)
  pro2_df = flowCore::exprs(pro2)
  pro3_df = flowCore::exprs(pro3)
  pro4_df = flowCore::exprs(pro4)
  pro5_df = flowCore::exprs(pro5)
  pro6_df = flowCore::exprs(pro6)
  pro7_df = flowCore::exprs(pro7)
  pro8_df = flowCore::exprs(pro8)
  pro9_df = flowCore::exprs(pro9)
  pro10_df = flowCore::exprs(pro10)
  pro11_df = flowCore::exprs(pro11)
  # pro12_df = flowCore::exprs(pro12)
  # pro13_df = flowCore::exprs(pro13)
  
  g1 = data.frame(pro1_df) 
  g2 = data.frame(pro2_df) 
  g3 = data.frame(pro3_df) 
  g4 = data.frame(pro4_df) 
  g5 = data.frame(pro5_df) 
  g6 = data.frame(pro6_df) 
  g7 = data.frame(pro7_df) 
  g8 = data.frame(pro8_df) 
  g9 = data.frame(pro9_df) 
  g10 = data.frame(pro10_df) 
  g11 = data.frame(pro11_df)
  # g12 = data.frame(pro12_df)
  # g13 = data.frame(pro13_df)
  
  ###############################################################################
  sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f",
                 "CD90", "CD44", "CD41", "CD235ab")
  
  #other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                  "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B") #, "BCL11a")
  
  #cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  #DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194", "Event_Length")
  
  isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di",
              
              #"Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Sm154Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di") #, "Sm147Di")
  
  #"Nd148Di", "Nd150Di", "Er168Di", "I127Di",
  
  #"Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di", "Event_length")
  
  pro_marker = c(sur_marker, tran_factor)
  
  seq=c(1:26)
  
  day_levels = c(2,4,6,8,10,11,12,14,16,18,20) # c(0,2,4,6,8,10,11,12,14,16,18,20,22)
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    #K1
    protein1 = bind_rows( g1[marker] %>% mutate(obstime = 2, col_id = 1:nrow(g1)),
                          g2[marker]%>% mutate(obstime = 4, col_id = 1:nrow(g2)),
                          g3[marker]%>% mutate(obstime = 6, col_id = 1:nrow(g3)),
                          g4[marker]%>% mutate(obstime = 8, col_id = 1:nrow(g4)),
                          g5[marker]%>% mutate(obstime = 10, col_id = 1:nrow(g5)),
                          g6[marker]%>% mutate(obstime = 11, col_id = 1:nrow(g6)),
                          g7[marker]%>% mutate(obstime = 12, col_id = 1:nrow(g7)),
                          g8[marker]%>% mutate(obstime = 14, col_id = 1:nrow(g8)),
                          g9[marker]%>% mutate(obstime = 16, col_id = 1:nrow(g9)),
                          g10[marker]%>% mutate(obstime = 18, col_id = 1:nrow(g10)),
                          g11[marker]%>% mutate(obstime = 20, col_id = 1:nrow(g11))) %>%
    
    #K2
    # protein1 = bind_rows( g1[marker] %>% mutate(obstime = 0, col_id = 1:nrow(g1)),
    #                       g2[marker]%>% mutate(obstime = 2, col_id = 1:nrow(g2)), 
    #                       g3[marker]%>% mutate(obstime = 4, col_id = 1:nrow(g3)),
    #                       g4[marker]%>% mutate(obstime = 6, col_id = 1:nrow(g4)),
    #                       g5[marker]%>% mutate(obstime = 8, col_id = 1:nrow(g5)), 
    #                       g6[marker]%>% mutate(obstime = 10, col_id = 1:nrow(g6)),
    #                       g7[marker]%>% mutate(obstime = 11, col_id = 1:nrow(g7)),
    #                       g8[marker]%>% mutate(obstime = 12, col_id = 1:nrow(g8)),
    #                       g9[marker]%>% mutate(obstime = 14, col_id = 1:nrow(g9)),
    #                       g10[marker]%>% mutate(obstime = 16, col_id = 1:nrow(g10)),
    #                       g11[marker]%>% mutate(obstime = 18, col_id = 1:nrow(g11)),
    #                       g12[marker]%>% mutate(obstime = 20, col_id = 1:nrow(g12)),
    #                       g13[marker]%>% mutate(obstime = 22, col_id = 1:nrow(g13))) %>% 
    #   
      mutate(obstime = factor(obstime, levels = day_levels))
    
    protein1 = protein1 %>% mutate(prot_name = pro_marker[i])
    colnames(protein1) = c("Intensity", "obstime", "cell_id", "cell.type")
    prot_list = c(prot_list, list(protein1))
  }   
  return(prot_list)
}

get_pro_data <- function() {
  
  #K1
  pro1 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c01_K1_02_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c02_K1_02_0_0_Day2_Ungated_Ungated.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c03_K1_02_0_0_Day4_Ungated_Ungated.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c04_K1_02_0_0_Day6_Ungated_Ungated.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c05_K1_02_0_0_Day8_Ungated_Ungated.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c06_K1_02_0_0_Day10_Ungated_Ungated.fcs", transformation=FALSE)
  pro7 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c07_K1_02_0_0_Day11_Ungated_Ungated.fcs", transformation=FALSE)
  pro8 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c08_K1_02_0_0_Day12_Ungated_Ungated.fcs", transformation=FALSE)
  pro9 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c09_K1_02_0_0_Day14_Ungated_Ungated.fcs", transformation=FALSE)
  pro10 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c10_K1_02_0_0_Day16_Ungated_Ungated.fcs", transformation=FALSE)
  pro11 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c15_K1_02_0_0_Day18_Ungated_Ungated.fcs", transformation=FALSE)
  pro12 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Kinetic/c16_K1_02_0_0_Day20_Ungated_Ungated.fcs", transformation=FALSE)
  # pro13 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Controls/c19_K1_02_0_0_MNCs_Ungated_Ungated.fcs", transformation=FALSE)
  # pro14 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic1/Ungated/Controls/c17_K1_02_0_0_Jurkat_Ungated_Ungated.fcs", transformation=FALSE)
  # 
  # # #K2
  # pro1 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c01_K2_01_0_0_Day0_Ungated_Ungated.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c02_K2_01_0_0_Day2_Ungated_Ungated.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c03_K2_01_0_0_Day4_Ungated_Ungated.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c04_K2_01_0_0_Day6_Ungated_Ungated.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c05_K2_01_0_0_Day8_Ungated_Ungated.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c06_K2_01_0_0_Day10_Ungated_Ungated.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c07_K2_01_0_0_Day11_Ungated_Ungated.fcs", transformation=FALSE)
  # pro8 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c08_K2_01_0_0_Day12_Ungated_Ungated.fcs", transformation=FALSE)
  # pro9 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c09_K2_01_0_0_Day14_Ungated_Ungated.fcs", transformation=FALSE)
  # pro10 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c10_K2_01_0_0_Day16_Ungated_Ungated.fcs", transformation=FALSE)
  # pro11 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c15_K2_01_0_0_Day18_Ungated_Ungated.fcs", transformation=FALSE)
  # pro12 <- read.FCS("Research/RBC/Data Analysis/normal_kinetics_09_17/Kinetic2/Ungated/Kinetic/c20_K2_01_0_0_Day20_Ungated_Ungated.fcs", transformation=FALSE)
  # 
  #tibble(desc = as.character(pro13@parameters@data$desc), name = names(pro13@parameters@data$desc)) %>% separate(desc, c("Isotope", "Protein"), "_") %>% View()
  
  pro1_df = exprs(pro1)
  pro2_df = exprs(pro2)
  pro3_df = exprs(pro3)
  pro4_df = exprs(pro4)
  pro5_df = exprs(pro5)
  pro6_df = exprs(pro6)
  pro7_df = exprs(pro7)
  pro8_df = exprs(pro8)
  pro9_df = exprs(pro9)
  pro10_df = exprs(pro10)
  pro11_df = exprs(pro11)
  pro12_df = exprs(pro12)
  # pro13_df = exprs(pro13)
  # pro14_df = exprs(pro14)
  
  g1 = data.frame(pro1_df)
  g2 = data.frame(pro2_df)
  g3 = data.frame(pro3_df)
  g4 = data.frame(pro4_df)
  g5 = data.frame(pro5_df)
  g6 = data.frame(pro6_df)
  g7 = data.frame(pro7_df)
  g8 = data.frame(pro8_df)
  g9 = data.frame(pro9_df)
  g10 = data.frame(pro10_df)
  g11 = data.frame(pro11_df)
  g12 = data.frame(pro12_df)
  # g13 = data.frame(pro13_df)
  # g14 = data.frame(pro14_df)
  ###############################################################################
  sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f", 
                 "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4")
  
  other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                  "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "MEF2C")
  
  cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")
  
  
  isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di", "Sm147Di", "Yb171Di",
              
              "Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Sm154Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di", "Nd146Di",
              
              "Nd148Di", "Nd150Di", "Er168Di", "I127Di",
              
              "Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di") 
  
  pro_marker = c(sur_marker, other_marker, tran_factor, cell_cycle, DNA)
  
  seq=c(1:41)
  
  #day_levels = paste0("", c("MNCs", "Jurkats")) #paste0("", c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16, 18 ,20, "MNCs", "Jurkats"))
  day_levels = paste0("", c(0, 2, 4, 6, 8, 10, 11, 12, 14, 16, 18 ,20))
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    protein1 = bind_rows( g1[marker] %>% mutate(Day = "0", col_id = 1:nrow(g1)), 
                          g2[marker]%>% mutate(Day = "2", col_id = 1:nrow(g2)),
                          g3[marker]%>% mutate(Day = "4", col_id = 1:nrow(g3)),
                          g4[marker]%>% mutate(Day = "6", col_id = 1:nrow(g4)),
                          g5[marker]%>% mutate(Day = "8", col_id = 1:nrow(g5)),
                          g6[marker]%>% mutate(Day = "10", col_id = 1:nrow(g6)),
                          g7[marker]%>% mutate(Day = "11", col_id = 1:nrow(g7)),
                          g8[marker]%>% mutate(Day = "12", col_id = 1:nrow(g8)),
                          g9[marker]%>% mutate(Day = "14", col_id = 1:nrow(g9)),
                          g10[marker]%>% mutate(Day = "16", col_id = 1:nrow(g10)), 
                          g11[marker]%>% mutate(Day = "18", col_id = 1:nrow(g11)),
                          g12[marker]%>% mutate(Day = "20", col_id = 1:nrow(g12))) %>% 
                          # g13[marker]%>% mutate(Day = "MNCs", col_id = 1:nrow(g13)),
                          # g14[marker]%>% mutate(Day = "Jurkats", col_id = 1:nrow(g14))) %>%
      mutate(Day = factor(Day, levels = day_levels))
    
    protein1 = protein1 %>% mutate(prot_name = pro_marker[i])
    colnames(protein1) = c("Intensity", "obstime", "cell_id", "cell.type")
    prot_list = c(prot_list, list(protein1))
    
  }  
  
  return(prot_list)
}


all_data = as_tibble(bind_rows(get_pro_data()))

all_data = all_data %>% mutate(cell = factor(paste(cell_id, obstime)))
all_data = all_data %>% mutate(capture = obstime)

express_data0 = all_data %>% spread(cell.type, Intensity)
#####################################################################
sur_marker = c("CD34", "CD71", "CD45RA", "CD123", "CD49f", 
               "CD90", "CD44", "CD41", "CD235ab", "CD33", "CXCR4",
               "CD36", "CD38")

other_marker = c("HBB", "HBA", "H3")

tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "TAL1", "GATA2",
                "RUNX1", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B", "Flag_Tag",
                "FLI1", "NFE2p45")

cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")

DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194")

gated_data = gene_data %>%  #express_data0 %>%
  mutate(Day = obstime) %>% 
  select(-cell_id, -cell, -capture, -obstime, -cell_cycle, -HBB, -H3, -DNA, -mean_Ir, -mean_Pt) %>% 
  mutate_if(is.numeric, ~ floor(.x)) %>% 
  filter((!Day %in% c("MNCs", "Jurkats")))

prot_name = gated_data %>% select(-Day) %>% colnames() %>% as.character()

prot_list = data.frame(prot_name = prot_name, prot_index = 1:length(prot_name)) %>% as.tibble() %>% 
  mutate(prot_name = as.character(prot_name))

gated_data %>% group_by(Day) %>% summarise(count=n())

cells_per_day = 3000
abv_max = gated_data %>% group_by(Day) %>% filter(n()<=cells_per_day) %>% ungroup()
sub_data = gated_data %>% group_by(Day) %>% filter(n()>cells_per_day) %>% sample_n(cells_per_day) %>% ungroup() %>%
  bind_rows(abv_max) %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% 
  arrange(Day) 

x_data   = sub_data %>% select(-Day) %>% as.matrix()
vec_d    = sub_data %>% mutate(new_indx = as.numeric(factor(as.character(Day), levels = unique(as.numeric(Day))))) %>% pull(new_indx)
indx_tp  = sub_data %>% group_by(Day) %>% summarise(num_cells = n()) %>% mutate(num_cells  = cumsum(num_cells)) %>% 
  mutate(from = c(0, num_cells[-n()]) + 1) %>% select(from, num_cells) %>% as.matrix()

x_indx = x_data %>% as_tibble() %>% mutate_all(~as.numeric(as.factor(.x)))
x_unique_aux = x_data %>% as_tibble() %>% map2(1:ncol(.), ~ as.factor(.x) %>% levels() %>% 
                                                 tibble(unique_count = ., prot =.y)) %>% bind_rows()
unique_prot_x = x_unique_aux %>% mutate(unique_count = as.numeric(unique_count)) %>% pull(unique_count)
indx_prot  = x_unique_aux %>% group_by(prot) %>% summarise(num_vals = n()) %>% mutate(num_vals  = cumsum(num_vals)) %>% 
  mutate(from = c(0, num_vals[-n()]) + 1) %>% select(from, num_vals) %>% as.matrix()

data_cyto = list(c = 25, k = ncol(x_data), N = nrow(x_data), x = x_data, d = max(vec_d), 
                 day_index = vec_d, indx_tp = indx_tp,
                 num_unique_counts = length(unique_prot_x),
                 unique_prot_x = unique_prot_x,
                 indx_prot = indx_prot,
                 x_indx = x_indx)

vec_d %>% unique()

library(rstan)

options(mc.cores = parallel::detectCores())
doMC::registerDoMC(cores = 14)
rstan_options(auto_write = TRUE)

model   = stan_model("~/mix_nb_model_Ed_vec_V4.stan") #V4
model   = stan_model("~/Software/cluster_model/mix_nb_model_Ed_vec_V4.stan") 

#vb a function in stan model
runs = c(1:4)

library(doParallel)

foreach(number = runs)%dopar%{
  start_time <- Sys.time()
  map_est4   = optimizing(model, data = data_cyto, iter = 3000, as_vector = F)
  #map_est4   = vb(model, data = data_cyto, algorithm = "fullrank")
  #print(start_time - Sys.time())
  
  direct = paste("Results_by_date/Dec2018/Oct/K2/", "run_", number, ".rda", sep="") #Desktop/Sep2018/Controls_Oct/Oct_K1/", "K1_mncs_nn_"
  save(map_est4, file = direct)
  
  return(start_time - Sys.time())
}

save(sub_data, file = "Results_by_date/Dec2018/Oct/K3/sampled_cell_data.rda")

###############################################################################
load(file = "Results_by_date/Dec2018/May/K3/run_1.rda")
#############if data from elsewhere###############
c = 25

day_vals = sub_data$Day %>% unique()
day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)

run_fit = map_est4

cluster_plot <- function(run_fit){
  
  vals = run_fit$par 
  
  time = tibble(Day = day_index_0$day_index,  real_time= day_index_0$Day) 
  
  cluster = vals$theta %>% as.tibble() %>% #t()
    pull(value) %>% as.tibble() %>%
    rename(vals = value) %>%
    mutate(cluster = 1:c)
  
  # cluster = vals$theta %>% t() %>% as.tibble() %>%
  # rename(vals = V1) %>%
  # mutate(cluster = 1:c)
  
  phi = vals$phi_s_cp %>% as.tibble() %>% 
    mutate(cluster = 1:c) %>% 
    gather(protein, std, -cluster) %>%
    mutate(protein = as.double(str_remove_all(protein, "V")))
  
  signal = vals$mu_s_cp %>% as.tibble() %>% 
    mutate(cluster = 1:c) %>% 
    gather(protein, vals, -cluster) %>%
    mutate(protein = as.numeric(str_remove_all(protein, "V"))) %>% 
    left_join(prot_list, by = c("protein" = "prot_index")) %>% 
    left_join(phi) 
  
  signal = signal %>% 
    group_by(prot_name) %>% 
    mutate(max_vals = max(vals), min_vals = min(vals), mid_vals = mean(vals), sd_vals = sd(vals)) %>% 
    ungroup()
  
  day_vals = gated_data$Day %>% unique()
  day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)
  
  mu_c_fit = vals$mu_s_cp %>% t()
  
  theta0_c_fit = vals$theta0_s_cp %>% t()
  
  phi_c_fit = vals$phi_s_cp %>% t()
  
  mu_bg_fit = vals$mu_bg %>% as.vector()
  
  st_fit = vals$scaling_t
  st_fit = c(1, st_fit)
  
  theta_fit = rep(cluster$vals, length(day_vals)) %>% 
    matrix(nrow = c)
  
  data_x = gated_data %>%
    as.tibble() %>% 
    mutate(Day = as.numeric(as.character(Day))) %>% 
    left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 
  
  data_x1 = data_x %>% select(-Day, -day_index) %>% as.matrix() #Ir191, Ir193, IdU, p_Rb
  
  day_indx = data_x %>%  pull(day_index)
  
  all_cells = predict_cluster(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx) # 
  
  #print(table(all_cells)) %>% sum()
  cluster_prob = all_cells %>% as.tibble()
  
  full_clusters = all_cells %>% as.tibble() %>% unique() %>% pull(V1)
  
  p1 = signal %>% mutate(valss = (vals-mid_vals)/(sd_vals)) %>% 
    #filter(cluster %in% full_clusters) %>% 
    #filter(prot_name %in% key_proteins) %>% 
    ggplot(aes(x = prot_name, y= vals)) + 
    geom_point(aes(colour = prot_name)) +
    scale_y_log10()+
    geom_text(aes(label = prot_name), size=3.0)+
    facet_wrap(~cluster, nrow = 2) +
    theme(legend.position = "none")
  
  p2 = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
    mutate(prob = num/sum(num)) %>% ungroup() %>% 
    select(-num) %>% 
    spread(cluster, prob, fill = 0) %>% 
    gather(cluster, prob, -Day) %>% 
    mutate(cluster = as.numeric(cluster)) %>% 
    ggplot(aes(x = Day, y = prob))+
    geom_point()+
    facet_wrap(~cluster, nrow=2)
  
  gridExtra::grid.arrange(p1, p2, ncol=1)
}

cluster_plot(map_est4)

############################################################
############################
all_cells = predict_cluster(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx)
all_cells_prob = predict_cluster_prob(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx)

#find the empty clusters
low_cluster = table(all_cells) %>% as.tibble() %>% rename(cluster = all_cells) %>% 
  mutate(cluster = as.numeric(cluster)) %>% filter(n<=200) %>% pull(cluster)

#find the missing clusters
missing_cluster = base::setdiff(1:c, unique(all_cells) %>% as.tibble() %>% pull(V1))

empty_cluster = c(low_cluster, missing_cluster)

cluster_prob = all_cells_prob %>% as.matrix()

cluster_prob = cluster_prob[, -(empty_cluster)]

colnames(cluster_prob) <- paste0("p", 1:(c-length(empty_cluster)))

old_cluster = 1:c

old_cluster = old_cluster[-empty_cluster]

new_cluster = data.frame(cluster = old_cluster, new_cluster = 1:length(old_cluster)) %>% as.tibble()

order_cluster = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>%
  group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  group_by(cluster) %>%
  mutate(x = max(prob)) %>%
  filter(prob == x) %>%
  ungroup() %>% select(Day, cluster) %>% arrange(., Day) %>%
  select(cluster) %>% 
  filter(cluster %in% old_cluster) %>% 
  mutate(cluster = as.integer(cluster)) %>% 
  left_join(new_cluster) %>% 
  unique() %>% as.data.frame() %>% pull(cluster) %>% 
  as.integer()


###########remove empty clusters#######
vals = run_fit$par 

time = tibble(Day = day_index_0$day_index,  real_time= day_index_0$Day) 

cluster = vals$theta %>% as.tibble() %>%
  pull(value) %>% as.tibble() %>%
  rename(vals = value) %>%
  mutate(cluster = 1:c)

phi = vals$phi_s_cp %>% as.tibble() %>% 
  mutate(cluster = 1:c) %>% 
  gather(protein, std, -cluster) %>%
  mutate(protein = as.double(str_remove_all(protein, "V")))

signal = vals$mu_s_cp %>% as.tibble() %>% 
  mutate(cluster = 1:c) %>% 
  gather(protein, vals, -cluster) %>%
  mutate(protein = as.numeric(str_remove_all(protein, "V"))) %>% 
  left_join(prot_list, by = c("protein" = "prot_index")) %>% 
  left_join(phi) 

signal = signal %>% 
  group_by(prot_name) %>% 
  mutate(max_vals = max(vals), min_vals = min(vals), mid_vals = mean(vals), sd_vals = sd(vals)) %>% 
  ungroup()

day_vals = gated_data$Day %>% unique()
day_index_0 = data.frame(day_index = 1:length(day_vals), Day = day_vals)

mu_c_fit = vals$mu_s_cp %>% t()

theta0_c_fit = vals$theta0_s_cp %>% t()

phi_c_fit = vals$phi_s_cp %>% t()

mu_bg_fit = vals$mu_bg %>% as.vector()

st_fit = vals$scaling_t
st_fit = c(1, st_fit)

theta_fit = rep(cluster$vals, length(day_vals)) %>% 
  matrix(nrow = c)

data_x = gated_data %>%
  as.tibble() %>% 
  mutate(Day = as.numeric(as.character(Day))) %>% 
  left_join(day_index_0 %>% as.tibble() %>% mutate(Day= as.numeric(as.character(Day))), by = c("Day" = "Day")) 

data_x1 = data_x %>% select(bCatenin:Thy1) %>% as.matrix()

day_indx = data_x %>%  pull(day_index)

all_cells = predict_cluster(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx)
all_cells_prob = predict_cluster_prob(data_x1, theta_fit, theta0_c_fit, mu_c_fit, phi_c_fit, mu_bg_fit, st_fit, day_indx)

print(table(all_cells)) %>% sum()

#find the empty clusters
low_cluster = table(all_cells) %>% as.tibble() %>% rename(cluster = all_cells) %>% 
  mutate(cluster = as.numeric(cluster)) %>% filter(n<=200) %>% pull(cluster)

#find the missing clusters
missing_cluster = base::setdiff(1:c, unique(all_cells) %>% as.tibble() %>% pull(V1))

empty_cluster = c(low_cluster, missing_cluster)

cluster_prob = all_cells_prob %>% as.matrix()

cluster_prob = cluster_prob[, -(empty_cluster)]

colnames(cluster_prob) <- paste0("p", 1:(c-length(empty_cluster)))

old_cluster = 1:c

old_cluster = old_cluster[-empty_cluster]

new_cluster = data.frame(cluster = old_cluster, new_cluster = 1:length(old_cluster)) %>% as.tibble()

order_cluster = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>%
  group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  group_by(cluster) %>%
  mutate(x = max(prob)) %>%
  filter(prob == x) %>%
  ungroup() %>% select(Day, cluster) %>% arrange(., Day) %>%
  select(cluster) %>% 
  filter(cluster %in% old_cluster) %>% 
  mutate(cluster = as.integer(cluster)) %>% 
  left_join(new_cluster) %>% 
  unique() %>% as.data.frame() %>% pull(cluster) %>% 
  as.integer()

ordered_cluster = tibble(cluster = as.integer(order_cluster), state = 1:(c-length(empty_cluster)))

#print(table(all_cells)) %>% sum()
key_proteins = c("DNA-1", "Oct4", "Klf4", "Sox2", "Ki67", "Nanog", "CD54", "CD73", "CD24", "CD140a")

p1 = signal %>% mutate(valss = (vals-mid_vals)/(sd_vals)) %>% 
  filter(cluster %in% old_cluster) %>% 
  left_join(new_cluster) %>% 
  left_join(ordered_cluster) %>% 
  #filter(prot_name %in% key_proteins) %>% 
  ggplot(aes(x = prot_name, y= vals)) + 
  geom_point(aes(colour = prot_name)) +
  scale_y_log10()+
  geom_text(aes(label = prot_name), size=3.0)+
  facet_wrap(~state, nrow = 2) +
  theme(legend.position = "none")

p2 = data_x %>% mutate(cluster = all_cells) %>% group_by(Day, cluster) %>% summarise(num= n()) %>% group_by(Day) %>% 
  mutate(prob = num/sum(num)) %>% ungroup() %>% 
  select(-num) %>% 
  spread(cluster, prob, fill = 0) %>% 
  gather(cluster, prob, -Day) %>% 
  mutate(cluster = as.numeric(cluster)) %>% 
  filter(cluster %in% old_cluster) %>% 
  left_join(new_cluster) %>% 
  left_join(ordered_cluster) %>% 
  ggplot(aes(x = Day, y = prob))+
  geom_point()+
  facet_wrap(~state, nrow=2)

gridExtra::grid.arrange(p1, p2, ncol=1)
############################################################
express_data0 = data_x %>% mutate(cluster = all_cells) %>%
  mutate(cluster = as.numeric(cluster)) %>% 
  filter(cluster %in% old_cluster) %>% 
  left_join(new_cluster) %>% 
  left_join(ordered_cluster) 
##############SPRING plot data generation###################
express_data1 = express_data0 %>% filter(Day != 0, Day != 18, Day != 20, Day != 22) %>% group_by(Day) %>%
  sample_n(2500) %>% ungroup()

express_data2 = express_data0 %>% filter(Day == 18 | Day == 20 | Day == 22)

express_data3 = express_data0 %>% filter(Day == 0) %>% filter(CD34 >= sinh(4.5))

express_data = rbind(express_data1, express_data2, express_data3)

express_data$Day <- NULL
express_data$day_index <- NULL
express_data$cluster <- NULL
express_data$new_cluster <- NULL
express_data$state <- NULL

x1=express_data %>% asinh()
write.table(x1, file='Software/SPRING_new_version_master/datasets/May_K3/prot_data.csv',quote=FALSE, sep='\t',col.names = FALSE, row.names = FALSE)

#cell groupings according to days
day_data = rbind(express_data1, express_data2, express_data3)
y <- day_data$Day
z <- day_data$state

x2 = data.frame(obstime = y, cluster = z)
x2 = t(x2)
rownames(x2) <- c('Day', 'Cluster') 
write.table(x2, file='Software/SPRING_new_version_master/datasets/May_K3/classification.csv',quote=FALSE, sep=',',col.names = FALSE)

all_prot = names(express_data)
write(all_prot, file='Software/SPRING_new_version_master/datasets/May_K3/prot_name.csv', ncolumns = length(all_prot), sep='\t')

flag_data = data_x %>% select(Flag_Tag, CD41, CD235ab) %>% asinh()
write.table(flag_data, file='Software/SPRING_new_version_master/datasets/Flag_FLI1/flag.csv', sep='\t', col.names = TRUE, row.names = FALSE)

