library(flowCore)
library(tidyverse)
library(stringr)

#tibble(desc = as.character(pro1@parameters@data$desc), name = names( pro1@parameters@data$desc)) %>% separate(desc, c("Isotope", "Protein"), "_") %>% View()

get_pro_data <- function() {
  
  # pro1 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c01_sample-1_01_0_Day2_live.fcs", transformation=FALSE)
  # pro2 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c02_sample-1_01_0_Day4_live.fcs", transformation=FALSE)
  # pro3 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c03_sample-1_01_0_Day6_live.fcs", transformation=FALSE)
  # pro4 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c04_sample-1_01_0_Day8_live.fcs", transformation=FALSE)
  # pro5 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c05_sample-1_01_0_Day10_live.fcs", transformation=FALSE)
  # pro6 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c16_sample-1_01_0_Day12_live.fcs", transformation=FALSE)
  # pro7 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c17_sample-1_01_0_Day14_live.fcs", transformation=FALSE)
  # pro8 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c18_sample-1_01_0_Day16_live.fcs", transformation=FALSE)
  # pro9 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c19_sample-1_01_0_Day18_live.fcs", transformation=FALSE)
  # pro10 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c20_sample-1_01_0_Day20_live.fcs", transformation=FALSE)
  
  pro1 <- read.FCS("Research/T3-titration3-15-5-2017/raw-data-NOT-gated-on-live-cells/c12_d1-1_01_0_MNCs_1:1_MNCs_Ungated.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/T3-titration3-15-5-2017/raw-data-NOT-gated-on-live-cells/c12_d1-2_01_0_MNCs_1:2_MNCs_Ungated.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/T3-titration3-15-5-2017/raw-data-NOT-gated-on-live-cells/c12_d1-4_01_0_MNCs_1:4_MNCs_Ungated.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/T3-titration3-15-5-2017/raw-data-NOT-gated-on-live-cells/c13_d1-1_01_0_Day 0_1:1_Day 0_Ungated.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/T3-titration3-15-5-2017/raw-data-NOT-gated-on-live-cells/c13_d1-2_01_0_Day 0_1:2_Day 0_Ungated.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/T3-titration3-15-5-2017/raw-data-NOT-gated-on-live-cells/c13_d1-4_01_0_Day 0_1:4_Day 0_Ungated.fcs", transformation=FALSE)

  pro1_df = exprs(pro1)
  pro2_df = exprs(pro2)
  pro3_df = exprs(pro3)
  pro4_df = exprs(pro4)
  pro5_df = exprs(pro5)
  pro6_df = exprs(pro6)
  # pro7_df = exprs(pro7)
  # pro8_df = exprs(pro8)
  # pro9_df = exprs(pro9)
  # pro10_df = exprs(pro10)
  # pro11_df = exprs(pro11)
  # pro12_df = exprs(pro12)
  # pro13_df = exprs(pro13)
  # pro14_df = exprs(pro14)
  
  g1 = data.frame(pro1_df) 
  g2 = data.frame(pro2_df) 
  g3 = data.frame(pro3_df) 
  g4 = data.frame(pro4_df) 
  g5 = data.frame(pro5_df) 
  g6 = data.frame(pro6_df) 
  # g7 = data.frame(pro7_df) 
  # g8 = data.frame(pro8_df) 
  # g9 = data.frame(pro9_df) 
  # g10 = data.frame(pro10_df) 
  # g11 = data.frame(pro11_df)
  # g12 = data.frame(pro12_df)
  # g13 = data.frame(pro13_df)
  # g14 = data.frame(pro14_df)
  ###############################################################################
  sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f",
                 "CD90", "CD44", "CD41", "CD235ab")
  
  #other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                  "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B")
  
  #cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194") #, "Event_Length")
  
  isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Yb174Di",
              
              #"Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Sm154Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di",
  
  #"Nd148Di", "Nd150Di", "Er168Di", "I127Di",
  
              "Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di") #, "Event_length")
  
  pro_marker = c(sur_marker, tran_factor, DNA)
  
  seq=c(1:31)
  
 day_levels = paste0("", c("MNCs1", "MNCs2", "MNCs3", "Day0_1", "Day0_2", "Day0_3"))
 
  prot_list = list()

  for(i in seq) {
  
    marker = isotops[i]
    
    protein1 = bind_rows( g1[marker] %>% mutate(Day = "MNCs1", col_id = 1:nrow(g1)), 
                          g2[marker]%>% mutate(Day = "MNCs2", col_id = 1:nrow(g2)), 
                          g3[marker]%>% mutate(Day = "MNCs3", col_id = 1:nrow(g3)), 
                          g4[marker]%>% mutate(Day = "Day0_1", col_id = 1:nrow(g4)), 
                          g5[marker]%>% mutate(Day = "Day0_2", col_id = 1:nrow(g5)), 
                          g6[marker]%>% mutate(Day = "Day0_3", col_id = 1:nrow(g6))) %>% 
                          # g7[marker]%>% mutate(Day = "11", col_id = 1:nrow(g7)),  
                          # g8[marker]%>% mutate(Day = "12", col_id = 1:nrow(g8)),  
                          # g9[marker]%>% mutate(Day = "14", col_id = 1:nrow(g9)),
                          # g10[marker]%>% mutate(Day = "16", col_id = 1:nrow(g10)), #) %>% 
                          # g11[marker]%>% mutate(Day = "18", col_id = 1:nrow(g11)),
                          # g12[marker]%>% mutate(Day = "20", col_id = 1:nrow(g12)),
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
#all_data = all_data %>% mutate(Intensity = asinh(Intensity))
#all_data = all_data %>% mutate(Intensity = hill_asinh(1e6, Intensity, 5))

hill_asinh <- function(Kd=1e6, L, n =5){

  y = L^n / (Kd + L^n)
  
  vals = y*asinh(L/150)
  
}

all_data = all_data %>% mutate(Intensity = hill_asinh(1e6, Intensity, 5))
#all_data %>% group_by(Day, Cell_id) %>% summarise(num_prot = n(), sum_prot = sum(Intensity), num_above = sum(Intensity > 2500), na_num = sum(is.na(Intensity)))
##################################################################################################

#stat_prot <- all_data %>% group_by(Day, Cell_id) %>% summarise(num_prot = n(), num_above = sum(Intensity > 2500), num_below = sum(Intensity < 1e-2))

# stat_prot <- all_data %>% group_by(Day, Cell_id) %>% summarise(num_below = sum(Intensity < 1e-2))

stat_prot %>% ggplot() +
  geom_line(mapping = aes(x=Cell_id, y = num_above)) +
  facet_wrap(~Day, ncol = 2)

stat_prot %>% ggplot() +
  geom_line(mapping = aes(x=Cell_id, y = num_below)) +
  facet_wrap(~Day, ncol = 2)

pic=list()
pic[[1]] = ggplot(all_data, aes(x = Day) )+ 
  geom_histogram(stat="count", bins = 100) +
  facet_wrap(~ Day, ncol = 2)

typeof(pic)

ggplot(all_data, aes(x = Day)) + 
  geom_bar() +
  facet_wrap(~Day, ncol = 2)

ggplot(data = all_data, aes(x = all_data$Day)) + 
  geom_histogram(stat="count") #+

####################################################################
prot_levels = c("CD41", "IdU", "Ce140", "CD235ab", "CD45RA", "HBB", "CEBPa", "MEF2C", "CD33", "p_HH3", "CD34", 
                "p_Rb", "CD123", "MAFG", "CD44", "NFE2p45", "CD36", "GATA1", "RUNX1", "IKZF1", "ATRX", "CD90", 
                "TAL1", "GATA2", "CD49f", "KLF1", "KAT3B", "PU1", "cyclin_B1", "FLI1", "HBA", "CXCR4", "CD38", 
                "c_Jun", "H3", "CD71", "c_Myc", "Ir191", "Ir193", "Pt195", "Pt194")

express_data = all_data %>% spread(prot_name, Intensity) %>% order_by(prot_levels)

gated_data = express_data %>% filter(Ce140 <= 5.0)

data_gating = function(Intensity){
  
  return(Intensity >= 3.75)
}

compare_data = gated_data %>% mutate(Pt_hi = data_gating(Pt195)) 

surface_data = compare_data %>% select(Day:Cell_id, CD123:CD90, CXCR4, Pt194, Pt_hi) %>% gather(prot_name, Intensity, CD123:Pt194)

Gb_Cc_DNA = compare_data %>% select(Day:Cell_id, H3:HBB, IdU, p_HH3, p_Rb, cyclin_B1, Ir191, Ir193, Pt194, Pt_hi) %>% gather(prot_name, Intensity, H3:Pt194)

Tf_data = compare_data %>% select(Day:Cell_id, CEBPa:c_Myc, ATRX, FLI1:GATA2, IKZF1, KAT3B:NFE2p45, PU1:TAL1, Pt194, Pt_hi) %>% gather(prot_name, Intensity, CEBPa:Pt194)

figure1 = gated_data %>% ggplot(aes(Ir191, Pt195, col = HBA)) + #fill = Pt_hi
  #geom_hex(bins = 80) +
  geom_point(alpha = 0.1)+
  facet_wrap(~Day,  ncol = 2)
  #alpha = 0.1
  #geom_density_2d(aes(colour = "red"))
  #stat_density_2d(aes(fill = ..level..), , geom = "polygon") +
png(filename = "Desktop/Oct24/test.png", width = 790, height = 780, units = "px", pointsize = 12, bg = "white")
plot(figure1)
dev.off()

compare_data %>% group_by(Day) %>% summarise(low = sum(Pt_hi == FALSE), high = sum(Pt_hi == TRUE))
                                             
compare_data %>% ggplot(aes(x=Pt195)) +
  geom_histogram(bins = 100) +
  #ylim(0,0.5)+
  facet_wrap(~Day)

# + #, scales = "free") 
  labs(x = "Lu175", y = "Ir191")
####################################################################
mean_irid = function(ir1, ir2){
    reg = lm(ir1~ir2)
    # reg = lm(data = gated_data, Ir191~Ir193)
    (predict(reg) + ir1)/2
  }
  
express_data = all_data %>% spread(prot_name, Intensity)

express_data1 = express_data %>% filter(Ce140 <= sinh(5.0)) %>% #, Day =="Jurkats" | Day == 11 | Day == 0 | Day =="MNCs") 
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) %>% 
  filter(mean_Pt >= sinh(5.0))

  
cell_data = express_data %>% filter(Ce140 <= sinh(5.0)) %>% #, Day =="Jurkats" | Day == 11 | Day == 0 | Day =="MNCs") 
 mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193))

gated_data1 = cell_data %>% filter(mean_Pt <= sinh(5.0) | CD235ab >= sinh(6.5))

gated_data = gated_data1 %>% filter(mean_Ir >= sinh(5.0) | CD235ab >= sinh(6.5), 
                                    mean_Ir <= sinh(8.0) | CD235ab >= sinh(6.5))

tukey_lim = function(x, k = 2){
  q1_q3 = quantile(x, c(0.75, 0.9))
  upper_lim = q1_q3[2] + k*(q1_q3[2]-q1_q3[1])
  upper_lim
}

gated_data = express_data %>% filter(FLI1 < tukey_lim(FLI1, 3) & KLF1 < tukey_lim(KLF1, 3))
gated_data = express_data %>% filter(GATA1 < tukey_lim(GATA1, 5) & PU1 < tukey_lim(PU1,15))

seq=c(1:42)

for(i in seq) {
  
marker = pro_marker[i]

figure = express_data %>% #filter(Day == 20) %>% 
  ggplot(aes(asinh(CD34), asinh(CD235ab))) +
  geom_point(alpha = 0.1)+
  geom_density_2d(aes(colour = "red"))+
  #ylim(0,1000)+
  #xlim(0,1000)+
  facet_wrap(~Day,  ncol = 4)

filename = paste("Desktop/Jan2018/PU1_GATA1_over.png")
png(filename, width = 1500, height = 1200, units = "px", pointsize = 20, bg = "white")
plot(figure)
dev.off()

}

gated_data %>% select(-Ce140, -Ir191, -Ir193, -Pt194, -Pt195) %>% gather(prot_name, Intensity, ATRX:TAL1) %>% mutate(Intensity = asinh(Intensity)) %>% 
  ggplot(aes(x = Day, y = Intensity, fill = Day)) +
  geom_violin(size=0.2, scale = "width", bw = 0.1)  +
  facet_wrap(~prot_name, ncol = 5, scales = "free")


express_data %>% ggplot(aes(x=CD235ab)) +
  geom_histogram(bins = 10) +
  ylim(0,350) +
  facet_wrap(~Day, scales = "free_y")

correct_data = all_data %>% spread(prot_name, Intensity) %>% filter(Ce140 <= hill_asinh(sinh(5.0))) #, Day =="Jurkats" | Day =="MNCs") #Day == 11 | Day == 0 |

select_protein = correct_data %>% select(ATRX:CD90, CEBPa:cyclin_B1, FLI1:IKZF1, KAT3B:p_Rb, PU1:TAL1)

sum_prot_data = correct_data %>% mutate(sum_prot = rowSums(select_protein)) %>% mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193)) 

gated_data1 = sum_prot_data %>% filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

gated_data = gated_data1 %>% filter(mean_Ir >= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5), mean_Ir <= sinh(8.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))

save(gated_data, file="Desktop/Nov9/Data_K1.rda")

#gather_data = gated_data %>% gather(prot_name, Intensity, ATRX:mean_Ir)
#gather_data = gather_data %>% mutate(Intensity = asinh(Intensity))
G0_data = gated_data %>% group_by(Day) %>% summarise(G0_num = sum(p_Rb <= sinh(5.0))/n()) 
G1_data = gated_data %>% group_by(Day) %>% summarise(G1_num = sum(p_Rb > sinh(5.0) & IdU < sinh(3.5) & mean_Ir <= sinh(6.15))/n())
S_data = gated_data %>% group_by(Day) %>% summarise(S_num = sum(p_Rb > sinh(5.0) & IdU >= sinh(3.5))/n())
G2_M_data = gated_data %>% group_by(Day) %>%  summarise(G2_M_num = sum(IdU < sinh(3.5) & mean_Ir > sinh(6.15) & p_Rb > sinh(5.0))/n())

cycle_data = bind_cols(G0_data, G1_num = G1_data$G1_num, S_num = S_data$S_num, G2_M_num = G2_M_data$G2_M_num) %>% gather(dates, percentage, G0_num:G2_M_num)

cycle_data1 = cycle_data %>% select(G1_num:G2_M_num) 

cycle_data = cycle_data %>% mutate(sum_cycle = rowSums(cycle_data1), G1_cycle = G1_num/sum_cycle, G2_M_cycle = G2_M_num/sum_cycle, S_cycle = S_num/sum_cycle) %>% 
  select(Day,G1_cycle:S_cycle) %>% gather(dates, percentage, G1_cycle:S_cycle)

figure1 = express_data %>% filter(Cell_id < 10000) %>% 
  ggplot(aes(x=Ir191, y=IdU)) +# +, colour = dates, group=dates)) +
  geom_point(size=2.5, alpha=0.5) +
  #geom_line(linetype=2)

#  lines(x=Day, y=percentage)
  facet_wrap(~Day,  nrow = 4, scale = "free")
#alpha = 0.1
#stat_density_2d(aes(fill = ..level..), , geom = "polygon") +
png(filename = "Desktop/11.png", width = 1000, height = 600, units = "px", pointsize = 25, bg = "white")
plot(figure1)
dev.off()

gated_data %>% select(contains("Pt")) %>% 
  View()
#****************************#
get_max_dens = function(x){
  d = density(x, n=1e6)
  i = which.max(d$y)
  d$x[i]
}

gated_data = sum_prot_data %>% mutate(mean_Pt = (Pt194 + Pt195)/2, mean_Ir = mean_irid(Ir191, Ir193)) 

hlin_vals = gated_data %>% group_by(Day) %>% summarise(peak_val =get_max_dens(mean_Ir), double_val = 2*peak_val)
#****************************#

figure1 = gated_data %>% #filter(Day == "Jurkats") %>% 
  
  ggplot(aes(mean_Ir, sum_prot))+ # col = Ir193
    #geom_hex(bins = 80) +
    geom_point(alpha = 0.1)+
    geom_density_2d(colour = "red")+
    xlim(0,1250)+
    ylim(0,5000)+
    #scale_y_log10() + 
    #scale_x_log10()+
    geom_vline(data =hlin_vals, aes(xintercept = peak_val)) + 
    geom_vline(data= hlin_vals, aes(xintercept = double_val)) + 
    facet_wrap(~Day,  nrow = 2)
  #alpha = 0.1
  #stat_density_2d(aes(fill = ..level..), , geom = "polygon") +
png(filename = "Desktop/Oct30/test.png", width = 1000, height = 1000, units = "px", pointsize = 25, bg = "white")
plot(figure1)
dev.off()
  
gated_data %>% group_by(Day) %>% summarise(low = sum(Ir191 <=6.25)/n(), high = sum(Ir191 >6.25)/n())

data_trans = gated_data %>% group_by(Day) %>% mutate(mean_Ir = mean_irid(Ir191, Ir193)) %>% ungroup()

#*******************************#
seq=c(1:41)

for(i in seq) {
  
  marker = prot_levels[i]
  
  figure = gated_data  %>% 
    ggplot(aes(x=gated_data[marker])) + #ceiling(
    geom_histogram(bins = 200) +
    #xlim(1, 250)+
    #scale_x_log10()+
    facet_wrap(~Day,  ncol = 5, scales = "free")
  
  filename = paste("Desktop/Nov22/hill_transform/hist_", marker, ".png")
  png(filename, width = 1500, height = 1200, units = "px", pointsize = 20, bg = "white")
  plot(figure)
  dev.off()
  
}

gated_data %>% filter(Day == "Jurkats") %>% ggplot(aes(x=mean_Ir))+
    geom_histogram(bins = 200) +
  #geom_vline(data= hlin_vals, aes(xintercept = peak_val)) + 
  #geom_vline(data= hlin_vals, aes(xintercept = double_val)) + 
  #ylim(0,250)+
  #xlim(0, 2000)+
  facet_wrap(~Day) # + #, scales = "free") labs(x = "Lu175", y = "Ir191")

get_max_dens = function(x){
  d = density(x, n=1e6)
  i = which.max(d$y)
  d$x[i]
}

gated_data%>% filter(IdU < sinh(3.5)) %>% group_by(Day) %>% summarise(peak_val =get_max_dens(Ir191), double_val = 2*peak_val)

gated_data %>% filter(IdU < sinh(3.5), Ir191 < 750) %>% mutate(mean_Irid = (Ir191 + Ir193)/2)%>% ggplot(aes(x=mean_Irid)) +
  +     geom_histogram(bins = 100) +
  +     xlim(0,1000)+
  +   geom_vline(data= hlin_vals, aes(xintercept = peak_val)) + 
  +   geom_vline(data= hlin_vals, aes(xintercept = double_val)) + 
  +   #ylim(0,0.5)+
  +     facet_wrap(~Day) # + #, scales = "free") labs(x = "Lu175", y = "Ir191")

data_trans %>% filter(IdU < sinh(3.5), mean_Ir < 750) %>% ggplot(aes(x=Ir191)) +
  geom_histogram(bins = 200) +
  geom_vline(data= hlin_vals, aes(xintercept = peak_val)) + 
  geom_vline(data= hlin_vals, aes(xintercept = double_val)) + 
  #ylim(0,0.5)+
  facet_wrap(~Day) # + #, scales = "free") labs(x = "Lu175", y = "Ir191")

value_5_95 <- function(x) { 
  
  value = quantile(x, probs = c(0.05, 0.95))
  
  return(value)}

value_5_95(c(1,2,3,4,5))
############################################################
gated_data_k1 %>% select(Day, Cell_id, CD235ab) %>% mutate(CD235ab = asinh(CD235ab)) %>% ggplot(aes(x = Day, y = CD235ab, fill = Day)) +
  geom_violin(size=0.1, scale = "width", bw = 0.1) 

new_gated_data %>% mutate(obstime = factor(obstime, levels = c(0,2,4,6,8,10,11,12,14,16,18,20))) %>% #mutate(obstime = as.character(obstime)) %>%  
  ggplot(aes(x = obstime, y = Intensity, group = obstime, fill=obstime)) +
  geom_violin(size=0.1, scale = "width", bw = 0.1)  +
  facet_wrap(~prot_name, ncol = 4, scales = "free") 

pic = gated_data %>% gather(prot_name, Intensity, CD235ab) %>% mutate(Intensity = asinh(Intensity)) %>% 
  filter(Day == "2"|Day == "6"| Day == "12"| Day == "14"| Day =="20") %>% 
  ggplot(aes(x = Day, y = Intensity, fill = Day)) +
  geom_violin(size=0.1, scale = "width", bw = 0.05)  +
  #scale_y_log10() +  
  #stat_summary(fun.y=mean, geom="point", size=2, colour="red") +
  #geom_boxplot(width=0.1) +
  facet_wrap(~prot_name, ncol = 2, scales = "free") #ATRX, CD235ab, CD34, HBA, KLF1, c_Jun

png("Desktop/Nov22/hill_transform/vio.png",width = 1400, height = 1300, units = "px", pointsize = 12, bg = "white")
plot(pic)
dev.off()
####################################################################
#KL Divergence Analysis
install.packages("entropy")
library("entropy")

express_data = all_data %>% spread(prot_name, Intensity)

gated_data= express_data %>% filter(Ce140 <= (sinh(5.0))) %>% gather(prot_name, Intensity, ATRX:TAL1) 

entropy_disc = function(x, range_vals = c(0, 12), bins = 100)
{
  y1 = discretize(x, numBins=bins, r=range_vals)
  entropy(y1)
}

entro_data = gated_data %>% group_by(Day, prot_name) %>% summarise(entro = entropy_disc(Intensity)) %>% ungroup() %>% 
  
  ggplot(aes(x = Day, y = entro)) +
  geom_point(size=2)+
  facet_wrap(~prot_name, ncol = 4, scales = "free") 

entro_data = gated_data %>% group_by(Day, prot_name) %>% summarise(mean = mean(Intensity)) %>% ungroup() #%>% 

entro_data = entro_data %>% filter(Day == 4, prot_name == "CD123" | prot_name == "CD33" | prot_name == "CD36" | prot_name == "CD41" | prot_name == "CD49f" | prot_name == "CD90" | prot_name == "HBB" | prot_name == "c_Jun" | prot_name == "cyclin_B1") %>% 
  mutate(vol = c(2.2, 5.5, 1.1, 1.1, 1.1, 4.4, 5.5, 6.6, 1.1)) 

entro_data %>% ggplot(aes(x = vol, y = mean)) +
  geom_point(size=2) #+
  facet_wrap(~prot_name, ncol = 4, scales = "free") 

figure = gather_data %>% ggplot(aes(x = Day, y = Intensity, fill=Day)) +
  geom_violin(size=0.1, scale = "width", trim = FALSE, bw = 0.1)  +
  #ylim(0,1000)+
  #scale_y_log10() +  
  facet_wrap(~prot_name, ncol = 4, scales = "free_y") 
####################################################################
pic1 = ggplot(all_data, aes(x = Day, y = Intensity, group = prot_name, colour = prot_name)) +
  scale_y_log10() +  
  stat_summary(fun.y=mean, geom="line")

pdf("/home/qcheng/Research/cytof-erythroid/cytof-erythroid/K1_analysis/K1_mean1.pdf")
plot(pic1)
dev.off()
####################################################################
cd34_data = filter(all_data, prot_name == "CD34")
cd34_data

mean_value = mean(cd34_data$Intensity)
std_value = sd(cd34_data$Intensity)

ggplot(cd34_data, aes(x = (Intensity - mean_value)/std_value)) + 
  geom_histogram(aes(y=..density..), bins=125) +
# scale_x_log10() +
  facet_wrap(~Day, ncol = 3)

ggplot(cd34_data, aes(x = (Intensity - mean_value)/std_value)) + 
  geom_histogram(bins=100) +
  # scale_x_log10() +
  facet_wrap(~Day, ncol = 3) +
  xlim(c(-10,10))

ggplot(cd34_data, aes(x = Intensity)) + 
  geom_freqpoly(aes(y=..density..), bins=125) +
  scale_x_log10() +
  facet_wrap(~Day, ncol = 3)
####################################################################
library(hexbin)

cd34_data = filter(all_data, prot_name == "CD34", Day == "Day12")
cd34_data
cd235_data = filter(all_data, prot_name == "CD235ab", Day == "Day12")
cd235_data

cd38_data = filter(all_data, prot_name == "CD38" , Day == "Day12")
cd38_data

cd41_data = filter(all_data, prot_name == "CD41")
cd41_data


df <- data.frame(log10(cd34_data[1]), log10(cd38_data[1]))
df

df1 = filter(df, Intensity >0, Intensity.1 > 0)
df1
# Create hexbin object and plot
h <- hexbin(df1)
plot(h) 

library(gplots)

# Default call
h2 <- hist2d(df1)
####################################################################
  labs(title = pro_marker[i],
       x = "Day",
       y = "Logarithmic Density",
       colour = "Day") 

pic = ggplot(all_data, aes(x = Day, y = Intensity)) +
  geom_point()  + facet_wrap(~prot_name)+ 
#######################################################
pdf("K1.pdf")
plot(pic)
dev.off()

  p = protein1 %>% gather(Day, value)
  pic = ggplot(p, aes(x = Day, y = value)) +
  geom_violin()  + scale_y_log10() +  
  labs(title = pro_marker[i],
       x = "Day",
       y = "Logarithmic Density",
       colour = "Day")
  #ggsave(pro_marker[i], device = "tiff", path = "/home/qcheng/Research/cytof-erythroid/cytof-erythroid/K1_analysis/" )
  plot_list[[i]] = pic
  
  tiff(pro_marker[i])
  print(plot_list[[i]])
  dev.off()

for (i in 1:26) {
  
  tiff(file_name)
 
}
p
pdf("plots.pdf")
for (i in 1:26) {
  print(plot_list[[i]])
}
dev.off()

#################################################################################
a_cd34 = c(mean(g1$Sm149Di), mean(g2$Sm149Di), mean(g3$Sm149Di), mean(g4$Sm149Di), mean(g5$Sm149Di), 
           mean(g6$Sm149Di), mean(g7$Sm149Di), mean(g8$Sm149Di), mean(g9$Sm149Di), mean(g10$Sm149Di), mean(g11$Sm149Di))

a_cd38 = c(mean(g1$Yb172Di), mean(g2$Yb172Di), mean(g3$Yb172Di), mean(g4$Yb172Di), mean(g5$Yb172Di),
           mean(g6$Yb172Di), mean(g7$Yb172Di), mean(g8$Yb172Di), mean(g9$Yb172Di), mean(g10$Yb172Di), mean(g11$Yb172Di))

a_cd71 = c(mean(g1$Lu175Di), mean(g2$Lu175Di), mean(g3$Lu175Di), mean(g4$Lu175Di), mean(g5$Lu175Di),
           mean(g6$Lu175Di), mean(g7$Lu175Di), mean(g8$Lu175Di), mean(g9$Lu175Di), mean(g10$Lu175Di), mean(g11$Lu175Di))

a_cd36 = c(mean(g1$Gd155Di), mean(g2$Gd155Di), mean(g3$Gd155Di), mean(g4$Gd155Di), mean(g5$Gd155Di),
           mean(g6$Gd155Di), mean(g7$Gd155Di), mean(g8$Gd155Di), mean(g9$Gd155Di), mean(g10$Gd155Di), mean(g11$Gd155Di))

a_cd235 = c(mean(g1$Pr141Di), mean(g2$Pr141Di), mean(g3$Pr141Di), mean(g4$Pr141Di), mean(g5$Pr141Di),
           mean(g6$Pr141Di), mean(g7$Pr141Di), mean(g8$Pr141Di), mean(g9$Pr141Di), mean(g10$Pr141Di), mean(g11$Pr141Di))

t_day = c("Day2", "Day4", "Day6", "Day8", "Day10", "Day11", "Day12", "Day14", "Day16", "Day18", "Day20")

pro_name1 = c("149Sm_CD34", "149Sm_CD34", "149Sm_CD34", "149Sm_CD34", "149Sm_CD34", "149Sm_CD34", "149Sm_CD34", "149Sm_CD34", "149Sm_CD34", "149Sm_CD34", "149Sm_CD34")
pro_name2 = c("172Yb_CD38", "172Yb_CD38", "172Yb_CD38", "172Yb_CD38", "172Yb_CD38", "172Yb_CD38", "172Yb_CD38", "172Yb_CD38", "172Yb_CD38", "172Yb_CD38", "172Yb_CD38")
pro_name3 = c("175Lu_CD71", "175Lu_CD71", "175Lu_CD71", "175Lu_CD71", "175Lu_CD71", "175Lu_CD71", "175Lu_CD71", "175Lu_CD71", "175Lu_CD71", "175Lu_CD71", "175Lu_CD71")
pro_name4 = c("155Gd_CD36", "155Gd_CD36", "155Gd_CD36", "155Gd_CD36", "155Gd_CD36", "155Gd_CD36", "155Gd_CD36", "155Gd_CD36", "155Gd_CD36", "155Gd_CD36", "155Gd_CD36")
pro_name5 = c("141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab", "141Pr_CD235ab")

value = c(a_cd34, a_cd38, a_cd71, a_cd36, a_cd235)
Day = c(t_day, t_day, t_day, t_day, t_day)
Pro_name = c(pro_name1, pro_name2, pro_name3, pro_name4, pro_name5)

value

level1 = data.frame(Day, value, Pro_name) %>% mutate(Day = factor(Day, levels = day_levels))
level1 

level1 %>% ggplot(aes(x = Day, y = value, group = Pro_name, colour = Pro_name)) +
  geom_line() +
  scale_y_log10() +  
  labs(title = "Protein level changes vs time",
       x = "Day",
       y = "Logarithmic Density")

#################################################################################
  geom_jitter(shape=16, alpha = 0.1) + 
  scale_y_log10() +
 
  # geom_boxplot()
  
data.frame(pro1_df) %>% ggplot(aes(x= BCKG190Di)) +
  geom_histogram(bins=20)

data.frame(pro1_df) %>% mutate(cell_id = 1:n()) %>% gather(prot_name, value, -cell_id) %>% 
  ggplot(aes(x= value)) + geom_histogram(bins = 200) + facet_wrap(~prot_name, scales = "free")

data.frame(x_df) %>% mutate(cell_id = 1:n()) %>% gather(prot_name, value, -cell_id) %>% 
  group_by(cell_id) %>% summarise(sum_protein = sum(value)) %>% 
  ggplot(aes(x= sum_protein)) + geom_histogram(bins = 200) + xlim(0, 10000)

# a_cd34_2 = mean(g1$Sm149Di)
# a_cd34_4 = mean(g2$Sm149Di)
# a_cd34_6 = mean(g3$Sm149Di)
# a_cd34_8 = mean(g4$Sm149Di)
# a_cd34_10 = mean(g5$Sm149Di)
# a_cd34_11 = mean(g6$Sm149Di)
# a_cd34_12 = mean(g7$Sm149Di)
# a_cd34_14 = mean(g8$Sm149Di)
# a_cd34_16 = mean(g9$Sm149Di)
# a_cd34_18 = mean(g10$Sm149Di)
# a_cd34_20 = mean(g11$Sm149Di)

# a_cd38_2 = mean(g1$Yb172Di)
# a_cd38_4 = mean(g2$Yb172Di)
# a_cd38_6 = mean(g3$Yb172Di)
# a_cd38_8 = mean(g4$Yb172Di)
# a_cd38_10 = mean(g5$Yb172Di)
# a_cd38_11 = mean(g6$Yb172Di)
# a_cd38_12 = mean(g7$Yb172Di)
# a_cd38_14 = mean(g8$Yb172Di)
# a_cd38_16 = mean(g9$Yb172Di)
# a_cd38_18 = mean(g10$Yb172Di)
# a_cd38_20 = mean(g11$Yb172Di)
uu2 = all_data %>% group_by(prot_name, ) %>% summarise(num_cells = n(), num_zero = sum(Intensity == 0), prop = num_zero/num_cells)
uu = all_data %>% filter(Intensity >0) %>% group_by(prot_name, Day) %>% summarise(min_val = min(Intensity))

