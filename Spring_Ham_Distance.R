install.packages('DeLorean')
install.packages("htmlwidgets")
if(!require("devtools")) install.packages("devtools")
devtools::install_github("bwlewis/rthreejs")
###############################
library("shiny")
library("cytofkit") 

library(flowCore)
library(tidyverse)
library(stringr)

get_pro_data <- function() {
  
  pro1 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c01_sample-1_01_0_Day2_live.fcs", transformation=FALSE)
  pro2 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c02_sample-1_01_0_Day4_live.fcs", transformation=FALSE)
  pro3 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c03_sample-1_01_0_Day6_live.fcs", transformation=FALSE)
  pro4 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c04_sample-1_01_0_Day8_live.fcs", transformation=FALSE)
  pro5 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c05_sample-1_01_0_Day10_live.fcs", transformation=FALSE)
  pro6 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c16_sample-1_01_0_Day12_live.fcs", transformation=FALSE)
  pro7 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c17_sample-1_01_0_Day14_live.fcs", transformation=FALSE)
  pro8 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c18_sample-1_01_0_Day16_live.fcs", transformation=FALSE)
  pro9 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c19_sample-1_01_0_Day18_live.fcs", transformation=FALSE)
  pro10 <- read.FCS("Research/RBC/Data Analysis/cytof-erythroid/cytof-erythroid/K1/export_c20_sample-1_01_0_Day20_live.fcs", transformation=FALSE)

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

  ###############################################################################
  sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f",
                 "CD90", "CD44", "CD41", "CD235ab")

  #other_marker = c("HBB", "HBA", "H3")
  
  tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                  "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B")
  
  #cell_cycle = c("p_HH3", "p_Rb", "cyclin_B1", "IdU")
  
  #DNA = c("Ir191", "Ir193", "Ce140", "Pt195", "Pt194", "Event_Length")
  
  isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di", "Dy164Di",
              "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di",

              #"Nd144Di", "Er170Di", "Yb174Di",
              
              "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
              "Gd158Di", "Sm154Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di", "Er166Di")
              
              #"Nd148Di", "Nd150Di", "Er168Di", "I127Di",
              
              #"Ir191Di",  "Ir193Di", "Ce140Di", "Pt195Di", "Pt194Di", "Event_length")
  
  pro_marker = c(sur_marker, tran_factor)
  
  seq=c(1:26)

  day_levels = c(2,4,6,8,10,12,14,16,18,20)
  
  prot_list = list()
  
  for(i in seq) {
    
    marker = isotops[i]
    
    protein1 = bind_rows( g1[marker] %>% mutate(obstime = 2, col_id = 1:nrow(g1)),
                          g2[marker]%>% mutate(obstime = 4, col_id = 1:nrow(g2)), 
                          g3[marker]%>% mutate(obstime = 6, col_id = 1:nrow(g3)),
                          g4[marker]%>% mutate(obstime = 8, col_id = 1:nrow(g4)),
                          g5[marker]%>% mutate(obstime = 10, col_id = 1:nrow(g5)), 
                          g6[marker]%>% mutate(obstime = 12, col_id = 1:nrow(g6)),
                          g7[marker]%>% mutate(obstime = 14, col_id = 1:nrow(g7)),
                          g8[marker]%>% mutate(obstime = 16, col_id = 1:nrow(g8)),
                          g9[marker]%>% mutate(obstime = 18, col_id = 1:nrow(g9)),
                          g10[marker]%>% mutate(obstime = 20, col_id = 1:nrow(g10))) %>% 
      
      mutate(obstime = factor(obstime, levels = day_levels))
    
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

mean_irid = function(ir1, ir2){
  reg = lm(ir1~ir2)
  (predict(reg) + ir1)/2
}

filter_data = all_data %>% filter(cell_id > 0, cell_id <= 1000)

express_data = all_data %>% spread(cell.type, Intensity)

cell_data = express_data %>% filter(Ce140 <= sinh(5.0)) %>%
  mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193))

gated_data1 = cell_data %>% filter(mean_Pt <= sinh(5.0) | CD235ab >= sinh(6.5))

gated_data = gated_data1 %>% filter(mean_Ir >= sinh(5.0) | CD235ab >= sinh(6.5), 
                                    mean_Ir <= sinh(8.0) | CD235ab >= sinh(6.5))

filter_data = gated_data %>% select(-Ce140, -Ir191, -Ir193, -Pt194, -Pt195, -mean_Ir, -mean_Pt) %>% 
  gather(cell.type, Intensity, ATRX:TAL1) %>% 
  filter(obstime == "MNCs1" | obstime == "MNCs2" | obstime == "MNCs3")


###########The same threshold for all days#############
tfs = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
        "RUNX1", "NFE2p45", "IKZF1", "MAFG", "c_Jun", "CEBPa", "KAT3B",
        
        "CD34", "CD36", "CD38", "CD41", "CD44", 
        "CD71", "CD123", "CD235ab", "CD45RA", "CD49f", "CD90")

sth = c(10, 6, 10, 9, 20, 4, 7, 15,
        10, 11, 20, 20, 11, 13, 12.5,
        
        7, 50, 25, 15, 20, 
        48, 10, 55, 5, 8, 5)

thre_tfs = data.frame(sth = sth) 
rownames(thre_tfs) <- tfs

th_prot = data.frame(cell.type = tfs, threshold = sth)

bin_data = filter_data %>% left_join(th_prot, by = c("cell.type" = "cell.type")) %>% 
  mutate(prot_on = Intensity > threshold)

bin_data = bin_data %>% select(obstime, cell_id, cell.type, cell, prot_on) %>% mutate(prot_on = as.integer(prot_on))

save(bin_data, file = "Desktop/Apr2018/MNCs_bin.RData")
load(file = "Desktop/Apr2018/MNCs_bin.RData") 

cell_code = bin_data %>%
  spread(cell.type, prot_on)

x0 = cell_code %>% select(ATRX:TAL1)

y0 = do.call("paste", x0)

y0 %>% unique()

count_data = cell_code %>% mutate(comb_code = y0) %>% select(obstime, cell, comb_code) 

count0 = count_data %>% group_by(comb_code) %>% summarise(count = n()) 
count1 = count0[order(count0$count, decreasing = T),]

count1 %>% mutate(per = count/sum(count)*100) %>% filter(per >=1)




protein_list = names(cell_code)[-1*(1:3)] 

test0 = paste(as.data.frame(protein_list), sep = " ")
typeof(protein_list) c(1, 2, 3)

count2 = count2 %>% spread(obstime, count)
count2[is.na(count2)] <- 0



x = count1 %>% ungroup() %>% as.tibble() %>% select(-comb_code) %>% rowSums(.)

count1 = count1 %>% ungroup() %>% mutate(sum = x)

count1 %>% filter(sum >= 50)

count0 = count_data %>% select(obstime, comb_code) %>% group_by(comb_code) %>% summarise(count = n()) 
count1 = count0[order(count0$count, decreasing = T),]
count1 = count1 %>% mutate(per = count/sum(count)*100)

count1 %>% ggplot(aes(x= comb_code, y=count))+geom_point()
typeof(count0)


x1 = cell_code %>% select(CD235ab, CD34, CD41)

y1 = do.call("paste", x1)

y1 %>% unique()

test_data = cell_code %>% mutate(comb_code = y1) %>% select(obstime, cell, comb_code) 

test2 = test_data %>% group_by(comb_code) %>% group_by(obstime, add = T) %>% summarise(count = n()) 

test2 = test2 %>% spread(obstime, count)
test2[is.na(test2)] <- 0

test1 = test2[order(test2$`2`, decreasing = T),]

count1 = test_data %>% select(obstime, comb_code) %>%  as.tibble() %>% table() #filter(comb_code != "1 1 0" & comb_code != "0 0 1") %>%
count2 = count1 %>% as.tibble() %>% group_by(obstime) %>% mutate(prob = n/sum(n)) %>% ungroup() %>% mutate(weight = prob) #%>% ggplot(aes(x=as.numeric(obstime), y=n/total)) + geom_line() + facet_wrap(~comb_code)

count3 = count2 %>% as.tibble() %>% mutate(id = 1:n(), obstime = as.numeric(obstime)) %>% select(obstime, comb_code, weight) %>% spread(obstime, weight)










protein = names(gated_data)[-1*(1:4)] 

protein = protein[-c(13, 21, 22, 27, 28, 32, 33)]

library(rlang)

new_gated = filter_data



tukey_lim = function(x, k = 3){
  q1_q3 = quantile(x, c(0.75, 0.9))
  upper_lim = q1_q3[2] + k*(q1_q3[2]-q1_q3[1])
  #if upper_lim 
  return(upper_lim)
}

for(i in protein){
  print(i)
  gated_data = express_data %>% rename(prot_count = !!!sym(i)) %>% group_by(obstime) %>% 
    filter(prot_count<tukey_lim(prot_count, 3)) %>% ungroup()
  
  colnames(gated_data)[which(names(gated_data) == "prot_count")] <- i
}

gated_data %>% group_by(obstime) %>% summarise(count = n())

check_cell_cycle <- function(mean_Ir, p_Rb, IdU){
  
  cycle_state = 0
  
  if(p_Rb <= sinh(5.0) & mean_Ir <= sinh(6.15)) cycle_state = 0
  if(p_Rb > sinh(5.0) & IdU < sinh(3.5) & mean_Ir <= sinh(6.15)) cycle_state = 1
  if(p_Rb > sinh(5.0) & IdU >= sinh(3.5)) cycle_state = 2
  if(IdU < sinh(3.5) & mean_Ir > sinh(6.15) & p_Rb > sinh(5.0)) cycle_state = 3
  
  return(cycle_state)
}

gated_data = gated_data3 %>% mutate(cycle_state = mapply(check_cell_cycle, mean_Ir, p_Rb, IdU)) %>% select(cell, obstime, cycle_state, tfs) #GATA1, PU1, ATRX, c_Myc, KLF1, FLI1, TAL1, GATA2, RUNX1, NFE2p45, IKZF1, MAFG, c_Jun, CEBPa, KAT3B, MEF2C) #%>% gather(protein, Intensity, ATRX:mean_Ir) #%>% group_by(protein)



###########Different thresholds for different days#############
th_prot<- read.csv(file="Desktop/th_prot.csv", header=TRUE, sep=",")

th_prot = th_prot %>% as.tibble()

long_th_prot = th_prot %>% gather(cell.type, threshold, -Day)
gated_data = all_data %>% left_join(long_th_prot, by = c("obstime" = "Day", "cell.type" = "cell.type")) %>% 
  mutate(prot_on = Intensity > threshold)

binary_name = gated_data %>% select(obstime, cell.type, cell, prot_on) %>% mutate(prot_on = as.integer(prot_on))

#expression data
# #%>% filter(Ce140 <= sinh(5.0)) %>% #remove bead-cell doublets
#  # mutate(mean_Pt = mean_irid(Pt194, Pt195), mean_Ir = mean_irid(Ir191, Ir193))
# #*************************************************#
# #for igraph data
# # 
# # gated_data1 = express_data %>% filter(mean_Pt <= hill_asinh(1e6, sinh(5.0), 5) | HBA >= hill_asinh(1e6, sinh(6.5), 5) | CD235ab >= hill_asinh(1e6, sinh(6.5), 5))
# # 
# # gated_data = gated_data1 %>% filter(mean_Ir >= hill_asinh(1e6, sinh(5.0), 5) | HBA >= hill_asinh(1e6, sinh(6.5), 5) | CD235ab >= hill_asinh(1e6, sinh(6.5), 5),
# #                                     mean_Ir <= hill_asinh(1e6, sinh(8.0), 5) | HBA >= hill_asinh(1e6, sinh(6.5), 5) | CD235ab >= hill_asinh(1e6, sinh(6.5), 5)) %>% 
# #   select(cell, obstime, ATRX:CD90, CEBPa:CXCR4, FLI1:HBB, IKZF1, KAT3B:NFE2p45, PU1:TAL1) #%>% gather(protein, Intensity, ATRX:mean_Ir) #%>% group_by(protein)
# 
# gated_data1 = express_data %>% filter(mean_Pt <= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5))
# 
# gated_data2 = gated_data1 %>% filter(mean_Ir >= sinh(5.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5),
#                                     mean_Ir <= sinh(8.0) | HBA >= sinh(6.5) | CD235ab >= sinh(6.5)) 

protein = names(express_data)[-1*(1:2)] 

protein = protein[-c(1, 2)]

library(rlang)

tukey_lim = function(x, k = 3){
  q1_q3 = quantile(x, c(0.75, 0.9))
  upper_lim = q1_q3[2] + k*(q1_q3[2]-q1_q3[1])
  #if upper_lim 
  return(upper_lim)
}

for(i in protein){
  print(i)
  gated_data = express_data %>% rename(prot_count = !!!sym(i)) %>% group_by(obstime) %>% 
    filter(prot_count<tukey_lim(prot_count, 3)) %>% ungroup()
  
  colnames(gated_data)[which(names(gated_data) == "prot_count")] <- i
}

gated_data %>% group_by(obstime) %>% summarise(count = n())

check_cell_cycle <- function(mean_Ir, p_Rb, IdU){
  
  cycle_state = 0
  
  if(p_Rb <= sinh(5.0) & mean_Ir <= sinh(6.15)) cycle_state = 0
  if(p_Rb > sinh(5.0) & IdU < sinh(3.5) & mean_Ir <= sinh(6.15)) cycle_state = 1
  if(p_Rb > sinh(5.0) & IdU >= sinh(3.5)) cycle_state = 2
  if(IdU < sinh(3.5) & mean_Ir > sinh(6.15) & p_Rb > sinh(5.0)) cycle_state = 3
  
  return(cycle_state)
}

gated_data = gated_data3 %>% mutate(cycle_state = mapply(check_cell_cycle, mean_Ir, p_Rb, IdU)) %>% select(cell, obstime, cycle_state, tfs) #GATA1, PU1, ATRX, c_Myc, KLF1, FLI1, TAL1, GATA2, RUNX1, NFE2p45, IKZF1, MAFG, c_Jun, CEBPa, KAT3B, MEF2C) #%>% gather(protein, Intensity, ATRX:mean_Ir) #%>% group_by(protein)

new_gated = all_data %>% mutate(obstime = as.numeric(as.character(obstime))) #gated_data %>% gather(prot_name, Intensity, GATA1:HBA)

cor_signal <- c()

for(i in 1:nrow(new_gated)) {
  
  
  correct_vals = remove_noise(new_gated[i,]$Intensity, new_gated[i,]$cell.type, new_gated[i,]$obstime)
  
  cor_signal = append(cor_signal, correct_vals)
  
}

new_gated1 = new_gated %>% mutate(Bi_Inten = cor_signal) 
new_gated = new_gated %>% mutate(Bi_Inten = cor_signal) %>% select(-cycle_state) #%>% spread(prot_name, Intensity)


remove_noise <- function(intensity, prot_name){
  
  signal_th = thre_tfs[prot_name,]
  
  if(intensity < signal_th) vals = 0
  else vals = 1 #intensity
  
  return(vals) 
}

cor_signal <- c()

for(i in 1:nrow(new_gated)) {
  
  correct_vals = remove_noise(new_gated[i,]$Intensity, new_gated[i,]$cell.type)
  
  cor_signal = append(cor_signal, correct_vals)
  
}

new_gated1 = new_gated %>% mutate(Bi_Inten = cor_signal) 
new_gated = new_gated %>% mutate(Bi_Inten = cor_signal) %>% select(-cycle_state) #%>% spread(prot_name, Intensity)


# gated_data %>% ggplot(aes(x = Ir193)) +
#   geom_histogram(bins=150)+
#   geom_vline(aes(xintercept = 305.4092), colour = "red") + 
#   geom_vline(aes(xintercept = 577.7555), colour = "red") + 
#   facet_wrap(~cycle_state, ncol = 1, scales = "free_y") +
#   xlim(0, 900)
# 
# xx = gated_data %>% filter(cycle_state==0.0) 
# mean(xx$Ir193)

new_gated %>% ggplot(aes(x=obstime, y = asinh(Intensity))) + geom_point(aes(colour = factor(Bi_Inten)), size = 0.5) + facet_wrap(~prot_name) 

new_gated %>% ggplot(aes(x=obstime, y = asinh(Intensity))) + geom_violin(aes(colour = factor(Bi_Inten)), size = 0.5) + facet_wrap(~prot_name) 

new_gated %>% select(-Bi_Inten) %>% spread(prot_name, Intensity) %>% 
  ggplot(aes(x= asinh(CD71), y = asinh(CD41)))  + geom_point() + facet_wrap(~obstime) #+ geom_count(aes(size = ..prop..))

 new_gated %>% select(-Intensity) %>% spread(prot_name, Bi_Inten) %>% select(obstime, CD235ab, HBA) %>% 
  group_by(obstime) %>% summarise(sum(HBA==1&CD235ab==0)/n()*100)
#######################################################

cell_code = binary_name %>%
  spread(cell.type, prot_on) %>% select(-KAT3B) #, -cycle_state)

x0 = cell_code %>% select(ATRX:TAL1)

y0 = do.call("paste", x0)

y0 %>% unique()

count_data = cell_code %>% mutate(comb_code = y0) %>% select(obstime, cell, comb_code) 

count2 = count_data %>% group_by(comb_code) %>% group_by(obstime, add = T) %>% summarise(count = n()) 

count2 = count2 %>% spread(obstime, count)
count2[is.na(count2)] <- 0

count1 = count2[order(count2$`2`, decreasing = T),]

x = count1 %>% ungroup() %>% as.tibble() %>% select(-comb_code) %>% rowSums(.)

count1 = count1 %>% ungroup() %>% mutate(sum = x)

count1 %>% filter(sum >= 50)

count0 = count_data %>% select(obstime, comb_code) %>% group_by(comb_code) %>% summarise(count = n()) 
count1 = count0[order(count0$count, decreasing = T),]
count1 = count1 %>% mutate(per = count/sum(count)*100)

count1 %>% ggplot(aes(x= comb_code, y=count))+geom_point()
typeof(count0)


x1 = cell_code %>% select(CD235ab, CD34, CD41)

y1 = do.call("paste", x1)

y1 %>% unique()

test_data = cell_code %>% mutate(comb_code = y1) %>% select(obstime, cell, comb_code) 

test2 = test_data %>% group_by(comb_code) %>% group_by(obstime, add = T) %>% summarise(count = n()) 

test2 = test2 %>% spread(obstime, count)
test2[is.na(test2)] <- 0

test1 = test2[order(test2$`2`, decreasing = T),]




count1 = test_data %>% select(obstime, comb_code) %>%  as.tibble() %>% table() #filter(comb_code != "1 1 0" & comb_code != "0 0 1") %>%
count2 = count1 %>% as.tibble() %>% group_by(obstime) %>% mutate(prob = n/sum(n)) %>% ungroup() %>% mutate(weight = prob) #%>% ggplot(aes(x=as.numeric(obstime), y=n/total)) + geom_line() + facet_wrap(~comb_code)

count3 = count2 %>% as.tibble() %>% mutate(id = 1:n(), obstime = as.numeric(obstime)) %>% select(obstime, comb_code, weight) %>% spread(obstime, weight)

x1 = c(0,1,1)
x2 = c(1,0,1)
hamming.distance(x1,x2)

install.packages("e1071")
library(e1071)
library(igraph)
library(Rtsne)
library('visNetwork')
A = y
B = y

vectorized_pdist <- function(A,B){
  
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
  m = nrow(A)
  n = nrow(B)
  
  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B))
}

vectorized_hdist <- function(A,B){
  
  m = nrow(A)
  n = nrow(B)
  
  tmp = matrix(0, ncol = n, nrow=m)

  for(i in 1:m) {
    for(j in 1:n){
  
      x1 = A[i,] 
      x2 = B[j,] 
      tmp[i,j] =  hamming.distance(x1, x2)
      
    }
    
  }
  return(tmp)
}

comb_code = count3 %>% select(comb_code) %>% data.frame()
comb_code[1,]
#gated_data = count_data %>% mutate(Day_index = apply(data.frame(obstime), 1, check_day))
# gated_vals = count %>% select(comb_code) %>% data.frame()
# rownames(gated_vals) <- 1:nrow(gated_vals)
# new_vals = gated_vals %>% as.matrix()

x = x0 %>% unique() %>% as.tibble() %>% as.matrix()
x = data.frame(c(0, 0, 0), c(0, 1, 0), c(0, 1, 1), c(1, 0, 0)) %>% as.matrix() %>% t()

dist_matrix1 = vectorized_hdist(x, x)
dist_matrix1 = dist_matrix1 %>% as.matrix()
diag(dist_matrix1) = 1e16

y = count3 %>% select(`2`:`20`) %>% as.matrix()
dist_matrix2 = vectorized_pdist(y, y)
dist_matrix2 = dist_matrix1 %>% as.matrix()
diag(dist_matrix2) = 1e16

colnames(dist_matrix1) <- 1:4 #nrow(count3)
comb_index = c(1:4)

uu2       =  as_tibble(dist_matrix1) %>% map(function(x) order(x, na.last = NA)[1])
knn_piars = bind_cols(uu2) %>% gather(node_from, nn) %>% mutate(node_from = factor(node_from)) %>% filter(node_from != nn)

cc_data = data.frame(x)
nodes     = bind_cols(cell_node = 1:4,comb_index=comb_index) #nrow(count3), comb_code = count3$comb_code, )

cell_id = select(new_gated, cell) %>% mutate(node_id = 1:nrow(new_gated)) %>%  data.frame()
nrow(cell_id)

net = graph_from_data_frame(d=knn_piars, vertices=nodes, directed=TRUE)
net[]

# visNetwork(nodes, data.frame(l), width="100%", height="400px", main="Cell Trajectory") %>% 
# visNodes(shape = "square", title = "node") %>% 
# visEdges(smooth = list(enabled = TRUE, type = "diagonalCross"))
# 
# clp <- cluster_label_prop(net)

# V(net)$size <-2.5
# V(net)$frame.color <- "white"
# V(net)$arrow.mode <- 0
# V(net)$community
# E(net)
# neigh.nodes <- neighbors(net, V(net)[146:214], mode="out")
# vcol <- rep("grey40", vcount(net))
# vcol[V(net)[146:214]] <- "red"
day_levels = x %>% unique()

l = layout_with_fr(net, dim=2)

V(net)$cell_id[V(net)$Day_index == "0"]

colrs <- c("blue", "firebrick1", "burlywood", "lightgoldenrod1") #, "chartreuse", "darkorchid") #"deeppink", "hotpink1", "olivedrab1", "cyan", "darkorange", "mediumorchid2")

dim(colrs)
V(net)$color <- colrs[V(net)$comb_index]

plot(net, layout = l, edge.arrow.size=.2, vertex.label=NA, vertex.size = 4.0, edge.size = 0.4)#, vertex.color=vcol)

legend(x=-1.5, y=1.1, comb_index, pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1) #%>% levels(day_levels)

new_data = count_data %>% group_by(obstime) %>% count(comb_code) %>% spread(comb_code, n, fill = 0) %>% ungroup()

library(pheatmap)

new_data1 = new_data %>% select(-obstime) %>% as.matrix() %>% magrittr::set_rownames(new_data$obstime) %>% t()#remove_rownames() %>% column_to_rownames(var = "obstime")

pheatmap(new_data1, cluster_cols  =FALSE)

new_data1 = new_data %>% as.matrix() %>% t() %>% as.data.frame() 

#colnames(new_data1) <- data.frame(new_data1[1,])

new_data1 = new_data1 %>% tibble::rownames_to_column()

new_data2 = new_data1[-1,] %>% as.tibble()

new_data2 = arrange(new_data2, desc(new_data2$`1`))

count_data %>% ggplot() + geom_bin2d(aes(x=obstime, y=comb_code))

####################
norm_signal <- c()

for(i in 1:nrow(new_gated1)) {
  
  norm_vals = data_normalise(new_gated1[i,]$Intensity, new_gated1[i,]$prot_name)
  
  norm_signal = append(norm_signal, norm_vals)
  
}

new_gated0 = new_gated1 %>% mutate(Intensity = norm_signal) %>% spread(prot_name, Intensity)

data_normalise <- function(intensity, prot_name){
  
  x = new_gated[prot_name] %>% unlist()
  upper_lim = tukey_lim(x, 3)
  
  vals = (intensity-0)/(upper_lim - 0)
  
  return(vals) 
}

new_gated_data = new_gated0 %>% gather(prot_name, Intensity, ATRX:TAL1)

new_gated_data = new_gated_data %>% filter(Intensity <= 1.0)


new_gated = new_gated_data %>% spread(prot_name, Intensity) #_data

new_gated = new_gated[complete.cases(new_gated), ] #%>% mutate(obstime = factor(obstime, levels = c(0,2,4,6,8,10,11,12,14,16,18,20)))


new_gated %>% ggplot(aes(x = CD34))+ 
  geom_histogram(bins = 100) +
  scale_y_log10() + 
  facet_wrap(~ obstime, ncol = 4, scales = "free")

new_gated %>% spread(prot_name, Intensity) %>% group_by(obstime) %>% summarise(cd34 = mean(CD34), klf1 = mean(KLF1), fli1 = mean(FLI1), gata1 = mean(GATA1), pu1 = mean(PU1), CD235 = mean(CD235ab)) %>% ungroup() %>%  
  #ggplot(aes(x = obstime, y = cd34))+
  ggplot(aes(x = obstime))+
  geom_point(aes(y = fli1),size=2, colour="red")+
  geom_point(aes(y = cd34),size=2, colour="green")+
  geom_point(aes(y = klf1),size=2, colour="blue")+
geom_point(aes(y = gata1),size=2, colour="yellow")+
geom_point(aes(y = pu1),size=2, colour="black")+
  legend()


+
  geom_point(aes(y = CD235),size=2, colour="grey")


  
  ggplot(aes(x = Day, y = entro)) +
  geom_point(size=2)+
  facet_wrap(~prot_name, ncol = 4, scales = "free") 


new_gated %>% spread(prot_name, Intensity) %>% select(cell, obstime, KLF1, FLI1, GATA1, PU1, CD34, HBA, CD235ab) %>% gather(prot_name, Intensity, KLF1:CD235ab) %>%  #%>% mutate(obstime = factor(obstime, levels = c(0,2,4,6,8,10,11,12,14,16,18,20))) %>% #mutate(obstime = as.character(obstime)) %>%  
  ggplot(aes(x = obstime, y = arcsinh(Intensity), group = obstime, fill=obstime)) +
  geom_violin(size=0.1, scale = "width", bw = 0.001)  +
  facet_wrap(~prot_name, ncol = 1, scales = "free") 

cell_id = select(new_gated, cell) %>% mutate(node_id = 1:nrow(new_gated)) %>%  data.frame()

gated_vals = new_gated %>% select(ATRX:TAL1) %>% data.frame()
rownames(gated_vals) <- 1:nrow(gated_vals)

new_vals = gated_vals %>% as.matrix()

library(tidyverse)
library(igraph)
library(Rtsne)
#install.packages('visNetwork')
library('visNetwork')

vectorized_pdist <- function(A,B){
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))
  
  m = nrow(A)
  n = nrow(B)
  
  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  sqrt( tmp - 2 * tcrossprod(A,B))
}

day_vals = unique(new_gated$obstime) 

day_vals = c(0,2,4,6,8,10,11,12,14,16,18,20)

check_day <- function(d) {
  dd = which(day_vals == d)
  dd
}

check_day(2)

gated_data = new_gated %>% mutate(Day_index = apply(data.frame(obstime), 1, check_day))

dist_matrix = vectorized_pdist(new_vals, new_vals)
diag(dist_matrix) = 1e16

uu2       =  as_tibble(dist_matrix) %>% map(function(x) order(x, na.last = NA)[1:5])
knn_piars = bind_cols(uu2) %>% gather(node_from, nn) %>% mutate(node_from = factor(node_from)) %>% filter(node_from != nn)

cc_data = data.frame(new_vals)
nodes     = bind_cols(cell_node = 1:nrow(cc_data), cell_id = gated_data$cell, Day = gated_data$obstime, 
                      Day_index = gated_data$Day_index,  cc_data)

nrow(cell_id)
net = graph_from_data_frame(d=knn_piars, vertices=nodes, directed=T)
net[]
# visNetwork(nodes, data.frame(l), width="100%", height="400px", main="Cell Trajectory") %>% 
# visNodes(shape = "square", title = "node") %>% 
# visEdges(smooth = list(enabled = TRUE, type = "diagonalCross"))
# 
# clp <- cluster_label_prop(net)

# V(net)$size <-2.5
# V(net)$frame.color <- "white"
# V(net)$arrow.mode <- 0
# V(net)$community
# E(net)
# neigh.nodes <- neighbors(net, V(net)[146:214], mode="out")
# vcol <- rep("grey40", vcount(net))
# vcol[V(net)[146:214]] <- "red"

l = layout_with_fr(net, dim=2)

V(net)$cell_id[V(net)$Day_index == "0"]

colrs <- c("blue", "firebrick1", "burlywood", "lightgoldenrod1", "chartreuse", "darkorchid",
           "deeppink", "hotpink1", "olivedrab1", "cyan", "darkorange", "mediumorchid2")

V(net)$color <- colrs[V(net)$Day_index]

plot(net, layout = l, edge.arrow.size=.2, vertex.label=NA, vertex.size = 4.0, edge.size = 0.4)#, vertex.color=vcol)
legend(x=-1.5, y=1.1, day_levels, pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=.8, bty="n", ncol=1)#%>% levels(c(0,2,4,6,8,10,11,12,14,16,18,20))

mean_distance(net, directed=T)

dist.from.STEM <- distances(net, v=V(net), to=V(net), weights=NA)

oranges <- colorRampPalette(c("dark red", "gold"))

col <- oranges(max(dist.from.STEM)+1)
col <- col[dist.from.STEM+1]

plot(net, vertex.color=col, vertex.label=dist.from.STEM, edge.arrow.size=.6, 
     
     vertex.label.color="white")


#c(0, 10, 12, 16, 18, 20,  4, 6,  8, 11,  2, 14),
unique(V(net)$Day_index)
V(net)$Day_index
plot(net, layout =l, edge.arrow.size=.2, vertex.label=NA, vertex.size = 4.0, edge.size = 0.4) #, vertex.color=vcol)

plot(net, layout = l, edge.arrow.size=.2, vertex.label=NA, vertex.size = 2.5, edge.size = 0.4, 
     vertex.color=c("blue",  "red")[1+(V(net)$HBA >= 0.15)]) #, vertex.color=vcol)

legend(x=-1.5, y=-1.1, c("High", "Low"), pch=21,
       col="#777777", pt.bg=c("red", "blue"), pt.cex=2, cex=.8, bty="n", ncol=1)

tkid <- tkplot(net)

tk_off()

library(threejs)

graphjs(net, layout=l, vertex.size = 0.8, edge.size = 0.4)#, vertex.color=c("red", "skyblue")[1+(V(net)$CD235ab >= 0.5)])

i <- cluster_optimal(net)

set.seed(1)
g <- sample_islands(3, 10, 5/10, 1)
i <- cluster_optimal(g)
graphjs(g, vertex.color=c("orange", "green", "blue")[i$membership], vertex.shape="sphere")

dev.off()
###########################

gene1 = express_data$cyclin_B1
gene2 = express_data$IdU
gene3 = express_data$p_HH3
gene4 = express_data$p_Rb

value_5_95 <- function(x) { 
  
  value = quantile(x, probs = c(0.05, 0.95))
  
  return(value)}

# express_data = express_data %>% mutate(cyclin_B1 = (cyclin_B1)/(value_5_95(gene1)[2]-value_5_95(gene1)[1]))
# express_data = express_data %>% mutate(IdU = (IdU)/(value_5_95(gene2)[2]-value_5_95(gene2)[1]))
# express_data = express_data %>% mutate(p_HH3 = (p_HH3)/(value_5_95(gene3)[2]-value_5_95(gene3)[1]))
# express_data = express_data %>% mutate(p_Rb = (p_Rb)/(value_5_95(gene4)[2]-value_5_95(gene4)[1]))

ggplot(express_data, aes(express_data$p_HH3, express_data$IdU)) +
  geom_point(alpha = 0.1)+
  geom_density_2d()+
  #stat_density_2d(aes(fill = ..level..), , geom = "polygon") +
  facet_wrap(~obstime) + #, scales = "free") 
  labs(x = "p HH3", y = "IdU")

#express_data = express_data %>% mutate(ratio1 = KLF1-FLI1, ratio2 = GATA1-PU1)

# express_data %>% ggplot() +
#   geom_jitter(mapping = aes(x=obstime, y = CD41), shape = 16, alpha = 0.5) +
#   geom_smooth(mapping = aes(x=obstime, y = CD41))
# 
# express_data %>% ggplot() +
#   geom_jitter(mapping = aes(x=obstime, y = ratio2), shape = 16, alpha = 0.5)
express_data = gated_data
express_data$obstime <- NULL
express_data$cell_id <- NULL
express_data$cell <- NULL
express_data$capture <- NULL

#x1=t(express_data)
x1=express_data
write.table(x1, file='SPRING/datasets/rbc/new_test.csv',quote=FALSE, sep='\t',col.names = FALSE, row.names = FALSE)

#cell groupings according to days
day_data = filter_data %>% spread(cell.type, Intensity)
y <- gated_data$obstime
x2 = t(y)
rownames(x2) <- 'Day'
#write.csv(x, "K3.csv", col.names = FALSE)
write.table(x2, file='SPRING/datasets/rbc/new_day.csv',quote=FALSE, sep=',',col.names = FALSE)

#################gene markers##########################
sur_marker = c("CD34", "CD36", "CD71", "CD38", "CD45RA", "CD123", "CD49f", "CD90", "CD44", "CD41", "CD235ab")

tran_factor = c("GATA1", "PU1", "ATRX", "c_Myc", "KLF1", "FLI1", "TAL1", "GATA2", 
                "RUNX1", "NFE2", "BACH1", "IKZF1", "MAFG", "C-Jun", "CEBPa")

isotops = c("Sm149Di", "Gd155Di", "Lu175Di", "Yb172Di", "Nd143Di", "Eu151Di",
            "Dy164Di", "Dy161Di", "Eu153Di", "Y89Di", "Pr141Di",
            "Gd156Di", "Er167Di", "Gd160Di", "Yb176Di", "Ho165Di", "Tm169Di", "Dy162Di", "Dy163Di",
            "Gd158Di", "Sm154Di", "Yb171Di", "Tb159Di", "Sm152Di", "Yb173Di", "Nd145Di")

pro_marker = data.frame(gene = c(sur_marker, tran_factor)) %>% as.tibble(.) %>%t(.)
pro_marker
sort_marker <- pro_marker[, order(pro_marker)]
sort_marker = data.frame(gene = sort_marker)
sort_marker

#pro_marker = data.frame(gene = c("CD34"))
########################################################
library(DeLorean)
library(dplyr)

# data(GuoDeLorean)
# class(colnames(guo.expr))
# guo.expr
# guo.cell.meta$cell# %>% group_by(cell, capture)
# guo.gene.meta
# guo.cell.meta
######################################################
# Limit number of cores to 2 for CRAN
options(DL.num.cores=min(default.num.cores(), 10))
all_data$cell_id
#filter_data = filter(all_data, all_data$Intensity != 0, all_data$Intensity < 20, all_data$obstime <10, all_data$cell_id <= 10)
filter_data1 = filter(all_data, all_data$cell_id > 0, all_data$cell_id <= 10)#, all_data$cell_id != 2) #%>% mutate(Intensity = log10(Intensity+1e-4))
filter_data1

# filter_data = filter(filter_data1, filter_data1$cell != "6 2", filter_data1$cell != "8 2", filter_data1$cell != "2 4",
#                      filter_data1$cell != "4 4", filter_data1$cell != "7 4", filter_data1$cell != "10 10",
#                      filter_data1$cell != "3 11", filter_data1$cell != "4 12", filter_data1$cell != "10 2",
#                      filter_data1$cell != "2 12", filter_data1$cell != "9 4", filter_data1$cell != "10 6",
#                      filter_data1$cell != "2 8", filter_data1$cell != "9 14", filter_data1$cell != "8 6",
#                      filter_data1$cell != "5 4", filter_data1$cell != "8 14", filter_data1$cell != "1 11",
#                      filter_data1$cell != "9 6", filter_data1$cell != "20 6", filter_data1$cell != "13 10", filter_data1$cell != "17 10",
#                      filter_data1$cell != "19 10", filter_data1$cell != "19 11", filter_data1$cell != "23 4",
#                      filter_data1$cell != "17 12", filter_data1$cell != "23 12", filter_data1$cell != "23 14",
#                      filter_data1$cell != "24 12", filter_data1$cell != "12 12", filter_data1$cell != "16 12",
#                      filter_data1$cell != "17 12", filter_data1$cell != "16 10", filter_data1$cell != "15 12",
#                      filter_data1$cell != "20 12", filter_data1$cell != "19 4", filter_data1$cell != "25 2",
#                      filter_data1$cell != "16 11", filter_data1$cell != "12 16", filter_data1$cell != "21 11",
#                      filter_data1$cell != "25 18", filter_data1$cell != "24 0", filter_data1$cell != "13 20",
#                      filter_data1$cell != "7 0", filter_data1$cell != "20 0", filter_data1$cell != "25 0" )

# filter_data = filter(filter_data1, filter_data1$cell != "20 6", filter_data1$cell != "13 10", filter_data1$cell != "17 10",
#                      filter_data1$cell != "19 10", filter_data1$cell != "19 11", filter_data1$cell != "23 4",
#                      filter_data1$cell != "17 12", filter_data1$cell != "23 12", filter_data1$cell != "23 14",
#                      filter_data1$cell != "24 12", filter_data1$cell != "12 12", filter_data1$cell != "16 12",
#                      filter_data1$cell != "17 12", filter_data1$cell != "16 10", filter_data1$cell != "15 12",
#                      filter_data1$cell != "20 12", filter_data1$cell != "19 4", filter_data1$cell != "25 2",
#                      filter_data1$cell != "16 11", filter_data1$cell != "12 16", filter_data1$cell != "21 11",
#                      filter_data1$cell != "25 18", filter_data1$cell != "24 0")

filter_data =filter_data1

express_data = filter_data %>% spread(cell.type, Intensity)
express_data
new_tbl = express_data %>% gather(prot, value, ATRX:TAL1) %>% select(cell, prot:value) %>% spread(cell, value)
new_tbl
new_mat = new_tbl %>% select(-prot)%>% as.matrix(.)
rownames(new_mat) = new_tbl$prot

cell_data = express_data
cell_data

cell_data = cell_data %>% select(obstime:capture) %>% mutate(cell.type = "NA", capture = factor(capture)) %>% 
  arrange(cell)
# order(cell_data$cell)
# cell_data <- cell_data[order(cell_data$cell),]
cell_data
# keeps1 <- c("obstime")
# keeps2 <- c("cell")
# exp1 = express_data[-keeps1]
# exp2 = express_data[keeps2]
# exp2 = t(exp2)
# 
# express_data1 = t(as.matrix(exp1, nrow= 1))
# colnames(express_data1) <- exp2
# rownames(express_data1) <- "CD34"
# exp1
#express_data = cd34_data %>% spread(Protein, Express) %>% mutate(cell = 1:n())
new_mat
pro_marker

nrow(new_mat)
nrow(filter_data)
colnames(new_mat)
cell_data$cell
rownames(new_mat)
ncol(new_mat)
nrow(filter_data)
cell_data$cell
x = is.na(cell_data)
View(x)

dl <- de.lorean(new_mat, sort_marker, cell_data)
#dl1 <- de.lorean(new_mat, sort_marker, cell_data)
dl
dl$gene
dl <- estimate.hyper(
  dl,
  sigma.tau=0.5,
  length.scale=20.0,
  model.name='exact',
  adjust.cell.sizes = FALSE
  )

# num.at.each.stage <- 5
# epi.sampled.cells <- guo.cell.meta %>%
#   filter(capture < "32C" |
#            "EPI" == cell.type |
#            "ICM" == cell.type) %>%
#   group_by(capture) %>%
#   do(sample_n(., num.at.each.stage))
# 
# dl <- filter_cells(dl, cells=epi.sampled.cells$cell)
# 
# dl <- aov.dl(dl)
# 
# head(dl$aov)
# 
# tail(dl$aov)
# 
# dl <- filter_genes(dl, genes=head(dl$aov, 20)$gene)

dl <- fit.dl(dl, method='vb')

dl <- examine.convergence(dl)

plot(dl)

plot(dl, type='pseudotime')

dl <- make.predictions(dl)

pp1 = plot(dl, type='profiles')
 
pp1 + facet_wrap(~gene, scales = "free") 

dl$predictions

+ geom_point(colour = cell_data$cell_id)
############################
filter_data %>% ggplot(aes(x = obstime, y = Intensity)) + geom_point() +facet_wrap(~cell.type)
filter_data %>% ggplot(aes(x = obstime, y = Intensity)) + geom_point(position = "jitter") +facet_wrap(~cell.type)
filter_data %>% ggplot(aes(x = obstime, y = Intensity)) + geom_point(position = "jitter", alpha = 0.4) +facet_wrap(~cell.type, scale = "free")

