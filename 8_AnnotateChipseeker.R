rm(list = ls())

library(ggplot2)
library(scales)
library(reshape2)

setwd("~/Documents/PCCR/MethylRAD/F9")

# ************************************ Create Plots of Changes for LSD1 ********************************************** #

CreateChangePlot_LSD1 <- function(subDat_fin, figPreName, type){
  colnames(subDat_fin) <- c('GeneName', 'UD', 'Trt', 'Trt_Min_UD', 'Trt_Min_UD_PerChange')
  
  subDat_fin$Pos <- ifelse(subDat_fin$Trt_Min_UD > 0, TRUE, FALSE)
  subDat_fin <- subDat_fin[order(subDat_fin$Trt_Min_UD, decreasing = T), ]
  
  subDat_fin$GeneName <- factor(subDat_fin$GeneName, levels = subDat_fin$GeneName[order(subDat_fin$Trt_Min_UD, decreasing = T)])
  subDat_fin$Pos <- as.factor(subDat_fin$Pos)
  
  # Waterfall Plot
  b <- ggplot(subDat_fin, aes(x = GeneName, y = Trt_Min_UD, color = Pos)) +
    scale_fill_discrete(name = "Positive Changes") + 
    theme_classic() + 
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(face="bold", angle = 90), legend.position = "none") 
  b <- b + geom_bar(stat = 'identity', position = position_dodge(width = 0.4)) +
    scale_color_manual(values =  c("#E69F00", "#009E73"))
  
  if (type == 'PG'){
    b <- b + labs(title = "Waterfall plot for changes in Methylation between F9_D4_PG and F9_UD", x = NULL, y = "Quantity of Change") +
      geom_hline(yintercept = c(100, -100), linetype = "dashed", color = "red")
  } 
  
  if (type == 'TCP'){
    b <- b + labs(title = "Waterfall plot for changes in Methylation between F9_D4_TCP and F9_UD", x = NULL, y = "Quantity of Change") +
      geom_hline(yintercept = c(100, -100), linetype = "dashed", color = "red")
  }
  
  if (type == 'None'){
    b <- b + labs(title = "Waterfall plot for changes in Methylation between F9_D4 and F9_UD", x = NULL, y = "Quantity of Change") +
      geom_hline(yintercept = c(100, -100), linetype = "dashed", color = "red")
  }
 
  pdf(paste("Figures/", figPreName, "_LSD1_waterfallPlot.pdf", sep = ""), width = 8, height = 6)
  print(b)  
  dev.off()
  
  # PieChart
  pos <- sum(subDat_fin$Trt_Min_UD > 0)
  neg <- sum(subDat_fin$Trt_Min_UD < 0)
  nch <- sum(subDat_fin$Trt_Min_UD == 0)
  df <- data.frame(group = c("Positive", "NoChange", "Negative"),
                   value = c(pos, nch, neg))
  df$perc <- df$value/sum(df$value)
  
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  pie <- ggplot(data = df, aes(x = "", y = value, fill = group)) + 
    geom_bar(stat = "identity") + coord_polar("y", start = 0) + blank_theme + 
    theme(axis.text.x = element_blank()) +
    geom_text(aes(label = percent(perc)), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = c("#E69F00", "#CC79A7", "#009E73"), name = "Changes")
  pie
  
  if (type == 'PG'){
    pie <- pie + labs(title = "Pie Chart of Changes between F9_D4_PG and F9_UD")
    
    df2 <- data.frame(group = c("Enhancers Not In F9_D4_PG + F9_UD", "Enhancers In F9_D4_PG + F9_UD"), 
                      value = c(3840 - nrow(subDat_fin), nrow(subDat_fin)))
    df2$perc <- df2$value/sum(df2$value)
    
    pie2 <- ggplot(data = df2, aes(x = "", y = value, fill = group)) + 
      geom_bar(stat = "identity") + coord_polar("y", start = 0) + blank_theme + 
      theme(axis.text.x = element_blank()) +
      scale_fill_discrete(name = "Changes") +
      geom_text(aes(label = percent(perc)), position = position_stack(vjust = 0.5), size = 4) + 
      labs(title = "Percentages of LSD1 Enhancers in both F9_D4_PG and F9_UD")
    
    pdf(paste("Figures/", figPreName, "_LSD1_PieChart_EnhNumbers.pdf", sep = ""), width = 8, height = 4)
    print(pie2)  
    dev.off()

  } 
  
  if (type == 'TCP'){
    pie <- pie + labs(title = "Pie Chart of Changes between F9_D4_TCP and F9_UD")
    
    df2 <- data.frame(group = c("Enhancers Not In F9_D4_TCP + UD", "Enhancers In F9_D4_TCP + F9_UD"), 
                      value = c(3840 - nrow(subDat_fin), nrow(subDat_fin)))
    df2$perc <- df2$value/sum(df2$value)
    
    pie2 <- ggplot(data = df2, aes(x = "", y = value, fill = group)) + 
      geom_bar(stat = "identity") + coord_polar("y", start = 0) + blank_theme + 
      theme(axis.text.x = element_blank()) +
      scale_fill_discrete(name = "Changes") +
      geom_text(aes(label = percent(perc)), position = position_stack(vjust = 0.5), size = 4) + 
      labs(title = "Percentages of LSD1 Enhancers in both F9_D4_TCP and F9_UD")
    
    pdf(paste("Figures/", figPreName, "_LSD1_PieChart_EnhNumbers.pdf", sep = ""), width = 8, height = 4)
    print(pie2)  
    dev.off()
    
  }
  
  if (type == 'None'){
    pie <- pie + labs(title = "Pie Chart of Changes between F9_D4 and F9_UD")
    
    df2 <- data.frame(group = c("Enhancers Not In F9_D4 + F9_UD", "Enhancers In F9_D4 + F9_UD"), 
                      value = c(3840 - nrow(subDat_fin), nrow(subDat_fin)))
    df2$perc <- df2$value/sum(df2$value)
    
    pie2 <- ggplot(data = df2, aes(x = "", y = value, fill = group)) + 
      geom_bar(stat = "identity") + coord_polar("y", start = 0) + blank_theme + 
      theme(axis.text.x = element_blank()) +
      scale_fill_discrete(name = "Changes") +
      geom_text(aes(label = percent(perc)), position = position_stack(vjust = 0.5), size = 4) + 
      labs(title = "Percentages of LSD1 Enhancers in both F9_D4 and F9_UD")
    
    pdf(paste("Figures/", figPreName, "_LSD1_PieChart_EnhNumbers.pdf", sep = ""), width = 8, height = 4)
    print(pie2)  
    dev.off()
  }
    
  pdf(paste("Figures/", figPreName, "_LSD1_PieChart_Changes.pdf", sep = ""), width = 8, height = 4)
  print(pie)  
  dev.off()
  
}

# *********************************************************************************** #
dat <- read.csv(file = "FC_allSamples_LSD1_counts.tsv", header = T, sep = "\t")

# Subset data to only treatments
subDat <- dat[, c('X', 'F9_UD', 'F9_D4', 'F9_D4_Min_UD', 'F9_D4_Min_UD_PerChange')]
sameInd <- which(subDat$F9_UD == 0.0001 & subDat$F9_D4 == 0.0001)
subDat_fin <- subDat[-sameInd, ]
colnames(subDat_fin)[1] <- 'GeneName'
nrow(subDat_fin)

subDat_fin <- subDat_fin[order(subDat_fin$F9_D4_Min_UD, decreasing = T), ]

CreateChangePlot(subDat_fin, figPreName = "F9_D4_Min_UD", type = "None")

# write.table(subDat_fin, file = "Figures/D4_Min_UD_LSD1_Changes.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# *********************************************************************************** #
subDat <- dat[, c('X', 'F9_UD', 'F9_D4_PG', 'F9_D4_PG_Min_UD', 'F9_D4_PG_Min_UD_PerChange')]
sameInd <- which(subDat$F9_UD == 0.0001 & subDat$F9_D4_PG == 0.0001)
subDat_fin <- subDat[-sameInd, ]
colnames(subDat_fin)[1] <- 'GeneName'
nrow(subDat_fin)

subDat_fin <- subDat_fin[order(subDat_fin$F9_D4_PG_Min_UD, decreasing = T), ]

CreateChangePlot(subDat_fin, figPreName = "F9_D4_PG_Min_UD", type = "PG")

# write.table(subDat_fin, file = "Figures/D4_PG_Min_UD_LSD1_Changes.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# *********************************************************************************** #
subDat <- dat[, c('X', 'F9_UD', 'F9_D4_TCP', 'F9_D4_TCP_Min_UD', 'F9_D4_TCP_Min_UD_PerChange')]
sameInd <- which(subDat$F9_UD == 0.0001 & subDat$F9_D4_TCP == 0.0001)
subDat_fin <- subDat[-sameInd, ]
colnames(subDat_fin)[1] <- 'GeneName'
nrow(subDat_fin)

subDat_fin <- subDat_fin[order(subDat_fin$F9_D4_TCP_Min_UD, decreasing = T), ]
nrow(subDat_fin)

CreateChangePlot(subDat_fin, figPreName = "F9_D4_TCP_Min_UD", type = "TCP")

# write.table(subDat_fin, file = "Figures/D4_TCP_Min_UD_LSD1_Changes.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# *************** Find length of LSD1 enhancers  for each treatment **************** # 
totFound <- nrow(dat)
ud <- dat[dat$F9_UD != 0.0001,]
nrow(ud)

D4 <- dat[dat$F9_D4 != 0.0001,]
nrow(D4)

D4_pg <- dat[dat$F9_D4_PG != 0.0001,]
nrow(D4_pg)

D4_tcp <- dat[dat$F9_D4_TCP != 0.0001, ]
nrow(D4_tcp)

# *************** Plot percent change for each treatment **************** # 
dat <- read.csv(file = "FC_allSamples_LSD1_counts.tsv", header = T, sep = "\t")

datChange <- dat[, 6:8]

datChange <- melt(datChange)
datChange$variable <- as.factor(datChange$variable)
datChange$pos <- ifelse(datChange$value > 0, "pos", "neg")

p <- ggplot(datChange, aes(x = value, color = variable)) + geom_density() + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  labs(title = "Density plot of Treatment minus UD") + xlab("Changes") + 
  ylab("Density")
p

pdf(paste("Figures/F9_LSD1_AllDensityChanges.pdf", sep = ""), width = 8, height = 4)
print(p)  
dev.off()

datVal <- dat[, 2:5]
datVal <- melt(datVal)
datVal$variable <- as.factor(datVal$variable)
datVal$pos <- ifelse(datVal$value > 0, "pos", "neg")

p <- ggplot(datVal, aes(x = value, color = variable)) + geom_density() + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  labs(title = "Density plot of Treatment minus UD") + xlab("Changes") + 
  ylab("Density")
p

# *********************************** Create Plots of Changes for ESC-J1 ********************************************** #
CreateChangePlot_ESC <- function(subDat_fin, figPreName, type){
  numTot <- 87603
  colnames(subDat_fin) <- c('GeneName', 'UD', 'Trt', 'Trt_Min_UD', 'Trt_Min_UD_PerChange')
  
  subDat_fin$Pos <- ifelse(subDat_fin$Trt_Min_UD > 0, TRUE, FALSE)
  subDat_fin <- subDat_fin[order(subDat_fin$Trt_Min_UD, decreasing = T), ]
  
  subDat_fin$GeneName <- factor(subDat_fin$GeneName, levels = subDat_fin$GeneName[order(subDat_fin$Trt_Min_UD, decreasing = T)])
  subDat_fin$Pos <- as.factor(subDat_fin$Pos)
  
  # Waterfall Plot
  b <- ggplot(subDat_fin, aes(x = GeneName, y = Trt_Min_UD, color = Pos)) +
    scale_fill_discrete(name = "Positive Changes") + 
    theme_classic() + 
    theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.y = element_text(face="bold", angle = 90), legend.position = "none") 
  b <- b + geom_bar(stat = 'identity', position = position_dodge(width = 0.4)) +
    scale_color_manual(values =  c("#E69F00", "#009E73"))
  
  if (type == 'PG'){
    b <- b + labs(title = "Waterfall plot for changes in Methylation between F9_D4_PG and F9_UD", x = NULL, y = "Quantity of Change") +
      geom_hline(yintercept = c(100, -100), linetype = "dashed", color = "red")
  } 
  
  if (type == 'TCP'){
    b <- b + labs(title = "Waterfall plot for changes in Methylation between F9_D4_TCP and F9_UD", x = NULL, y = "Quantity of Change") +
      geom_hline(yintercept = c(100, -100), linetype = "dashed", color = "red")
  }
  
  if (type == 'None'){
    b <- b + labs(title = "Waterfall plot for changes in Methylation between F9_D4 and F9_UD", x = NULL, y = "Quantity of Change") +
      geom_hline(yintercept = c(100, -100), linetype = "dashed", color = "red")
  }
  
  pdf(paste("Figures/", figPreName, "_ESC_J1_waterfallPlot.pdf", sep = ""), width = 8, height = 6)
  print(b)  
  dev.off()
  
  # PieChart
  pos <- sum(subDat_fin$Trt_Min_UD > 0)
  neg <- sum(subDat_fin$Trt_Min_UD < 0)
  nch <- sum(subDat_fin$Trt_Min_UD == 0)
  df <- data.frame(group = c("Positive", "NoChange", "Negative"),
                   value = c(pos, nch, neg))
  df$perc <- df$value/sum(df$value)
  
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  pie <- ggplot(data = df, aes(x = "", y = value, fill = group)) + 
    geom_bar(stat = "identity") + coord_polar("y", start = 0) + blank_theme + 
    theme(axis.text.x = element_blank()) +
    geom_text(aes(label = percent(perc)), position = position_stack(vjust = 0.5), size = 4) +
    scale_fill_manual(values = c("#E69F00", "#CC79A7", "#009E73"), name = "Changes")
  pie
  
  if (type == 'PG'){
    pie <- pie + labs(title = "Pie Chart of Changes between F9_D4_PG and F9_UD")
    
    df2 <- data.frame(group = c("Enhancers Not In F9_D4_PG + F9_UD", "Enhancers In F9_D4_PG + F9_UD"), 
                      value = c(numTot - nrow(subDat_fin), nrow(subDat_fin)))
    df2$perc <- df2$value/sum(df2$value)
    
    pie2 <- ggplot(data = df2, aes(x = "", y = value, fill = group)) + 
      geom_bar(stat = "identity") + coord_polar("y", start = 0) + blank_theme + 
      theme(axis.text.x = element_blank()) +
      scale_fill_discrete(name = "Changes") +
      geom_text(aes(label = percent(perc)), position = position_stack(vjust = 0.5), size = 4) + 
      labs(title = "Percentages of ESC_J1 Enhancers in both F9_D4_PG and F9_UD")
    
    pdf(paste("Figures/", figPreName, "_ESC_J1_PieChart_EnhNumbers.pdf", sep = ""), width = 8, height = 4)
    print(pie2)  
    dev.off()
    
  } 
  
  if (type == 'TCP'){
    pie <- pie + labs(title = "Pie Chart of Changes between F9_D4_TCP and F9_UD")
    
    df2 <- data.frame(group = c("Enhancers Not In F9_D4_TCP + UD", "Enhancers In F9_D4_TCP + F9_UD"), 
                      value = c(numTot - nrow(subDat_fin), nrow(subDat_fin)))
    df2$perc <- df2$value/sum(df2$value)
    
    pie2 <- ggplot(data = df2, aes(x = "", y = value, fill = group)) + 
      geom_bar(stat = "identity") + coord_polar("y", start = 0) + blank_theme + 
      theme(axis.text.x = element_blank()) +
      scale_fill_discrete(name = "Changes") +
      geom_text(aes(label = percent(perc)), position = position_stack(vjust = 0.5), size = 4) + 
      labs(title = "Percentages of ESC_J1 Enhancers in both F9_D4_TCP and F9_UD")
    
    pdf(paste("Figures/", figPreName, "_ESC_J1_PieChart_EnhNumbers.pdf", sep = ""), width = 8, height = 4)
    print(pie2)  
    dev.off()
    
  }
  
  if (type == 'None'){
    pie <- pie + labs(title = "Pie Chart of Changes between F9_D4 and F9_UD")
    
    df2 <- data.frame(group = c("Enhancers Not In F9_D4 + F9_UD", "Enhancers In F9_D4 + F9_UD"), 
                      value = c(numTot - nrow(subDat_fin), nrow(subDat_fin)))
    df2$perc <- df2$value/sum(df2$value)
    
    pie2 <- ggplot(data = df2, aes(x = "", y = value, fill = group)) + 
      geom_bar(stat = "identity") + coord_polar("y", start = 0) + blank_theme + 
      theme(axis.text.x = element_blank()) +
      scale_fill_discrete(name = "Changes") +
      geom_text(aes(label = percent(perc)), position = position_stack(vjust = 0.5), size = 4) + 
      labs(title = "Percentages of ESC_J1 Enhancers in both F9_D4 and F9_UD")
    
    pdf(paste("Figures/", figPreName, "_ESC_J1_PieChart_EnhNumbers.pdf", sep = ""), width = 8, height = 4)
    print(pie2)  
    dev.off()
  }
  
  pdf(paste("Figures/", figPreName, "_ESC_J1_PieChart_Changes.pdf", sep = ""), width = 8, height = 4)
  print(pie)  
  dev.off()
  
}

# ********************************************************************************** #

dat <- read.csv(file = "FC_allSamples_ESC_J1_counts.tsv", header = T, sep = "\t")

# Subset data to only treatments
subDat <- dat[, c('X', 'F9_UD', 'F9_D4', 'F9_D4_Min_UD', 'F9_D4_Min_UD_PerChange')]
sameInd <- which(subDat$F9_UD == 0.0001 & subDat$F9_D4 == 0.0001)
subDat_fin <- subDat[-sameInd, ]
colnames(subDat_fin)[1] <- 'GeneName'
nrow(subDat_fin)

subDat_fin <- subDat_fin[order(subDat_fin$F9_D4_Min_UD, decreasing = T), ]

CreateChangePlot_ESC(subDat_fin, figPreName = "F9_D4_Min_UD", type = "None")

# write.table(subDat_fin, file = "Figures/D4_Min_UD_LSD1_Changes.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# *********************************************************************************** #
subDat <- dat[, c('X', 'F9_UD', 'F9_D4_PG', 'F9_D4_PG_Min_UD', 'F9_D4_PG_Min_UD_PerChange')]
sameInd <- which(subDat$F9_UD == 0.0001 & subDat$F9_D4_PG == 0.0001)
subDat_fin <- subDat[-sameInd, ]
colnames(subDat_fin)[1] <- 'GeneName'
nrow(subDat_fin)

subDat_fin <- subDat_fin[order(subDat_fin$F9_D4_PG_Min_UD, decreasing = T), ]

CreateChangePlot_ESC(subDat_fin, figPreName = "F9_D4_PG_Min_UD", type = "PG")

# write.table(subDat_fin, file = "Figures/D4_PG_Min_UD_LSD1_Changes.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# *********************************************************************************** #
subDat <- dat[, c('X', 'F9_UD', 'F9_D4_TCP', 'F9_D4_TCP_Min_UD', 'F9_D4_TCP_Min_UD_PerChange')]
sameInd <- which(subDat$F9_UD == 0.0001 & subDat$F9_D4_TCP == 0.0001)
subDat_fin <- subDat[-sameInd, ]
colnames(subDat_fin)[1] <- 'GeneName'
nrow(subDat_fin)

subDat_fin <- subDat_fin[order(subDat_fin$F9_D4_TCP_Min_UD, decreasing = T), ]
nrow(subDat_fin)

CreateChangePlot_ESC(subDat_fin, figPreName = "F9_D4_TCP_Min_UD", type = "TCP")

# write.table(subDat_fin, file = "Figures/D4_TCP_Min_UD_LSD1_Changes.tsv", sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# *************** Find length of LSD1 enhancers  for each treatment **************** # 
totFound <- nrow(dat)
ud <- dat[dat$F9_UD != 0.0001,]
nrow(ud)

D4 <- dat[dat$F9_D4 != 0.0001,]
nrow(D4)

D4_pg <- dat[dat$F9_D4_PG != 0.0001,]
nrow(D4_pg)

D4_tcp <- dat[dat$F9_D4_TCP != 0.0001, ]
nrow(D4_tcp)

# *************** Plot percent change for each treatment **************** # 
datChange <- dat[, 6:8]

datChange <- melt(datChange)
datChange$variable <- as.factor(datChange$variable)
datChange$pos <- ifelse(datChange$value > 0, "pos", "neg")

p <- ggplot(datChange, aes(x = value, color = variable)) + geom_density() + 
  theme(legend.position = "bottom", legend.title = element_blank()) + 
  labs(title = "Density plot of Treatment minus UD") + xlab("Changes") + 
  ylab("Density")
p

pdf(paste("Figures/F9_ESC_J1_AllDensityChanges.pdf", sep = ""), width = 8, height = 4)
print(p)  
dev.off()
