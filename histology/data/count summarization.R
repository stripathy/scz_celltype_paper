library(ggplot2)

setwd("C:/Users/Dwight_Newton/Dropbox/PhD/RNAseq - Cell Specific/Aim 1 (Pitt Cohort)/Cell counting")
options(stringsAsFactors = FALSE)

groups <- read.csv("groups.csv")

counts <- read.csv("cell count subset.csv")
countSummary <- as.data.frame(matrix(nrow=19, ncol=9))
names(countSummary) <- c("Subject", "Group", "Area(um^2)", "PYR_23", "PYR_56", "SST", "PV", "VIP", "PYR_total")

seqindex <- seq(from=1, to=nrow(counts), by=40)

#Error w/ subject 683, need to figure that out
for (i in 4:19){
  j <- seqindex[i]
  #Get subject ID, as numeric
  countSummary[i,1] <- as.numeric(substr(counts$Subject[j], 1, (regexpr("-", counts$Subject[j])-2)))
  countSummary[i,2] <- groups$Subject.Group[groups$HU. == countSummary[i,1]]
  #Area
  countSummary[i,3] <- mean(counts$Area..um.2.[c(j, j+20)], na.rm = TRUE)
  #L2/3 PYR
  countSummary[i,4] <- mean(counts$X488.Channel[c(j:(j+9))], na.rm = TRUE)
  #L5/6 PYR
  countSummary[i,5] <- mean(counts$X488.Channel[c((j+10):(j+19))], na.rm = TRUE)
  #SST
  countSummary[i,6] <- mean(counts$X488.Channel[c((j+20):(j+40))], na.rm = TRUE)
  #PV
  countSummary[i,7] <- mean(counts$X568.Channel[c(j:(j+19))], na.rm = TRUE)
  #VIP
  countSummary[i,8] <- mean(counts$X568.Channel[c((j+20):(j+40))], na.rm = TRUE)
  #TotalPYR
  countSummary[i,9] <- mean(counts$X488.Channel[c(j:(j+20))], na.rm = TRUE)
}

#assume 683 was mis-labelled from 863 (MDD)
for (i in 3){
  j <- seqindex[i]
  #Get subject ID, as numeric
  countSummary[i,1] <- as.numeric(substr(counts$Subject[j], 1, (regexpr("-", counts$Subject[j])-2)))
  #countSummary[i,2] <- groups$Subject.Group[groups$HU. == countSummary[i,1]]
  #Area
  countSummary[i,3] <- mean(counts$Area..um.2.[c(j, j+20)], na.rm = TRUE)
  #L2/3 PYR
  countSummary[i,4] <- mean(counts$X488.Channel[c(j:(j+9))], na.rm = TRUE)
  #L5/6 PYR
  countSummary[i,5] <- mean(counts$X488.Channel[c((j+10):(j+19))], na.rm = TRUE)
  #SST
  countSummary[i,6] <- mean(counts$X488.Channel[c((j+20):(j+40))], na.rm = TRUE)
  #PV
  countSummary[i,7] <- mean(counts$X568.Channel[c(j:(j+19))], na.rm = TRUE)
  #VIP
  countSummary[i,8] <- mean(counts$X568.Channel[c((j+20):(j+40))], na.rm = TRUE)
  #TotalPYR
  countSummary[i,9] <- mean(counts$X488.Channel[c(j:(j+20))], na.rm = TRUE)
}
countSummary[3,2] <- "MDD"

write.csv(countSummary, "preliminary counts.csv")

#main effect of group?
Areares <- aov(`Area(um^2)` ~ Group, data=countSummary)
PYR23res <- aov(PYR_23 ~ Group, data=countSummary)
PYR56res <- aov(PYR_56 ~ Group, data=countSummary)
SSTres <- aov(SST ~ Group, data=countSummary)
PVres <- aov(PV ~ Group, data=countSummary)
VIPres <- aov(VIP ~ Group, data=countSummary)
AllPyrres <- aov(PYR_total ~ Group, data=countSummary)

summary(Areares)
summary(PYR23res)
summary(PYR56res)
summary(SSTres)
summary(PVres)
summary(VIPres)
summary(AllPyrres)

#All very n.s. except for L5/6 PYR cells -   SCHIZ had lower density v.s. controls and MDD
TukeyHSD(PYR56res)
TukeyHSD(AllPyrres)

#PYR figures
#weird value inserted, removed here
countSummary$PYR_total[8] <- NA

countSummary$Group <- factor(countSummary$Group, levels = c("Control", "MDD", "Bipolar", "SCHIZ"))
countSummary_Means <- aggregate(. ~ Group, mean, data=countSummary)

ggplot(countSummary, aes(x=Group, y=PYR_23, colour=Group)) +
  geom_jitter(width = 0.15, size=3) +
  geom_crossbar(data=countSummary_Means, aes(ymin = PYR_23, ymax = PYR_23),size=0.5, col="black", width = 0.3) +
  scale_y_continuous(expand=c(0,0), limits=c(0,25)) +
  ggtitle("L2/3 PYR Cells") +
  ylab(label = bquote("Cells / 100,000"~ mu*"M"^2)) +
  theme(axis.title.x = element_blank(), legend.position = "none")
  
ggplot(countSummary, aes(x=Group, y=PYR_56, colour=Group)) +
  geom_jitter(width = 0.15, size=3) +
  geom_crossbar(data=countSummary_Means, aes(ymin = PYR_56, ymax = PYR_56),size=0.5, col="black", width = 0.3) +
  scale_y_continuous(expand=c(0,0), limits=c(0,25)) +
  ggtitle("L5/6 PYR Cells") +
  ylab(label = bquote("Cells / 100,000"~ mu*"M"^2)) +
  theme(axis.title.x = element_blank(), legend.position = "none")

ggplot(countSummary, aes(x=Group, y=PYR_total, colour=Group)) +
  geom_jitter(width = 0.15, size=3) +
  geom_crossbar(data=countSummary_Means, aes(ymin = PYR_total, ymax = PYR_total),size=0.5, col="black", width = 0.3) +
  scale_y_continuous(expand=c(0,0), limits=c(0,25)) +
  ggtitle("Total PYR Cells") +
  ylab(label = bquote("Cells / 100,000"~ mu*"M"^2)) +
  theme(axis.title.x = element_blank())



#Interneurons
ggplot(countSummary, aes(x=Group, y=SST, colour=Group)) +
  geom_jitter(width = 0.15, size=3) +
  geom_crossbar(data=countSummary_Means, aes(ymin = SST, ymax = SST),size=0.5, col="black", width = 0.3) +
  scale_y_continuous(expand=c(0,0), limits=c(0,10)) +
  ggtitle("SST Cells") +
  ylab(label = bquote("Cells / 100,000"~ mu*"M"^2)) +
  theme(axis.title.x = element_blank(), legend.position = "none")

ggplot(countSummary, aes(x=Group, y=PV, colour=Group)) +
  geom_jitter(width = 0.15, size=3) +
  geom_crossbar(data=countSummary_Means, aes(ymin = PV, ymax = PV),size=0.5, col="black", width = 0.3) +
  scale_y_continuous(expand=c(0,0), limits=c(0,2)) +
  ggtitle("PV Cells") +
  ylab(label = bquote("Cells / 100,000"~ mu*"M"^2)) +
  theme(axis.title.x = element_blank(), legend.position = "none")

ggplot(countSummary, aes(x=Group, y=VIP, colour=Group)) +
  geom_jitter(width = 0.15, size=3) +
  geom_crossbar(data=countSummary_Means, aes(ymin = VIP, ymax = VIP),size=0.5, col="black", width = 0.3) +
  scale_y_continuous(expand=c(0,0), limits=c(0,2)) +
  ggtitle("VIP Cells") +
  ylab(label = bquote("Cells / 100,000"~ mu*"M"^2)) +
  theme(axis.title.x = element_blank())



#Number of cells/slide
library(reshape2)
library(dplyr)
CellNums <- countSummary
#333x333uM is the FOV area at 20x
CellNums[,4:9] <- (CellNums[,4:9]/333^2)*CellNums$`Area(um^2)`

Cells <- melt(CellNums[,c(2,4:9)])

#PYRs
ggplot(filter(Cells, variable=="PYR_23"|variable=="PYR_56"|variable=="PYR_total"), aes(x=variable, y=value, colour=Group)) +
  geom_jitter(width=0.2, size=3) +
  geom_hline(yintercept=125) +
  ggtitle(label="PYR-Cells") +
  scale_x_discrete(labels=c("PYR_23"="L2/3", "PYR_56"="L5/6", "PYR_total"="Total"))+
  ylab("Cells / section") +
  theme(axis.title.x = element_blank())

#Interneurons
ggplot(filter(Cells, variable=="SST"|variable=="PV"|variable=="VIP"), aes(x=variable, y=value, colour=Group)) +
  geom_jitter(width=0.2, size=3) +
  geom_hline(yintercept=125) +
  ggtitle(label="Interneurons") +
  theme(axis.title.x = element_blank())
