birthRates = as.data.frame(c(-2.51405, 1.31437, -1.79009, 0.655826, -1.84878, 
                             0.91377, -1.21294, 0.194515, -2.49839, 0.888437, -1.39956, 
                             0.162913, -2.07355, 0.447807, -2.59514, 0.991364, -3.31508, 
                             1.04925))

names(birthRates)<-"b"

birthRates$type<-c("Head and neck",
                   "Head and neck","Esophogaeal",
                   "Esophogaeal","Colorectal",
                   "Colorectal","Rectal",
                   "Rectal","Breast",
                   "Breast","Cervix",
                   "Cervix","Melanoma",
                   "Melanoma","SC Lung",
                   "SC Lung",
                   "D-LBCL",
                   "D-LBCL")

growthRates = as.data.frame(c(-5.3441, 0.874367, -1.98498, 0.477869, -7.00215, 
                              1.54467, -5.46545, 0.578148, -6.30794, 0.998794, -6.78071, 
                              1.06405, -5.04616, 0.410637, -4.85863, 0.246221, -5.4622, 
                              1.83844))

names(growthRates)<-"g"

growthRates$type<-birthRates$type

inVec<-c(2,4,6,8,10,12,14,16,18)
gRLN<-list()
for (i in inVec){
  gRLN[[growthRates$type[i]]]<-rlnorm(100, meanlog = growthRates$g[i-1], sdlog =  growthRates$g[i])
}

bRLN<-list()
for (i in inVec){
  bRLN[[birthRates$type[i]]]<-rlnorm(100, meanlog = birthRates$b[i-1], sdlog =  birthRates$b[i])
}

dOVERb<-list()
for (name in names(bRLN)){
  dOVERb[[name]]<-as.data.frame(1-(gRLN[[name]]/bRLN[[name]]))
}

for (name in names(dOVERb)){
  dOVERb[[name]][dOVERb[[name]]<0] <- 0
}

typeList<-list()
for (name in names(dOVERb)){
  typeList[[name]]<-as.data.frame(repmat(name,100,1))
}
type<-unlist(typeList)

dOVERbUNLIST<-as.data.frame(unlist(dOVERb))
dOVERbUNLIST$type<-type

names(dOVERbUNLIST)<-c("db", "type")

# read in lagging chromosome data
inp<-read.table("../data/Predicted_LaggingChrPct_TCGA.txt")
names(inp)<-inp[1,]
inp<-inp[-c(1),]
inp$Lagging.Chromosome<-as.numeric(inp$Lagging.Chromosome)
inp$Lagging.Chromosome[inp$Lagging.Chromosome>100]=100
inp$Lagging.Chromosome<-inp$Lagging.Chromosome/100
## change names in inp to match type
inp$group<-gsub("BRCA", "Breast", inp$group)
inp$group<-gsub("CESC", "Cervix", inp$group)
inp$group<-gsub("COAD", "Colorectal", inp$group)
inp$group<-gsub("DLBC", "D-LBCL", inp$group)
inp$group<-gsub("ESCA", "Esophogeal", inp$group)
inp$group<-gsub("HNSC", "Head and neck", inp$group)
inp$group<-gsub("SKCM", "Melanoma", inp$group)
inp$group<-gsub("READ", "Rectal", inp$group)
inp$group<-gsub("LUSC", "SC Lung", inp$group)
## remove mets
inp<-inp[inp$group!="met",]

dOVERbUNLIST$AdjDB<-(1-dOVERbUNLIST$db)

plot.x <- ggplot(dOVERbUNLIST) + geom_boxplot(aes(type, AdjDB)) + scale_y_continuous(trans="log")
plot.y <- ggplot(inp) + geom_boxplot(aes(group, Lagging.Chromosome))

grid.arrange(plot.x, plot.y, ncol=2) # visual verification of the boxplots

plot.x <- layer_data(plot.x)[,1:6]
plot.y <- layer_data(plot.y)[,1:6]
colnames(plot.x) <- paste0("x.", gsub("y", "", colnames(plot.x)))
colnames(plot.y) <- paste0("y.", gsub("y", "", colnames(plot.y)))
df <- cbind(plot.x, plot.y); rm(plot.x, plot.y)
df$category <- sort(unique(dOVERbUNLIST$type))

df.outliers <- df %>%
  dplyr::select(category, x.middle, x.outliers, y.middle, y.outliers) %>%
  data.table::data.table()
df.outliers <- df.outliers[, list(x.outliers = unlist(x.outliers), y.outliers = unlist(y.outliers)), 
                           by = list(category, x.middle, y.middle)]

ggplot(df, aes(fill = category, color = category)) +
  
  # 2D box defined by the Q1 & Q3 values in each dimension, with outline
  geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), alpha = 0.3) +
  geom_rect(aes(xmin = x.lower, xmax = x.upper, ymin = y.lower, ymax = y.upper), 
            color = "black", fill = NA) +
  
  # whiskers for x-axis dimension with ends
  geom_segment(aes(x = x.min, y = y.middle, xend = x.max, yend = y.middle)) + #whiskers
  geom_segment(aes(x = x.min, y = y.lower, xend = x.min, yend = y.upper)) + #lower end
  geom_segment(aes(x = x.max, y = y.lower, xend = x.max, yend = y.upper)) + #upper end
  
  # whiskers for y-axis dimension with ends
  geom_segment(aes(x = x.middle, y = y.min, xend = x.middle, yend = y.max)) + #whiskers
  geom_segment(aes(x = x.lower, y = y.min, xend = x.upper, yend = y.min)) + #lower end
  geom_segment(aes(x = x.lower, y = y.max, xend = x.upper, yend = y.max)) + #upper end
  
  # outliers
#  geom_point(data = df.outliers, aes(x = x.outliers, y = y.middle), size = 3, shape = 1) + # x-direction
#  geom_point(data = df.outliers, aes(x = x.middle, y = y.outliers), size = 3, shape = 1) + # y-direction
  
  xlab("AdjDB") + ylab("misSeg") +
  coord_cartesian(xlim = c(-9, 3), ylim = c(-.1, 1)) +
  theme_classic()
