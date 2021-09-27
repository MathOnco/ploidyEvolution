### 2DmisSegboxplots.R by Thomas Veith
### This script will produce 2D boxplots for any two sets of data you read into it
### However, it was originally intended to be used for Chromosomal Mis-segregation over 1-death/birth
### It owes a lot to the walk-through here: https://stackoverflow.com/questions/46068074/double-box-plots-in-ggplot2

### Read in birth rate medians and sdNorms to generate log norm distribution
birthRates = as.data.frame(
  c(
    -2.51405,
    1.31437,
    -1.79009,
    0.655826,
    -1.84878,
    0.91377,
    -1.21294,
    0.194515,
    -2.49839,
    0.888437,
    -1.39956,
    0.162913,
    -2.07355,
    0.447807,
    -2.59514,
    0.991364,
    -3.31508,
    1.04925
  )
)
### set names for birthRates dataframe
names(birthRates) <- "b"

### set cancer types
birthRates$type <- c(
  "Head and neck",
  "Head and neck",
  "Esophageal",
  "Esophageal",
  "Colorectal",
  "Colorectal",
  "Rectal",
  "Rectal",
  "Breast",
  "Breast",
  "Cervix",
  "Cervix",
  "Melanoma",
  "Melanoma",
  "SC Lung",
  "SC Lung",
  "D-LBCL",
  "D-LBCL"
)

### read in growth rate medians and sdNorms to generate log norm distribution
growthRates = as.data.frame(
  c(
    -5.3441,
    0.874367,
    -1.98498,
    0.477869,
    -7.00215,
    1.54467,
    -5.46545,
    0.578148,
    -6.30794,
    0.998794,
    -6.78071,
    1.06405,
    -5.04616,
    0.410637,
    -4.85863,
    0.246221,
    -5.4622,
    1.83844
  )
)
### set name for growthRates vector
names(growthRates) <- "g"

### add type to growthRates
growthRates$type <- birthRates$type

### create an index vector that will allow us to grab mean and sdlog two at a time (i.e., to grab the two variables we need by type)
inVec <- c(2, 4, 6, 8, 10, 12, 14, 16, 18)
### set up list of log norm distributed growth rate vectors
gRLN <- list()
for (i in inVec) {
  gRLN[[growthRates$type[i]]] <-
    rlnorm(100, meanlog = growthRates$g[i - 1], sdlog =  growthRates$g[i])
}
### set up list of log norm distributed birth rate vectors
bRLN <- list()
for (i in inVec) {
  bRLN[[birthRates$type[i]]] <-
    rlnorm(100, meanlog = birthRates$b[i - 1], sdlog =  birthRates$b[i])
}
### set up list of death/birth dataframes where d=1-g/b
dOVERb <- list()
for (name in names(bRLN)) {
  dOVERb[[name]] <- as.data.frame(1 - (gRLN[[name]] / bRLN[[name]]))
}
### replace negative values in d/b with 0
for (name in names(dOVERb)) {
  dOVERb[[name]][dOVERb[[name]] < 0] <- 0
}
### create vector to add type to d/b
typeList <- list()
for (name in names(dOVERb)) {
  typeList[[name]] <- as.data.frame(repmat(name, 100, 1))
}
type <- unlist(typeList)
### make dataframe for all d/b
dOVERbUNLIST <- as.data.frame(unlist(dOVERb))
dOVERbUNLIST$type <- type
### set names for d/b dataframe
names(dOVERbUNLIST) <- c("db", "type")

# read in lagging chromosome data
inp <- read.table("../data/Predicted_LaggingChrPct_TCGA.txt")
names(inp) <- inp[1, ]
inp <- inp[-c(1), ]
inp$Lagging.Chromosome <- as.numeric(inp$Lagging.Chromosome)
inp$Lagging.Chromosome[inp$Lagging.Chromosome > 100] = 100
inp$Lagging.Chromosome <- inp$Lagging.Chromosome / 100
## change names in inp to match type
inp$group <- gsub("BRCA", "Breast", inp$group)
inp$group <- gsub("CESC", "Cervix", inp$group)
inp$group <- gsub("COAD", "Colorectal", inp$group)
inp$group <- gsub("DLBC", "D-LBCL", inp$group)
inp$group <- gsub("ESCA", "Esophageal", inp$group)
inp$group <- gsub("HNSC", "Head and neck", inp$group)
inp$group <- gsub("SKCM", "Melanoma", inp$group)
inp$group <- gsub("READ", "Rectal", inp$group)
inp$group <- gsub("LUSC", "SC Lung", inp$group)
## remove mets
inp <- inp[inp$group != "met", ]

# inp<-(inp[inp$group==c("HNSC","ESCA"),])


### read in crit curves
critcurve40 = read.csv("~/Downloads/criticalcurve_n_40 (1).csv")
critcurve20 = read.csv("~/Downloads/criticalcurve_n_20 (1).csv")
critcurve10 = read.csv("~/Downloads/criticalcurve_n_10 (1).csv")
critcurve5 = read.csv("~/Downloads/criticalcurve_n_5 (1).csv")
### set names
names(critcurve40) <- c("x1", "y1")
names(critcurve20) <- c("x1", "y1")
names(critcurve10) <- c("x1", "y1")
names(critcurve5) <- c("x1", "y1")
### change x-values to 1-x for visualization later
critcurve40$x1 <- (1 - critcurve40$x1)
critcurve20$x1 <- (1 - critcurve20$x1)
critcurve10$x1 <- (1 - critcurve10$x1)
critcurve5$x1 <- (1 - critcurve5$x1)
### set any zero values to 10^-3
critcurve40$x1[critcurve40$x1 == 0] <- 10 ^ (-3)
critcurve20$x1[critcurve20$x1 == 0] <- 10 ^ (-3)
critcurve10$x1[critcurve10$x1 == 0] <- 10 ^ (-3)
critcurve5$x1[critcurve5$x1 == 0] <- 10 ^ (-3)
### have to add a category to the crit curves so it will work with ggplot down later
critcurve5$category <- "Head and neck"
critcurve10$category = "Head and neck"
critcurve20$category = "Head and neck"
critcurve40$category = "Head and neck"
### Add 1-d/b to dOVERb dataframe (Again, for visualization)
dOVERbUNLIST$one_minus_d_over_b <- (1 - dOVERbUNLIST$db)

##########################
### VISUALIZATION TIME ###
#########################

### Okay, to start, you gotta make two box plots (For your two dimensions later)
plot.x <-
  ggplot(dOVERbUNLIST) + geom_boxplot(aes(type, one_minus_d_over_b))  + coord_flip()
plot.y <-
  ggplot(inp) + geom_boxplot(aes(group, Lagging.Chromosome))  + coord_flip()
#plot.z <- ggplot(dOVERbUNLIST) + geom_boxplot(aes(type, db)) + coord_flip()#

### Now take a look at the two boxplots together (good for checking that your 2D output is right later)
# pdf("2D_boxplots_setting_xyticks.pdf", width = 7,height = 7)
grid.arrange(plot.x, plot.y, ncol = 2) # visual verification of the boxplots
# dev.off()

### We are going to "layer the data" (i.e., take out what we want from the boxplots so we can construct our geom_rect objects later)
plot.x <- layer_data(plot.x)[, 1:6]
plot.y <- layer_data(plot.y)[, 1:6]
colnames(plot.x) <- paste0("x.", gsub("y", "", colnames(plot.x)))
colnames(plot.y) <- paste0("y.", gsub("y", "", colnames(plot.y)))

### Set up df dataframe that will be used by ggplot
df <- cbind(plot.x, plot.y)
rm(plot.x, plot.y)
df$category <- sort(unique(dOVERbUNLIST$type))
df.outliers <- df %>%
  dplyr::select(category, x.middle, x.outliers, y.middle, y.outliers) %>%
  data.table::data.table()
df.outliers <-
  df.outliers[, list(x.outliers = unlist(x.outliers),
                     y.outliers = unlist(y.outliers)),
              by = list(category, x.middle, y.middle)]

### Now that your df is set up, limit the y-values of your crit curves so that they are in the range of your df values (otherwise plot is ugly)
critcurve5$y1[critcurve5$y1 < min(df$y.min)] <- min(df$y.min)
critcurve10$y1[critcurve10$y1 < min(df$y.min)] <- min(df$y.min)
critcurve20$y1[critcurve20$y1 < min(df$y.min)] <- min(df$y.min)
critcurve40$y1[critcurve40$y1 < min(df$y.min)] <- min(df$y.min)

### And here the magic happens
ggplot(df, aes(fill = category, color = category)) +
  geom_line(data = critcurve5,
            aes(x = x1, y = y1),
            size = 1,
            color = "black") + geom_line(
              data = critcurve10,
              aes(x = x1, y = y1),
              size = 1,
              color = "black",
              linetype = "dotted"
            ) + geom_line(
              data = critcurve20,
              aes(x = x1, y = y1),
              size = 1,
              color = "black",
              linetype = "dashed"
            ) + geom_line(
              data = critcurve40,
              aes(x = x1, y = y1),
              size =
                1,
              color = "black",
              linetype = "dotdash"
            ) +
  
  # 2D box defined by the Q1 & Q3 values in each dimension, with outline
  geom_rect(aes(
    xmin = x.lower,
    xmax = x.upper,
    ymin = y.lower,
    ymax = y.upper
  ),
  alpha = 0.3) +
  geom_rect(
    aes(
      xmin = x.lower,
      xmax = x.upper,
      ymin = y.lower,
      ymax = y.upper
    ),
    color = "black",
    fill = NA
  ) +
  
  # whiskers for x-axis dimension with ends
  geom_segment(aes(
    x = x.min,
    y = y.middle,
    xend = x.max,
    yend = y.middle
  )) + #whiskers
  geom_segment(aes(
    x = x.min,
    y = y.lower,
    xend = x.min,
    yend = y.upper
  )) + #lower end
  geom_segment(aes(
    x = x.max,
    y = y.lower,
    xend = x.max,
    yend = y.upper
  )) + #upper end
  
  # whiskers for y-axis dimension with ends
  geom_segment(aes(
    x = x.middle,
    y = y.min,
    xend = x.middle,
    yend = y.max
  )) + #whiskers
  geom_segment(aes(
    x = x.lower,
    y = y.min,
    xend = x.upper,
    yend = y.min
  )) + #lower end
  geom_segment(aes(
    x = x.lower,
    y = y.max,
    xend = x.upper,
    yend = y.max
  )) + #upper end
  
  # outliers
  # geom_point(data = df.outliers, aes(x = x.outliers, y = y.middle), size = 3, shape = 1) + # x-direction
  # geom_point(data = df.outliers, aes(x = x.middle, y = y.outliers), size = 3, shape = 1) + # y-direction
  
  xlab("1-d/b") + ylab("misSeg") +
  scale_x_continuous(trans = "log",
                     breaks = c(0.002478752, 0.018315638, 0.135335283, 1.00)) + scale_y_continuous(trans =
                                                                                                     "log",
                                                                                                   breaks = c(0.002478752, 0.018315638, 0.135335283, 1.00)) +
  # scale_x_continuous(limits=c(0, 2)) + scale_y_continuous(limits=c(0,1))+
  # coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  theme_classic() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(colour = "black", size = 4)
  )
