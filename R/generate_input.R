x <- read.csv("C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia/6069 Slides and Data_cells_xy.txt",sep="\t")

bounds_x <- c(min(x$Centroid.X.µm),max(x$Centroid.X.µm))
bounds_y <- c(min(x$Centroid.Y.µm),max(x$Centroid.Y.µm))

L <- max(diff(bounds_x),diff(bounds_y))
offset <- min(bounds_x[1],bounds_x[2])

x <- read.csv("C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia/6069_Slides_and_Data__vessels_xy.txt",sep="\t")
x$Centroid.X.µm <- x$Centroid.X.µm-offset
x$Centroid.Y.µm <- x$Centroid.Y.µm-offset

write.table(x,"C:/Users/4473331/Documents/projects/ploidyEvolution-main/Julia/6069_vessels.txt",sep=",",quote=F, row.names = F)