### COLOGNE WORKSHOP: GEOMETRIC MORPHOMETRICS FOR ARCHAEOLOGISTS ###
### AUTHOR: Dr. Christian Steven Hoggard (University of Southampton) ###

### SYSTEM INFORMATION ###
### R version 3.6.1 (2019-07-05) ###
### Platform: x86_64-w64-mingw32/x64 (64-bit) ###
### Running under: Windows 10 x64 (build 18362) ###

### ATTACHED BASE PACKAGES:
### [1] stats     graphics  grDevices utils  
### [6] methods   base

### EXERCISE 1: R ENVIRONMENT SETUP
### Download all files from the GitHub Repository: https://github.com/CSHoggard/-cologne_workshop
### Set the working directory using the "Session/Set Working Directory" function in RStudio or through the setwd() function
### Run the setup.R code through opening the file, or alternatively:

if(!require("Momocs")) install.packages('Momocs', repos='http://cran.us.r-project.org') 
if(!require("tidyverse")) install.packages('tidyverse', repos='http://cran.us.r-project.org')
if(!require("cowplot")) install.packages('cowplot', repos='http://cran.us.r-project.org')

install.packages("devtools")
devtools::install_github("geomorphR/geomorph", ref = "Develop", build_vignettes = TRUE)

install.packages("BiocManager")
library(BiocManager)
install("ggtree")

### EXERCISE 2: CREATING AND READING LANDMARK DATA #1
### Load Geomorph v.3.1.3 (Beta) in RStudio
### Import the skull_1.ply file into RStudio using the read.ply() function and label it "SK1"
### Use the buildtemplate() function to create a template for (23 LM / 200 SLM)
### Experiment with plotting landmarks for this or any other skull using digitsurface()
### Use the readmulti.nts() function to read in the already digitised files

library(geomorph)
library(Momocs)
library(ggtree)
library(tidyverse)
library(cowplot)

SK1 <- read.ply("skull_1.ply", ShowSpecimen = FALSE, addNormals = TRUE) ### import .PLY file
plot3d(SK1, col = "yellow") ### plots the 3d mesh model (colour is optional)
buildtemplate(SK1, 23, 200) ### build a template using the first .PLY file (23 landmarks, 200 surface semilandmarks)

### The specimen you choose for your template is your most complete and definable specimen
### The order of the landmarks you make on this template will need to be replicated onto the rest of your sample
### The semilandmarks that are then used in the template will be projected onto the rest of the sample

SK2 <- read.ply("skull_2.ply",ShowSpecimen = TRUE, addNormals = FALSE) ### import .PLY file
digitsurface(SK2, 23) ### digitising the surface (23 landmarks)

SK3 <- read.ply("skull_3.ply", ShowSpecimen = TRUE, addNormals = FALSE) ### import .PLY file
digitsurface(SK3, 23) ### digitising the surface (23 landmarks)

SK4 <- read.ply("skull_4.ply", ShowSpecimen = TRUE, addNormals = FALSE) ### import .PLY file
digitsurface(SK4, 23) ### digitising the surface (23 landmarks)

SK5 <- read.ply("skull_5.ply", ShowSpecimen = TRUE, addNormals = FALSE) ### import .PLY file
digitsurface(SK5, 23) ### digitising the surface (23 landmarks)

SK6 <- read.ply("skull_6.ply", ShowSpecimen = TRUE, addNormals = FALSE) ### import .PLY file
digitsurface(SK6, 23) ### digitising the surface (23 landmarks)

skull  <- readmulti.nts(c("SK1.nts","SK2.nts","SK3.nts","SK4.nts","SK5.nts","SK6.nts")) ### importing finished file

groups <- read.csv("skulls.csv", header=T, row.names=1) ### import the metadata
is.factor(groups$Sex) ### checks factor
is.factor(groups$Location) ### checks factor
is.character(groups$Code) ### checks character
groups$Code <- as.character(groups$Code) ### convert to character

surfslide<-read.csv("surfslide.csv", header=TRUE) ### sliding surface semilandmarks (configuration)
surfslide<-as.matrix(surfslide) ### convert to matrix

### EXERCISE 3: CREATING AND READING LANDMARK DATA #2
### Download tpsDig2 and tpsUtiltpsUtil: https://life.bio.sunysb.edu/morph/soft-utility.htmltpsDig2: https://life.bio.sunysb.edu/morph/soft-dataacq.html 
### Collate all the Ellensbanke and Eskebjerg tanged points into one .tps file in tpsUtil using the "Build tps file from images" - save the file as "tanged.tps"
### Load the new.tps file in TpsDig2 and place six landmarks at the extremities and shoulder mid-points
### Save the .tps file and import into R using the readland.tps function

tanged.points <- import_tps("tanged.tps") ### import the tps
tanged.data <- read.csv("tanged.csv", header = TRUE, row.names = 1) ### metadata
outline <- rbind(c(1,2),c(2,3),c(3,4),c(4,5),c(5,6),c(6,1)) ### outline linking
tanged.points <- Ldk(tanged.points$coo, fac = tanged.data, links = outline) ### convert to landmark file in Momocs (this procedure can also be done in Geomorph!)
inspect(tanged.points)

panel(tanged.points, names = TRUE) ### visualisation
panel(tanged.points, names = TRUE, fac = "Site", cex.names = 1, points.pch = 16, points.cex = 0.5) ### stylistic changes

library(GUImorph) ### load GUImorph
GUImorph() ### activate

### EXERCISE 4: GENERALISED PROCRUSTES ANALYSIS
### Using both the multi.nts file (skull) and the Ldk file (tanged.points) perform a Full Generalised Procrustes alignment for all examples.
### For the multi.nts file use the gpagen function in Geomorph (call ?gpagen)
### For the tanged.points file use the fgProcrustes function in Geomorph (call ?fgProcrustes)
### Explore visualisation styles for both gpa files (?stack) and Momocs (?plot)

gpatp <- fgProcrustes(tanged.points)
stack(gpatp, ldk_cex = 2,  title = "GPA: Tanged Points", meanshape = TRUE)

gpaskull <- gpagen(skull, Proj = TRUE, ProcD = TRUE, curves = NULL, surfaces = surfslide) ### generalised Procrustes Analysis
gpaskull ### calls the object
plot(gpaskull)
ref <- mshape(gpaskull$coords) ### calculates the mean shape of all skulls
ref ### calls the object
ref <- as.matrix(ref)	### converts the object to a matrix

### EXERCISE 5: VISUALISING SHAPE CHANGE
### Explore the tps_grid, tps_iso, tps_raw and tps_arr functions in Momocs to compare between two examples
### Note: these visualisations work better with outline data 

tps_grid(gpatp$coo$Ellensbanke_1, gpatp$coo$Eskebjerg_1_4, legend = FALSE, poly = FALSE, amp = 1, grid.size = 5)
tps_iso(gpatp$coo$Ellensbanke_1, gpatp$coo$Eskebjerg_1_4, shp = FALSE, amp = 1, grid.size = 5)
tps_raw(gpatp$coo$Ellensbanke_1, gpatp$coo$Eskebjerg_1_4, amp = 1, grid.size = 5)
tps_arr(gpatp$coo$Ellensbanke_1, gpatp$coo$Eskebjerg_1_5, shp = FALSE, amp = 1, grid.size = 5)

### EXERCISE 6: PRINCIPAL COMPONENT ANALYSIS
### Explore the main sources of theoretical shape variation and the distribution of points through the plotting of a principal component space for both datasets
### Produce a scree table and scree plot for the tanged points
### Customise where possible

pcatp <- PCA(gpatp, fac = tanged.data) ### creation of the principal component class item
pcatp$x ### principal component scores
database <- as_df(pcatp) ### convert to ggplot object (as a data frame)
ggplot(database, aes(PC1, PC2)) + geom_point() ### produces the ggplot object
scree(pcatp) ### scree table
scree_plot(pcatp) ### scree plot
plot(pcatp) ### plots the PCA graphic
plot(pcatp, xax = 2, yax = 5, cex = 1.5, fac = "Site", grid = FALSE, palette = col_autumn, chull.filled.alpha = 0.6, title = "Principal Component Analysis (PC1 vs. PC2)") ### customise the PCA graphic

pcasex <- plotTangentSpace(gpaskull$coords, axis1 = 1, axis2 = 2, warpgrids = TRUE, groups = groups$Sex, verbose = TRUE, label=groups$Code) ### principal component analysis
summary(pcasex) ### pca summary
pcasex$pc.shapes ### output (shape coordinates of the extreme ends of all PC axes)
pcasex$pc.shapes$PC1max ### e.g. PC1 max
pcasex$rotation ### rotation values

plotRefToTarget(ref, pcasex$pc.shapes$PC1max, method= "points") ### shape change from mean shape (ref) to max PC1 (as vector)

pcalocation <- plotTangentSpace(gpaskull$coords, axis1 = 1, axis2 = 2, warpgrids = TRUE, groups = groups$Location, verbose = TRUE, label=groups$Code) ### principal component analysis
summary(pcalocation) ### pca summary
pcalocation$pc.shapes ### output (shape coordinates of the extreme ends of all PC axes)
pcalocation$pc.shapes$PC1max ### e.g. PC1 max
pcalocation$rotation ### rotation values

### EXERCISE 7: DISCRIMINANT FUNCTION ANALYSIS
### Use the LDA() function in Momocs to create a discriminant analysis for the site variable
### Plot the resulting graphic and examine the LDA file to extract scores
### Customise the Lineal Discriminant plot

ldatp <- LDA(pcatp, fac = "Site") ### create the LDA file
ldatp ### call the LDA item
plot(ldatp)

### EXERCISE 8: MANOVA / Procrustes ANOVA
### Using the MANOVA() and procD.lm() functions test for difference between:
#### Tanged points vs. Site
#### Skull vs. sex
#### Skull vs. location
### For geomorph the data will need to be in a two-dimensional array format!

MANOVA(pcatp, fac = "Site", test = "Hotelling")

lmspecimensex <- procD.lm(two.d.array(gpaskull$coords) ~ groups$Sex, iter=99) ### Procrustes ANOVA (shape vs. sex)
lmspecimensex$aov.table ### anova table (summary)
lmspecimensex$call ### calls the code used
lmspecimensex$QR ### QR decompositions
lmspecimensex$fitted ### the fitted values
lmspecimensex$residuals ### the residuals (observed responses)
lmspecimensex$data ### the data frame for the model

lmspecimenlocation <- procD.lm(two.d.array(gpaskull$coords) ~ groups$Location, iter=99) ### Procrustes ANOVA (shape vs. location)
lmspecimenlocation$aov.table ### anova table (summary)
lmspecimenlocation$call ### calls the code used
lmspecimenlocation$QR ### QR decompositions
lmspecimenlocation$fitted ### the fitted values
lmspecimenlocation$residuals ### the residuals (observed responses)
lmspecimenlocation$data ### the data frame for the model

### EXERCISE 9: MEAN AND MEDIAN SHAPES 
### Produce the mean and median shapes for the tanged points (according to site)

tpms <- mshapes(gpatp, FUN = mean, "Site") ### change FUN for median
panel(Ldk(tpms$Coe), names = TRUE, cex.names = 0.75, points.pch = 16, points.cex = 1.5)

### EXERCISE 10: CLUSTER ANALYSIS
### Using the CLUST() function produce a phylo class object with the following conditions:
#### Distance method: "Euclidean"
#### Hclust method: "complete"
#### Site as a factor
### Visualise different tree types (cladogram, phylogram and fan!) and play the visual tools
### Optional: install ggtree and explore the phylo class object further

cluster <- CLUST(pcatp, fac = "Site", type = "phylogram", dist_method = "euclidean", hclust_method = "complete") ### hierarchical clustering (complete)
cluster <- CLUST(pcatp, fac = "Site", type = "fan", dist_method = "euclidean", hclust_method = "complete") ### hierarchical clustering (complete)

ggtree(cluster, layout="rectangular") + geom_nodepoint() + geom_tiplab(hjust = -0.1) + xlim(0,0.03)

### EXERCISE 11: IMPORTING OUTLINE DATA
### Inspect the outline file by opening the file in notepad
### Import the handaxe outline file through the import_tps() function
### Import the associated metadata through the read.csv() function
### Use the View() function to open the data table
### Use the print() function to view the raw coordinate data
### Convert the new handaxe function to an outline file (with the data) using the Out() function
### Explore the file using the functions in Momocs (?panel)

handaxe <- import_tps("handaxe.tps", curves = TRUE) ### .tps file with the outline data
handaxe.data <- read.csv("handaxe.csv", header = T, row.names = 1) ### handaxe dataset
View(handaxe.data)
print(handaxe)

handaxe.data$MIS = factor(handaxe.data$MIS, c("MIS 13", "MIS 11", "MIS 9", "MIS 7", "MIS 4/3")) ### reorder factor according to Marine Isotope Stage (MIS)
handaxe.data$Context = factor(handaxe.data$Context, c("Warren Hill", "Boxgrove", "Bowman's Lodge", "Elveden", "Swanscombe", "Broom", "Furze Platt", "Cuxton", "Pontnewydd", "Lynford")) ### reorder factor according to context
summary(handaxe.data$MIS) ### count data for the different MISs
summary(handaxe.data$Context) ### count data for the different archaeological contexts

outlinefile <- Out(handaxe$coo, fac = handaxe.data) ### creation of an outline file with the database supplying metadata
outlinefile ### call the outline file
panel(outlinefile) ### panel all examples
panel(outlinefile, fac = "MIS") ### panel coloured by MIS stages
panel(outlinefile, fac = "Context") ### panel coloured by context

### EXERCISE 12: PRIOR NORMALISATION
### Prior normalisation (EFA), all examples should be centred, scaled and closed (if open). These can be done using functions in Momocs.
### Use the coo_close(), coo_centre() and coo_scale() functions and normalise the outlines.
### Check with the stack() function.

outlinefile <- coo_close(outlinefile) ### ensure all outlines are closed
stack(outlinefile, title="") ### view outlines
outlinefile <- coo_center(outlinefile) ### centre outlines to a common centroid (0,0)
stack(outlinefile, title="") ### view outlines
outlinefile <- coo_scale(outlinefile) ### scale outlines to a common centroid size
stack(outlinefile, title="") ### stack for visual examination (see ?pile for further information)

### EXERCISE 13: HARMONICS AND EFA CREATION
### Use the following functions to determine the right levels of harmonics necessary for 99.9 harmonic power:
### calibrate_harmonicpower_efourier()
### calibrate_reconstructions_efourier()
### calibrate_deviations_efourier()
### Use the efourier() function to create the EFA class file and examine its contents

calibrate_harmonicpower_efourier(outlinefile) ### confirm how many harmonics equate to 99.9% harmonic power (may take some time!)
calibrate_reconstructions_efourier(outlinefile, range=1:20) ### confirm through reconstruction (of a random example)
calibrate_deviations_efourier(outlinefile) ### confirm through analysis of centroid deviations

efourierfile <- efourier(outlinefile, nb.h = 38, smooth.it = 0, norm = TRUE, start = FALSE) ### creation of EFA class (38 harmonics); normalisation is suitable in this instance given previous procedures.
efourierfile ### calls the file detailing the created OutCoe object (data and factors)

### EXERCISE 14: DATA EXPLORATION
### Using the EFA class object and the embedded data perform all of the following:
### Principal Component Analysis (PCA)
### Discriminant Analysis (LDA/CVA)
### MANOVA / Procrustes ANOVA
### Mean shape examination
### Regression (Lineal/Multiple)
### Cluster Analysis (Hierarchical Cluster Analysis etc.)

### Principal Component Analysis...

pca1 <- PCA(efourierfile, scale. = FALSE, center = TRUE, fac = Database) ### creation of PCA class
scores <- data.frame(pca1$x) ### creation of scores into a dataframe
scores <- scores[match(rownames(handaxe.data),rownames(scores)),] ### match to row ID on database
handaxe.data <- cbind(handaxe.data, scores) ### column bind to database
rm(scores) ### removes scores (no longer necessary)
View(handaxe.data) ### views the main database (again!)

scree(pca1) ### produces a tibble of the proportion and cumulative percentage for the PC inertia
scree_plot(pca1)
PCcontrib(pca1) ### visualises the main shape changes among all handaxe (PCs)
plot(pca1, xax = 1, yax = 2, points = FALSE, center.origin = FALSE, zoom = 1, grid = TRUE, pos.shp = "xy", size.shp = 0.4, ellipses = FALSE, chull = FALSE, chull.filled = FALSE, eigen = FALSE, rug = FALSE, title = "Principal Component Analysis (PC1 vs. PC2): XY Warps", labelsgroups=FALSE) ### produces a principal component plot (first two axes)
plot(pca1, pca1$MIS, xax = 1, yax = 2, points = FALSE, center.origin = FALSE, zoom = 1, grid = TRUE, morphospace = FALSE, ellipses = TRUE, conf.ellipses = 0.99, chull = FALSE, chull.filled = FALSE, eigen = FALSE, rug = FALSE, title = "Principal Component Analysis (PC1 vs. PC2): Confidence Ellipses (99%)", labelsgroups=TRUE, cex.labelsgroups=1) ### produces a principal component plot (first two axes)

ggplot(database, aes(MIS, PC1)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.4) +  coord_flip() + labs(x = "Marine Isotope Stage (MIS)", y = "Principal Component 1 (67.29%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) ### PC1 scores categorised by MIS
ggplot(database, aes(MIS, PC2)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.4) +  coord_flip() + labs(x = "Marine Isotope Stage (MIS)", y = "Principal Component 2 (11.46%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) ### PC1 scores categorised by context
ggplot(database, aes(Context, PC1)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.5) +  coord_flip() + labs(x = "Context", y = "Principal Component 1 (67.29%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) ### PC2 scores categorised by MIS
ggplot(database, aes(Context, PC2)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.5) +  coord_flip() + labs(x = "Context", y = "Principal Component 2 (11.46%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) ### PC2 scores categorised by context

figurea <- ggplot(database, aes(MIS, PC1)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.4) +  coord_flip() + labs(x = "Marine Isotope Stage (MIS)", y = "Principal Component 1 (67.29%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) ### PC1 scores categorised by MIS
figureb <- ggplot(database, aes(MIS, PC2)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.4) +  coord_flip() + labs(x = "Marine Isotope Stage (MIS)", y = "Principal Component 2 (11.46%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) ### PC1 scores categorised by context
figurec <- ggplot(database, aes(Context, PC1)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.5) +  coord_flip() + labs(x = "Context", y = "Principal Component 1 (67.29%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) ### PC2 scores categorised by MIS
figured <- ggplot(database, aes(Context, PC2)) + geom_boxplot(colour = "#E69F00", fill = "#ffd475", width = 0.5) +  coord_flip() + labs(x = "Context", y = "Principal Component 2 (11.46%)") + theme(text = element_text(size=9), axis.text.x = element_text(size=9), axis.text.y = element_text(size=9)) ### PC2 scores categorised by context
figure <- plot_grid(figurea, figurec, figureb, figured, labels= "AUTO", ncol = 2, align = 'v') #synthesis of the four figures
plot(figure) ### plots the figure
ggsave("Figure.tiff", plot = last_plot(), dpi = 400, units = "mm", height = 150, width = 250) ##saves the file

### Discriminant Analysis (LDA/CVA)...

lda1 <- LDA(pca1, fac = "MIS") ### creation of a discriminant analysis by MIS
lda1  ### details of the discriminant analysis
plot(lda1, xax = 1, yax = 2, points = TRUE, pch = 20, cex = 0.6, center.origin = FALSE, zoom = 1.8, grid = TRUE, pos.shp = "circle", size.shp = 0.6, ellipses = TRUE, ellipsesax = FALSE, conf.ellipses = 2/3, chull = FALSE, chull.filled = FALSE, eigen = FALSE, rug = FALSE) ### plots an LDA

### Statistical testing through a MANOVA...

MANOVA(pca1, "MIS", test = "Hotelling", retain = 0.99) ### MANOVA vs. MIS
mismanovapw <- MANOVA_PW(pca1, "MIS", retain = 0.99) ### pairwise values 
mismanovapw$stars.tab ### pairwise values represented by stars

MANOVA(pca1, "Context", test = "Hotelling", retain = 0.99) ### MANOVA against context
contextmanovapw <- MANOVA_PW(pca1, "Context", retain = 0.99) ### pairwise values
contextmanovapw$stars.tab ### pairwise values represented by stars

### Mean shapes...

handaxems <- mshapes(efourierfile, FUN = mean, "MIS") ### change FUN for median
coo_plot(handaxems$shp$`MIS 13`)
coo_draw(handaxems$shp$`MIS 11`, border='forestgreen')
coo_draw(handaxems$shp$`MIS 9`, border='red')
coo_draw(handaxems$shp$`MIS 7`, border='yellow')
coo_draw(handaxems$shp$`MIS 4/3`, border='blue')

tps_arr(handaxems$shp$`MIS 13`, handaxems$shp$`MIS 11`)
tps_iso(handaxems$shp$`MIS 13`, handaxems$shp$`MIS 11`)
tps_grid(handaxems$shp$`MIS 13`, handaxems$shp$`MIS 11`)

### Regression (Lineal/Multiple)...

ggplot(handaxe.data, aes(Length, PC1)) + geom_point(size = 1, pch = 16, alpha = 0.4, colour = "#E69F00", fill = "#ffd475") + geom_smooth(method = "lm", se = FALSE, colour = "grey") + theme(text = element_text(size=8), axis.text = element_text(size = 8)) + xlab("Length (mm)") + ylab("Principal Component 1 (67.29%)")
handaxelm <- lm(Length~PC1, data = handaxe.data)
cor.test(handaxe.data$Length, handaxe.data$PC1)

### Cluster analysis...

handaxe.cluster <- CLUST(pca1, fac = "MIS", type = "fan", dist_method = "euclidean", hclust_method = "complete") ### hierarchical clustering (complete)
ggtree(handaxe.cluster, layout = "fan", branch.length = "none") + theme_tree() + geom_nodepoint() + geom_point(size = 1) 

### Self-organising maps (SOM)...

install.packages('kohonen')
library(kohonen)
grid.som <- somgrid(xdim = 10, ydim = 10, topo = "hexagonal")
handaxex <- scale(pca1$x[,1:5])

set.seed(222)
handaxemap <- som(handaxex, grid = grid.som, alpha = c(0.05, 0.01), radius = 1)
handaxemap

plot(handaxemap, type = "changes")
plot(handaxemap,  main = "Self Organising Map (SOM) of PC variance")
plot(handaxemap, type = "count") # counts in each node
plot(handaxemap, main = "Self Organising Map (SOM) of PC variance", type = "mapping") # counts in each node

### Extracting other measures...

outlinefile %>% coo_elongation() #elongation

               