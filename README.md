RegeneronSTS Submission Code
===
###### Tags: `Single Cell` `small intestinal`  `IBD` 

*Editor date : 20210906*
*Update date : 20210906*

If used in a project, please cite "PerpetualOwl" - dm PerpetualOwl#5670 on Discord for correspondence.

Code may be deprecated, has not been checked.

---

# Single Cell RNA-seq


Dataset paper: https://pubmed.ncbi.nlm.nih.gov/29144463/

## Datasets description

Control: 4
Salmonella: 2
H. Poly Day 3: 2
H. Poly Day 10: 2


## Data Download (terminal)
To begin, you download the data from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92332
``` linux
wget 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92332/suppl/GSE92332_RAW.tar'
```
Unpack the folder, and make sure it is named correctly.
```
tar -xvf GSE92332_RAW.tar
```


## Setting Up R
Load the following libraries.
``` {r}
library(Seurat)
library(ggplot2)
library(grid)
library(sctransform)
library(dplyr)
library(stringr)
library(patchwork)
```


## Convert to Seurat
Following code changes file names and file organization to make the files 10X format.
```{r}
GEO_ID = "GSE92332"
GEO_ID_RAW = paste0(GEO_ID,"_RAW")
setwd(paste0("scRNAseq/",GEO_ID))
fs=list.files(paste0("./",GEO_ID_RAW,"/"), "^GSM")
fs

samples=str_split(fs,'_',simplify = T)[,1]

lapply(unique(samples),function(x){
	y=fs[grepl(x,fs)]
	folder=paste0(GEO_ID_RAW,"/", str_split(y[1],'_',simplify = T)[,1])
	dir.create(folder,recursive = T)
	file.rename(paste0(GEO_ID_RAW,"/",y[1]),file.path(folder,"barcodes.tsv.gz"))
	file.rename(paste0(GEO_ID_RAW,"/",y[2]),file.path(folder,"features.tsv.gz"))
	file.rename(paste0(GEO_ID_RAW,"/",y[3]),file.path(folder,"matrix.mtx.gz"))
})

samples=list.files(GEO_ID_RAW) 
samples

sceList = lapply(samples,function(pro){
  folder=file.path(GEO_ID_RAW,pro)
  CreateSeuratObject(counts = Read10X(folder),
                     project = pro )
})

alldata <- merge(sceList[[1]], c(sceList[[2]], sceList[[3]], sceList[[4]], sceList[[5]], sceList[[6]], sceList[[7]], sceList[[8]], sceList[[9]], sceList[[10]], sceList[[11]], sceList[[12]], sceList[[13]], sceList[[14]], sceList[[15]], sceList[[16]], sceList[[17]], sceList[[18]], sceList[[19]], sceList[[20]], sceList[[21]], sceList[[22]], sceList[[23]], sceList[[24]], sceList[[25]], sceList[[26]], sceList[[27]], sceList[[28]], sceList[[29]], sceList[[30]], sceList[[31]], sceList[[32]], sceList[[33]], sceList[[34]], sceList[[35]], sceList[[36]], sceList[[37]], sceList[[38]], sceList[[39]], sceList[[40]], sceList[[41]], sceList[[42]], sceList[[43]], sceList[[44]], sceList[[45]], sceList[[46]]), add.cell.ids = samples)

```
Store each independent Seurat object in a list
``` {r}
object <- vector(mode = "list", length = 60)
for (item in 51:60){
directory <- paste("~/GSE92332_RAW/GSM28394", as.character(item), sep="")
data <- Read10X(data.dir=directory)
object[[item]] <- CreateSeuratObject(counts=data)
}
```
Merge the objects we created earlier
```{r}
CtrlObject <- merge(Object[[51]], y = c(Object[[52]], Object[[53]], Object[[54]]))
SalmObject <- merge(Object[[59]], y = Object[[60]])
Hpoly3Object <- merge(Object[[57]], y = Object[[58]])
Hpoly10Object <- merge(Object[[55]], y = Object[[56]])
```


## Data
For each of the four objects, perform the following. Only the process for CtrlObject will be shown.
```{r}
object <- CtrlObject #to make it easier to change objects
```
Filter the mitochondrial-rna-rich-cells because they are often outliers
```{r}
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "mt-")
VlnPlot(object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```
Based on the nFeature_RNA and percent_mt, decide the parameters for the subset.

```{r}
object <- subset(object, subset = nFeature_RNA >  500 & nFeature_RNA <  2300 & percent.mt < 5)
#I chose between 500 and 2300 features, and less than 5% mitochondrial RNA
```
Normalize, filter out variable, scale, PCA:
```{r}
object <- NormalizeData(object)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(object)
object <- ScaleData(object, features = all.genes)
object <- RunPCA(object, features = VariableFeatures(object = object))
```
Run Jackstraw to see the different dimensionality (optional)
```{r}
object <- JackStraw(object, num.replicate = 100)
object <- ScoreJackStraw(object, dims = 1:20)
JackStrawPlot(object, dims = 1:15)
```
Run ElbowPlot
```{r}
ElbowPlot(object)
```
Based on the result, decide how many dimensions to continue with. Generally it is better to have too many rather than too few, look for when the drop levels out in the plot. I chose 15 for this example.
```{r}
d <- 15
object <- FindNeighbors(object, dims = 1:d)
object <- FindClusters(object, resolution = 0.5)
```
The below code is to create an individual UMAP plot for the condition
```{r}
object <- RunUMAP(object, dims = 1:d)
DimPlot(object, reduction = "umap")
```
Add idents to keep track of them, then create a merged object
``` {r}
CtrlObject$orig.ident <- "CtrlObject"
SalmObject$orig.ident <- "SalmObject"
Hpoly3Object$orig.ident <- "Hpoly3Object"
Hpoly10Object$orig.ident <- "Hpoly10Object"
mergeobject <- merge(x = CtrlObject, y = c(Hpoly3Object, Hpoly10Object, SalmObject), add.cell.ids = c("c", "h3", "h10", "s"))
```

## Cluster Annotation
Must be manually done. The following can be used to find the gene markers in a specific cluster, other more specific searches can be performed using parameters found on Seurat's Reference webpage.
``` {r}
cluster0.markers <- FindMarkers(object, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5) #change n to display more top markers
```

The following renames clusters, the ones I used are included
```{r}
new.cluster.ids <- c("T Cells", "Neuron", "Enterocyte Precursor", "Enterocyte", "Endothelial", "Goblet", "Paneth", "Tuft", "Macrophage", "B Cell", "NK Cell")
names(new.cluster.ids) <- levels(object)
object <- RenameIdents(object, new.cluster.ids)
```


## UMAP
Create a UMAP representation of all of the data
```{r}
mergeobject <- RunPCA(object, features = VariableFeatures(object = mergeobject))
mergeobject <- RunUMAP(mergeobject, dims = 1:15)
DimPlot(mergeobject, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.5) + NoLegend()
```
![](https://i.imgur.com/JBFoFNP.png)



## Heatmap for Cell-type signatures
For any Seurat object, the following will create a DO Heatmap. Shown is the creation using "CtrlObject".
```{r}
object <- CtrlObject
object.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- object.markers %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC)
DoHeatmap(object, features = top10$gene) + NoLegend()
```
![](https://i.imgur.com/kQvUHVS.png)


## DEG

Create a .csv file with logFC and adjusted P Value for each gene between two conditions. Shown below is an example between the Control and Salmonella conditions.

``` {r}
SalmvsCtrl <- FindMarkers(mergeobject, test.use = "MAST", group.by = "orig.ident", ident.1 = "SalmObject", ident.2 = "CtrlObject")
write.csv(SalmvsCtrl, "~/GSE92332_RAW/SalmvsCtrl.csv")
```
Using the "group.by" parameter and the reference page https://satijalab.org/seurat/reference/findmarkers, it is possible to test nearly any condition.
Manipulate this file using ggplot2 to create plots and other tools like Metascape for further analysis.
