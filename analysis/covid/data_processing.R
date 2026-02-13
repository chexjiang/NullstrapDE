base_dir='~/proj/NullstrapDE/covid/GSE158055'
setwd(base_dir)
{
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(magrittr)
  library(stringr)
  library(readr)
  fts <- read_tsv(paste0(base_dir,'/GSE158055_covid19_features.tsv.gz'))
  annot <- read_csv(paste0(base_dir,'/GSE158055_cell_annotation.csv.gz'))
  mdat <- readxl::read_xlsx(paste0(base_dir,'/GSE158055_sample_metadata.xlsx'))
  # which(is.na(mdat$`# High-throughput sequencing metadata template (version 2.1).`))
  meta.dat <- data.frame(mdat[21:304,2:ncol(mdat)])
  colnames(meta.dat) <- mdat[20,2:ncol(mdat)]
  rownames(meta.dat) <- mdat$`# High-throughput sequencing metadata template (version 2.1).`[21:304]
}

{
  meta.dat.2 <- meta.dat %>% dplyr::filter(`characteristics: Sample time` %in% c('control','progression'))%>%
    dplyr::filter(`characteristics: CoVID-19 severity` %in% c('mild/moderate', 'severe/critical'))
  meta.dat.2%<>%dplyr::filter(`characteristics: Sample type` %in% c('frozen PBMC', 'fresh PBMC'))%>%
    dplyr::filter(`characteristics: Single cell sequencing platform`=="10X 5'")
  annot2<-annot%>%dplyr::filter(sampleID %in% meta.dat.2$title)%>%
    dplyr::filter(majorType %in% c('Mono'))

  #Change the file format of GSE158055_covid19_counts.mtx.gz to be read by Read10X as described in https://github.com/satijalab/seurat/issues/4030.
  #Original data is too large to be read by Read10X.
  #Both folder part 1 and part 2 contain three files (barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz).
  #GSE158055_covid19_barcodes.tsv.gz is renamed to barcodes.tsv.gz of part1/part2.
  #GSE158055_covid19_features.tsv.gz is renamed to features.tsv.gz of part1/part2.
  #GSE158055_covid19_counts.mtx.gz is divided into matrix.mtx.gz of part1/part2.
  c1 <- Read10X('part1/', gene.column=1)
  c1.1 <- c1[,colnames(c1)%in%annot2$cellName]
  c2 <- Read10X('part2/', gene.column=1)
  c2.1 <- c2[,colnames(c2)%in%annot2$cellName]
  #c3 <- c1.1 + c2.1
  c3 <- cbind(c1.1, c2.1)

  # Create Seurat object
  dt <- CreateSeuratObject(
    counts = c3,
    project = "GSE158055",
    min.cells = 3,        # keep genes expressed in >= 3 cells
    min.features = 200    # keep cells with >= 200 detected genes
  )

  # Filter out cells with >10% MT reads
  dt[["percent.mt"]] <- PercentageFeatureSet(dt, pattern = "^MT-")
  dt <- subset(dt, subset = percent.mt <= 10)

  # Add sample- and cell-level metadata
  # Patient, Group, Batch based on annot2/meta.dat.2
  annot2$Patient <- meta.dat.2$Patients[match(annot2$sampleID, meta.dat.2$title)]

  cellinfo <- annot2 %>%
    transmute(
      cellName,
      sampleID,
      celltype,
      majorType,
      Group   = gsub("/.*", "", meta.dat.2$`characteristics: CoVID-19 severity`[match(sampleID, meta.dat.2$title)]),
      Batch   = meta.dat.2$`characteristics: Sex`[match(sampleID, meta.dat.2$title)],
      Patient = Patient
    ) %>%
    distinct(cellName, .keep_all = TRUE)

  rownames(cellinfo) <- cellinfo$cellName
  meta_aligned <- cellinfo[match(colnames(dt), rownames(cellinfo)), , drop = FALSE]
  #cellinfo$cellName <- NULL
  dt <- AddMetaData(dt, metadata = meta_aligned)


  # Output per-cell info
  write.table(cellinfo, "cellinfo.txt", sep="\t", quote=FALSE, row.names=FALSE)

  # Pseudobulk seurat
  keep <- !is.na(dt$Patient) & !is.na(dt$celltype)
  dt_pb <- subset(dt, cells = colnames(dt)[keep])
  pb_list <- AggregateExpression(
    object         = dt_pb,
    assays         = "RNA",
    group.by       = c("Patient"),
    slot           = "counts",                  # raw counts
    fun            = "sum",
    return.seurat  = FALSE                      # get matrices directly
  )
  pb_counts <- pb_list$RNA                       # genes x pseudobulk samples (matrix)
  write.table(pb_counts, "pb_counts.txt")

  dim(pb_counts)

  patients_vec <- colnames(pb_counts)
  meta2_reord <- meta.dat.2[match(patients_vec, meta.dat.2$Patients), , drop = FALSE]
  write.table(meta2_reord, "pb_meta.txt")


  # # Output count
  # counts <- GetAssayData(dt, slot = "counts")
  # write.table(counts,'counts.txt')


  # writeMM(obj = c3, file=paste0("COVID19.raw.txt"))
  # DN <- dimnames(c3)
  #
  # count.matrix <- readMM("COVID19.raw.txt")
  # dimnames(count.matrix) = DN
  # count.matrix%<>%as.matrix%<>%as.data.frame()
  #
  #
  # annot2$Patients <- meta.dat.2$Patients[match(annot2$sampleID,meta.dat.2$title)]
  #
  # #generate cellinfo
  # Cell <- annot2$cellName
  # Cell <- gsub(Cell,pattern='[-]',replacement='.')
  # Patients <- sapply(annot2$sampleID,FUN=function(x){return(meta.dat.2$`characteristics: Sex`[which(meta.dat.2$title==x)])})%>%as.vector()
  # Batch <- sapply(annot2$sampleID,FUN=function(x){return(meta.dat.2$`characteristics: Sex`[which(meta.dat.2$title==x)])})%>%as.vector()
  # Group <- sapply(annot2$sampleID,FUN=function(x){return(meta.dat.2$`characteristics: CoVID-19 severity`[which(meta.dat.2$title==x)])})%>%as.vector()%>%gsub(pattern='[/].*',replacement='')
  # cellinfo <- data.frame(Cell=Cell,Batch=Batch, Group=Group, Patient=Patients)
  #
  # #filter COVID19 data
  # dt <- CreateSeuratObject(counts=count.matrix, min.cells = 0.05*ncol(count.matrix))
  # dt[["percent.mt"]] <- PercentageFeatureSet(dt, pattern = "^MT-")
  # dt <- subset(dt, subset = percent.mt <= 10)
  # # g1 <- dt@assays$RNA@layers$counts%>%rownames()
  # g1 <- rownames(LayerData(dt, layer = "counts"))
  #
  # write.table(count.matrix[which(rownames(count.matrix) %in% g1),colnames(dt)],'counts.txt')
  # cellinfo%<>%dplyr::filter(Cell%in%gsub(pattern='[-]',replacement='.',colnames(dt)))
  # write.table(cellinfo[,c('Cell','Group','Batch')],'cellinfo.txt')
  #
  # #We aggregate cells for each patients in COVID19 data pseudobulk analysis
  # cellinfo$Batch=cellinfo$Patient
  # write.table(cellinfo[,c('Cell','Group','Batch')],'cellinfo_pseudobulk.txt')


}
