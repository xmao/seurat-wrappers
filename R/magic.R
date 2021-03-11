#' @include internal.R
#' @importFrom Seurat SetAssayData
#' @importFrom Rmagic magic
#'
NULL

#' Run Magic
#'
#' @param object Seurat object
#' @param assay Assay used as input to MAGIC, default RNA
#' @param new.assay.name New created assay from MAGIC output
#' @param genes Genes imputed for, all genes by default
#' @param knn knn number for MAGIC algorithm
#' @param solver solver used in MAGIC
#'
#' @return Seurat object with new assay named MAGIC from MAGIC output
#'
#' @export
#'
#' @seealso \link[Rmagic package]{https://cran.r-project.org/web/packages/Rmagic/index.html}
#'
RunMAGIC <- function(
  object, assay = "RNA", new.assay.name = 'MAGIC',
  genes = "all_genes", knn = 15, solver = c("exact", "approximate"),
  verbose = TRUE
) {
  CheckPackage(package = "Rmagic", respository = "CRAN")
  solver <- match.arg(solver)

  assay <- assay %||% DefaultAssay(object = object)
  umi <- GetAssayData(object = object, assay = assay, slot = "counts")
  umi.norm <- GetAssayData(object = object, assay = assay, slot = "data")
  magic.out <- Rmagic::magic(t(as.matrix(umi.norm)), genes=genes, knn = knn, solver = solver)

  assay.obj <- CreateAssayObject(counts = umi)
  assay.obj <- SetAssayData(assay.obj, slot = "data", new.data = t(magic.out$result))
  object[[new.assay.name]] <- assay.obj
  DefaultAssay(object) <- new.assay.name

  object <- LogSeuratCommand(object = object)
  return(object)
}
