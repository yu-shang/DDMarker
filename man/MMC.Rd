\name{MMC}
\alias{MMC}
\alias{mmc}
\alias{MMC-method}
\alias{Path}
\title{
Minimal Metabolize Circulation
}
\description{
Minimal Metabolize Circulation for Diagnose and Detect Markers in Extracellular Circulating \cr \cr
results = DDMarkerMMC(gene, \ldots); \cr \cr
}
\arguments{
\item{gene}{ Data vector. }
\item{IDType}{ -IDType shoulde be ( "Entrez", "GENE", "TCGA" ), "Entrez" denotes the entrez IDs, "GENE" denotes the gene IDs, "TCGA" denotes the TCGA IDs. \cr default: "Entrez" }
\item{PV}{ -PV should be a double number from 0.00 to 1.00. Telling which metabolize circulation should be written out. If -PV is 0.05, the DDMarkerMMC will only write the statistical significant results. Else if -PV is 1.00, the DDMarkerMMC will write all the metabolize circulations. \cr default: 0.05 }
\item{out}{ -out should be a integer number greater than 0, a variable of -Path, DDMarkerMMC writting the minimal quantity between -out and the minimal metabolize circulation, or use "FALSE" to ignore.  \cr default: "F" }
\item{mirna}{ -mirna shoulde be ( "T", "F" ). If -mirna is "TRUE", the DDMarkerMMC will not only do MMC function, but also write out the regulating relationship among Micro RNAs, no matter how less the statistical significant is. \cr default: "F" }
\item{pvalue}{ -pvalue should be a vector, only working significant when DDMarker is called by DDMarker. \cr default: "F" }
}
\value{
The R function, \code{DDMarkerMMC} returns an object of \code{list}, comes from the packages pathview: \cr 
\tabular{ll}{
plot.data.gene \tab data.frame returned by node.map function for rendering mapped gene nodes, including node name, type, positions (x, y), sizes (width, height), and mapped gene.data. This data is also used as input for pseduo-color coding through node.color function. Default plot.data.gene=NULL. \cr 
plot.data.cpd \tab same as plot.data.gene function, except for mapped compound node data. d plot.data.cpd=NULL. Default plot.data.cpd=NULL. Note that plot.data.gene and plot.data.cpd can't be NULL simultaneously. \cr 
}
The R function, \code{DDMarkerMMC} will write the result images in files named by the minimal metabolize circulations. The images will be in the \code{Working Directory} . \cr
}
\details{
results = DDMarkerMMC(gene = c(), IDType = "Entrez", PV = 0.05, out = FALSE, mirna = FALSE, pvalue = FALSE);
}
\author{
Yu Shang (JLU & UGA) \email{yushang@uga.edu} \cr 
Qiong Yu (JLU & UGA) \email{yuqiong@uga.edu} \cr 
Huansheng Cao (UGA) \email{hshcao@uga.edu} \cr 
Guoqing Liu (IMUST & UGA) \email{gqliu@uga.edu} \cr 
Wei Du (JLU & UGA) \email{weidu@uga.edu} \cr 
Yan Wang (JLU & UGA) \email{wy6868@hotmail.com} \cr 
Hao Wu (BIT & UGA) \email{wuhao@uga.edu} \cr 
Jingyu Yang (GSUT & UGA) \email{jyyyjy@uga.edu} \cr 
Xiufeng Liu (GZUCM & UGA) \email{xfliu@uga.edu} \cr 
Yuan Guo (SAHCQMU & UGA) \email{yuanguo@uga.edu} \cr 
Ying Xu (JLU & UGA) \email{xyn@bmb.uga.edu} \cr \cr 
Maintainer: Yu Shang (JLU & UGA) \email{yushang@uga.edu} \cr 
}
\references{
citation("DDMarkerMMC");
}
\keyword{ DDMarker }
\keyword{ MMC }
\keyword{ EC }
\keyword{ pathway }
\keyword{ metabolism }
\keyword{ gene }
\keyword{ protein }
\keyword{ isoform }
\keyword{ miRNA }
\keyword{ entrez }
\seealso{
\code{\link{DDMarker-package}} \code{\link{DDMarker-method}} \code{\link{EC}} \code{\link{MMC}} \code{\link{demo}} \cr
}
\examples{

data(data0);
# load the demo of DDMarker

results = DDMarkerMMC();
# run MMC with default parameters
# use a demo data named MMC.rda comes from the packages DDMarkerMMC

results = DDMarkerMMC(data0$demo.pm, mirna = "TRUE");
# run MMC uses a demo MiRNA data

results = DDMarkerMMC(data0$MI, mirna = "TRUE");
# run MMC uses another demo MiRNA data

# results is a variable of list comes from the packages pathview,
#   $plot.data.gene, data.frame returned by node.map function for rendering mapped gene nodes, including node name, type, positions (x, y), sizes (width, height), and mapped gene.data. This data is also used as input for pseduo-color coding through node.color function. Default plot.data.gene=NULL.
#   $plot.data.cpd, same as plot.data.gene function, except for mapped compound node data. d plot.data.cpd=NULL. Default plot.data.cpd=NULL. Note that plot.data.gene and plot.data.cpd can't be NULL simultaneously.
# The results returned by keggview.native and codekeggview.graph are both a list of graph plotting parameters. These are not intended to be used externally.

}
