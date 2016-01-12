\name{DDMarker}
\alias{ddmarker}
\alias{DDMarker-method}
\alias{DDMarker}
\title{
Diagnose and Detect Markers
}
\description{
Diagnose and Detect Markers in Extracellular Circulating is a homo sapiens deductive system solving the markers in extracellular circulating. It entails the symbols of markers, like the genes, the proteins, the micro RNAs, and the isoforms, whether can be diagnose and detect in extracellular circulating, especially the blood serum and the urine for the biological and medicine significance. With the help of a homo sapiens annotation database in DDMarkerData package, DDMarker can even diagnose and detect the sequence among the genes, the proteins, the micro RNAs, and the isoforms. There are two main function in this package, \code{the ddmarker}, and \code{the DDMarkerMMC}, short for Minimal Metabolize Circulation. \code{DDMarkerMMC} entails the markers among the minimal metabolize circulation. \cr \cr
results = ddmarker(data, \ldots); \cr \cr
}
\arguments{
\item{data}{ Data vector. }
\item{pre}{ -pre should be ( "gene", "protein", "MiRNA", "entrez", "isoform", "mix" ), is a data type variable, telling DDMarker under which constraint DDMarker will do the diagnose and detect. If -pre is "mix", DDMarker will do the diagnose and detect on all constraints. \cr default: "mix" }
\item{file}{ -file should be a filename to output the results if want to save, or use "FALSE" to ignore. \cr default: "FALSE" }
\item{type}{ -type should be ( "both", "blood", "urine", "all" ), is a predict type variable, telling DDMarker which extracellular circulating DDMarker will do the diagnose and detect. If -pre is "all", DDMarker will do the diagnose and detect on all extracellular circulatings. \cr default: "all" }
\item{Seq}{ -Seq should be ( "T", "F" ). If -Seq is True, the DDMarker will predic the sequence first, then doing DDMarker under the results of that. The sequence can come from the genes, the proteins, the micro RNAs, and the isoforms, and the data vector must be the sequence vector. \cr default: "F" }
\item{Path}{ -Path should be ( "T", "F" ). If -Path is True, The DDMarker will do Minimal Metabolize Circulation, write the images into the files under the results of the diagnose and detect in ddmarker. \cr default: "F" }
\item{PV}{ -PV should be a double number from 0.00 to 1.00. It is a variable of -Path, telling which metabolize circulation should be written out. If -PV is 0.05, the DDMarker will only write the statistical significant results. Else if -PV is 1.00, the DDMarker will write all the metabolize circulations. \cr default: 0.05 }
\item{out}{ -out should be a integer number greater than 0, a variable of -Path, DDMarker writting the minimal quantity between -out and the minimal metabolize circulation, or use "FALSE" to ignore.  \cr default: "F" }
\item{pvalue}{ -pvalue should be a vector, is a vector of path. If -pvalue is not "False", DDMarker not only writes the minimal metabolize circulation, but also entails the degree of the markers. \cr default: "F" }
}
\value{
The R function, \code{ddmarker} returns an object of \code{list}: \cr 
\tabular{ll}{
BLOOD \tab An object of \code{list}, having \code{Gene}, \code{Entrez}, \code{Isoform} and \code{Protein}. \cr 
URINE \tab An object of \code{list}, having \code{Gene}, \code{Entrez}, \code{Isoform} and \code{Protein}. \cr 
MIRNA \tab An object of \code{list}, having \code{RNA_name}, \code{RNA_class}, \code{RNA_type} and \code{Journal}. \cr 
Gene \tab An entailment, entails the symbol of the genes which can be diagnosed and detected in \code{Extracellular Circulating} or \code{BLOOD} and \code{URINE}. \cr 
Entrez \tab A entailment, entails the symbol of the entrez which can be diagnosed and detected in \code{Extracellular Circulating} or \code{BLOOD} and \code{URINE}. \cr 
Isoform \tab A entailment, entails the symbol of the isoform which can be diagnosed and detected in \code{Extracellular Circulating} or \code{BLOOD} and \code{URINE}. \cr 
Protein \tab A entailment, entails the the symbol of the protein which can be diagnosed and detected in \code{Extracellular Circulating} or \code{BLOOD} and \code{URINE}. \cr 
RNA_name \tab A entailment, entails the symbol of the micro RNAs which can be diagnosed and detected in \code{Extracellular Circulating}. \cr 
RNA_class \tab A entailment, entails the class of the micro RNAs. \cr 
RNA_type \tab A entailment, entails the type of the micro RNAs, simply which EC part it can be diagnosed and detected. \cr 
Journal \tab A entailment, if existed, entails in which scientific research, it being studied. \cr 
}
The R function, \code{ddmarker} will write the results in a file named by -file, IFF -file is not FALSE by default. The file will be in the \code{Working Directory} . \cr
The R function, \code{ddmarker} will write the minimal metabolize circulation results in the files named by the MMC names, IFF -Path is not FALSE by default. The file will be in the \code{Working Directory} . \cr
}
\details{
data(data0);
results = ddmarker(data = data0$demo.ddm, pre = "Mix", file = "FALSE", type = "All", Seq = "F", Path = "FALSE", PV = 0.05, out = "FALSE", pvalue = "FALSE");
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
citation("DDMarker");
\cr 
\code{[1]} Yan Wang, et al. (2009) \emph{DMarker: A Bio-Marker Inference System for Human Diseases based
on Microarray Gene Expression Data} 2009 \cr
\code{[2]} Juan Cui, et al. (2008) \emph{Computational prediction of human proteins that can be secreted into the bloodstream} BIOINFORMATICS, Vol.24 no. 20 2008 pages 2370-2375 \cr
\code{[3]} Jiaxin Wang, et al. (2013) \emph{Computational Prediction of Human Salivary Proteins from Blood Circulation and Application to Diagnostic Biomarker Identification} PLoS ONE, DOI: 10.1371/journal.pone.0080211, 2013 \cr
\code{[4]} S Hong, et al. (2011) \emph{A Computational Method for Prediction of Excretory Proteins and Application to Identification of cancer markers in urine and application to gastric cancer} PLoS ONE,6(2):e16875, 2011 \cr
\code{[5]} \url{http://bioinfosrv1.bmb.uga.edu/DMarker/} \cr
}
\keyword{ DDMarker }
\keyword{ blood }
\keyword{ unine }
\keyword{ EC }
\keyword{ pathway }
\keyword{ gene }
\keyword{ protein }
\keyword{ isoform }
\keyword{ miRNA }
\keyword{ entrez }
\seealso{
\code{\link{DDMarker-package}} \code{\link{MMC-method}} \code{\link{EC}} \code{\link{demo}} \cr
}
\examples{

data(data0);
# load the demo of DDMarker

results = ddmarker(data0$demo.ddm);
# run dmarker with default parameters
# results is a variable of list with,
#   $BLOOD,
#   $BLOOD[,"Gene"], entails the symbol of the genes which can be diagnosed and detected in BLOOD,
#   $BLOOD[,"Entrez"],entails the symbol of the entrez which can be diagnosed and detected in BLOOD,
#   $BLOOD[,"Isoform"], entails the symbol of the isoform which can be diagnosed and detected in BLOOD,
#   $BLOOD[,"Protein"], entails the the symbol of the protein which can be diagnosed and detected in BLOOD,
#   $URINE,
#   $URINE[,"Gene"], entails the symbol of the genes which can be diagnosed and detected in URINE,
#   $URINE[,"Entrez"], entails the symbol of the entrez which can be diagnosed and detected in URINE,
#   $URINE[,"Isoform"], entails the symbol of the isoform which can be diagnosed and detected in URINE,
#   $URINE[,"Protein"], entails the the symbol of the protein which can be diagnosed and detected in URINE,
#   $MIRNA,
#   $MIRNA[,"RNA_name"], entails the symbol of the micro RNAs which can be diagnosed and detected in extracellular circulating,
#   $MIRNA[,"RNA_class"], entails the class of the micro RNAs,
#   $MIRNA[,"RNA_type"], entails the type of the micro RNAs, simply which EC part it can be diagnosed and detected,
#   $MIRNA[,"Journal"], if existed, entails in which scientific research, it being studied,
# or results = ddmarker(data = data0$demo.ddm, pre = "Mix", file = "FALSE", type = "All", Seq = "F", Path = "FALSE", PV = 0.05, out = "FALSE", pvalue = "FALSE");

results = ddmarker(data0$demo.ddm, file = "DDMarkerResults.txt");
# Saving the results both in R variable results and the file, named DMarkerResults.txt by -file.

results = ddmarker(data0$demo.seq, Seq = "TRUE");
# Diagnose and detect the sequences of markers in extracellular circulating.

results = ddmarker(data0$demo.path, path = "TRUE", pvalue = data0$demo.pm);
# Doing both DDMarker and minimal metabolize circulation with the degree of the markers.

results = ddmarker(data = data0$demo.pro, pre = "protein", type = "blood");
# Proteins are diagnosed and detected whether in blood.

results = ddmarker(data = data0$demo.gen, pre = "gene", type = "urine");
# Genes are diagnosed and detected whether in urine.

results = ddmarker(data = data0$demo.iso, pre = "isoform", type = "both");
# Isoforms are diagnosed and detected whether in blood or urine

results = ddmarker(data = data0$demo.ent, pre = "entrez", type = "all");
# Entrezs are diagnosed and detected whether in extracellular circulating.

results = ddmarker(data = data0$demo.mir, pre = "MiRNA");
# Micro RNAs are diagnosed and detected whether in extracellular circulating.

}
