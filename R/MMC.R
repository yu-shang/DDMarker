MMC = function(gene = c(), IDType = "Entrez", PV = 0.05, out = FALSE, mirna = FALSE, pvalue = FALSE)
{
  data("MMC");
  if(!require(DDMarkerData))
  {
    require(DDMarkerData);
    data(MMC0);
    if(MMC0$V > MMC$V)
      MMC = MMC0;
  }
  if(mode(gene) == "NULL")  gene = MMC$gene;
  if(is.null(dim(gene)))
  {
    if(length(gene) < 1)
    {
      return (FALSE);
    }
  }
  if(length(gene) > 0)
  {
    if(!is.null(dim(gene)))
    {
      gene2 = cbind(as.numeric(gene[,2]), as.numeric(gene[,3]));
      gene = gene[,1];
      names(gene) = gene;
    }
  }
  genes = names(gene);
  library("GOstats");
  library("pathview");
  Marker2KEGG = GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params", geneSetCollection = MMC$DDMarkerKEGG, geneIds = genes, universeGeneIds = MMC$DDMarkerUN, pvalueCutoff = PV, testDirection = "over");
  Marker2KEGG = hyperGTest(Marker2KEGG);
  #Marker2KEGG = summary(Marker2KEGG);
  if(length(Marker2KEGG@oddsRatios[which(Marker2KEGG@pvalues<PV)]) > 0)
  {
    Marker2KEGG = data.frame("KEGGID" = names(Marker2KEGG@oddsRatios[which(Marker2KEGG@pvalues<PV)]), "nROW" = as.numeric(0));
  }
  #print(Marker2KEGG)
  cd = rainbow(30);
  path.out = list();
  gene[] = 0.8;
  gene = cbind(gene, gene2);
  names(gene) = genes;
  if(out){
    out = nrow(Marker2KEGG);
  }else{
    out = min(out, nrow(Marker2KEGG));
    if(out == 0)  out = max(out, nrow(Marker2KEGG));
  }
  i = 0;
  if(out > 0)
  {
    for(i in 1:out)
    {
      path.out[[i]] = pathview(gene.data = gene, pathway.id = Marker2KEGG$KEGGID[i], limit = list(gene = 1), bins = list(gene = 12), low = list(gene = cd[11:15], cpd = "blue"), mid = list(gene = c(cd[6:10]), cpd = "gray"), high = list(gene = cd[1:5], cpd = "yellow"), species = "hsa", out.suffix = Marker2KEGG$KEGGID[i], kegg.native = TRUE, plot.col.key = T);
    }
  }
  if(mirna)
  {
    data(data0);
    temp = names(gene)[grep("has-miR", names(gene))];
    for(i in 1:length(data0$MP))
    {
      tt = names(gene)[grep(names(data0$MP)[i], names(gene))];
      if(length(tt) > 0)
      {
        tt = c(1, 0.8, 0.8);
        tt = t(as.matrix(tt));
        rownames(tt) = data0$MP[i];
        gene = rbind(gene, tt);
      }
    }
    path.out[[i+1]] = pathview(gene.data = gene, pathway.id = as.character("05206"), limit = list(gene = 1), bins = list(gene = 12), low = list(gene = cd[11:15], cpd = "blue"), mid = list(gene = c(cd[6:10]), cpd = "gray"), high = list(gene = cd[1:5], cpd = "yellow"), species = "hsa", out.suffix = "MicroRNAs in cancer", kegg.native = TRUE, plot.col.key = T);
  }
  return (path.out);
}

