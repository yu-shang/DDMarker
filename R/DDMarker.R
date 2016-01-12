.onLoad <- function(lib, pkg)
{
  library.dynam("DDMarker", pkg, lib);
}

ddmarker = function(data = c("ZZZ3","O75019","PCDHA3", "uc003pjd.1", "345456", "hsa-miR-126"), pre = "Mix", file = FALSE, type = "All", Seq = "F", Path = FALSE, PV = 0.05, out = FALSE, pvalue = FALSE)
{
  ktype = 0;
  kgorp = 0;
  kf = 0;
  if( (type == "both") || (type == "Both") || (type == "BOTH") )
    ktype = 1;
  if( (type == "B") || (type == "blood") || (type == "Blood") || (type == "BLOOD") || (type == "b") )
    ktype = 2;
  if( (type == "U") || (type == "u") || (type == "Urine") || (type == "urine") || (type == "URINE") )
    ktype = 3;
  if( (type == "A") || (type == "a") || (type == "All") || (type == "all") || (type == "ALL") )
    ktype = 4;
  if( (pre == "mix") || (pre == "Mix") || (pre == "MIX") )
    kgorp = 1;
  if( (pre == "G") || (pre == "g") || (pre == "gene") || (pre == "Gene") || (pre == "GENE") )
    kgorp = 2;
  if( (pre == "P") || (pre == "p") || (pre == "protein") || (pre == "Protein") || (pre == "PROTEIN") )
    kgorp = 3;
  if( (pre == "M") || (pre == "m") || (pre == "mi") || (pre == "Mi") || (pre == "MiRNA") || (pre == "mirna") )
    kgorp = 4;
  if( (pre == "E") || (pre == "e") || (pre == "Entrez") || (pre == "entrez") || (pre == "ENTREZ") )
    kgorp = 5;
  if( (pre == "I") || (pre == "i") || (pre == "Isoform") || (pre == "isoform") || (pre == "ISOFORM") )
    kgorp = 6;
  if(file != "FALSE")
  {
    kf = 1;
    cat("The results of prediction will write to the file : ");
    cat(file);
    cat(" .\n\n");
  }
  if(ktype == 0)
  {
    cat("The data type can not be recognised.\n ");
    ktype = 1;
  }
  if(kgorp == 0)
  {
    cat("The prediction type can not be recognised.\n");
    kgorp = 1;
  }
  cat("DDMarker will predict");
  if(kgorp == 2)
    cat(" the GENE ");
  if(kgorp == 3)
    cat(" the PROTEIN ");
  if(kgorp == 4)
    cat(" the MiRNA ");
  if(kgorp == 4)
    cat(" the ENTREZ ");
  if(kgorp == 4)
    cat(" the ISOFORM ");
  if(kgorp == 1)
    cat(" all the symbols in the deductive system, ");
  cat("whether can be detected");
  if(ktype == 1)
    cat(" both in blood and urine .\n");
  if(ktype == 2)
    cat(" in blood .\n");
  if(ktype == 3)
    cat(" in urine .\n");
  if(ktype == 4)
    cat(" in all extracellular circulating .\n");
  data = as.character(data);
  if( (Seq == "t") || (Seq == "T") || (Seq == "True") || (Seq == "TRUE") || (Seq == "true") )
  {
    if(!require(DDMarkerData))
    {
      cat("DDMarker needs DDMarkerData package to predict the sequence.\nNow trying to install DDMarkerData ... \n");
      if(!require(devtools))
      {
        cat("DDMarkerData package require devtools package.\nNow installing ... \n");
        install.packages("devtools");
        if(!require(devtools))
        {
          cat("devtools can not be installed.\nDDMarker can not predict the sequence.\nPlease use the symbols of the genes, the proteins, the micro RNAs, and the isoforms instead.\nBlast can be used for the annotation of sequence to name.\n");
          return (FALSE);
        }
      }
      require(devtools);
      install_github('yu-shang/DDMarkerData');
      if(!require(DDMarkerData))
      {
        cat("DDMarkerData can not be installed.\nDDMarker can not predict the sequence.\nPlease use the symbols of the genes, the proteins, the micro RNAs, and the isoforms instead.\nBlast can be used for the annotation of sequence to name.\n");
        return (FALSE);
      }
    }
    require(DDMarkerData);
    data(EC0);
    temp = c(-1, -1, -1);
    datami = c();
    results = c();
    for(i in 1:length(data))
    {
      tt = data[1];
      tc = grep(tt, EC0$MS);
      if(length(tc) > 0)
      {
        datami = c(datami, tc);
      }
      tc = grep(tt, EC0$BP);
      if(length(tc) > 0)
      {
        for(j in 1:length(tc))
        {
          temp = rbind(temp, EC0$BP[tc[i],1], 1, 0);
        }
      }
      tc = grep(tt, EC0$BG);
      if(length(tc) > 0)
      {
        for(j in 1:length(tc))
        {
          temp = rbind(temp, EC0$BG[tc[i],1], 1, 0);
        }
      }
      tc = grep(tt, EC0$UP);
      if(length(tc) > 0)
      {
        for(j in 1:length(tc))
        {
          temp = rbind(temp, EC0$UP[tc[i],1], 0, 1);
        }
      }
      tc = grep(tt, EC0$UG);
      if(length(tc) > 0)
      {
        for(j in 1:length(tc))
        {
          temp = rbind(temp, EC0$UG[tc[i],1], 0, 1);
        }
      }
      tc = grep(tt, EC0$HSAD[,8]);
      results = c(results, tc);
      tc = grep(tt, EC0$HSAD[,9]);
      results = c(results, tc);
    }
    datami = unique(datami);
    datami = datami[order(datami)];
    results = unique(results);
    results = results(order(results));
    data = results;


    data = as.character(temp);
  }else{
    n = length(data);
    data(EC);
    results = list();
    if(!require(DDMarkerData))
    {
      require(DDMarkerData);
      data(EC1);
      if(EC1$V > EC$V)
        EC = EC1;
    }
    EC$HSAD = as.matrix(EC$HSAD);
    datami = c();
    if(kgorp == 2)
    {
      rownames(EC$HSAD) = EC$HSAD[,2];
      data = as.numeric(EC$HSAD[intersect(rownames(EC$HSAD), data),1]);
    }
    if(kgorp == 3)
    {
      rownames(EC$HSAD) = EC$HSAD[,3];
      data = as.numeric(EC$HSAD[intersect(rownames(EC$HSAD), data),1]);
    }
    if(kgorp == 5)
    {
      rownames(EC$HSAD) = EC$HSAD[,4];
      data = as.numeric(EC$HSAD[intersect(rownames(EC$HSAD), data),1]);
    }
    if(kgorp == 6)
    {
      rownames(EC$HSAD) = EC$HSAD[,5];
      data = substr(data, 1, 8);
      data = as.numeric(EC$HSAD[intersect(rownames(EC$HSAD), data),1]);
    }
    if(kgorp == 1)
    {
      rownames(EC$HSAD) = EC$HSAD[,2];
      temp = as.numeric(EC$HSAD[intersect(rownames(EC$HSAD), data),1]);
      rownames(EC$HSAD) = EC$HSAD[,3];
      temp = c(temp, as.numeric(EC$HSAD[intersect(rownames(EC$HSAD), data),1]));
      rownames(EC$HSAD) = EC$HSAD[,4];
      temp = c(temp, as.numeric(EC$HSAD[intersect(rownames(EC$HSAD), data),1]));
      rownames(EC$HSAD) = EC$HSAD[,5];
      datami = as.numeric(EC$Mi[intersect(rownames(EC$Mi), data),1]);
      data = substr(data, 1, 8);
      data = c(temp, as.numeric(EC$HSAD[intersect(rownames(EC$HSAD), data),1]));
      data = unique(data);
      data = data[order(data)];
    }
    if(kgorp == 4)
    {
      datami = as.numeric(EC$Mi[intersect(rownames(EC$Mi), data),1]);
    }
    temp = FALSE;
  }
  n = length(data);
  data(EC);
  if(!require(DDMarkerData))
  {
    require(DDMarkerData);
    data(EC1);
    if(EC1$V > EC$V)
      EC = EC1;
  }
  EC$HSAD = as.matrix(EC$HSAD);
  rownames(EC$HSAD) = EC$HSAD[,1];
  tall = EC$HSAD[data,];
  tall2 = EC$Mi[datami,];
  if(temp[1])
  {
    tall = as.matrix(tall);
    for(i in 1:nrow(temp))
    {
      tc = c("", temp[1], temp[1], "", "", temp[2], temp[3], "", "");
      tall = rbind(tall, tc);
    }
    rownames(tall) = c();
  }
  results = list();
  if(ktype == 2)
  {
    results$BLOOD = tall[which(tall[,6] == "1"), c(2,3,4,5)];
  }
  if(ktype == 3)
  {
    results$URINE = tall[which(tall[,7] == "1"), c(2,3,4,5)];
  }
  if(ktype == 1)
  {
    results$BLOOD = tall[which(tall[,6] == "1"), c(2,3,4,5)];
    results$URINE = tall[which(tall[,7] == "1"), c(2,3,4,5)];
  }
  if(ktype == 4)
  {
    results$MIRNA = tall2[, c(2,3,5,6)];
  }
  if(!is.null(results$BLOOD))
  {
    if(nrow(results$BLOOD) == 0)
    {cat("Under the Constraint of the ");cat(pre);cat(". Nothing can be detected in BLOOD");cat(".\n");}
    cat("\n");
    if(nrow(results$BLOOD) > 0)
    {
      cat("Under the Constraint of the ");cat(pre);cat(".\n");cat(nrow(results$BLOOD));cat(" can be detected in BLOOD");cat(".\n");
      #print(results$BLOOD);
    }
  }
  if(!is.null(results$URINE))
  {
    if(nrow(results$URINE) == 0)
    {cat("Under the Constraint of the ");cat(pre);cat(". Nothing can be detected in URINE");cat(".\n");}
    cat("\n");
    if(nrow(results$URINE) > 0)
    {
      cat("Under the Constraint of the ");cat(pre);cat(".\n");cat(nrow(results$URINE));cat(" can be detected in URINE");cat(".\n");
      #print(results$URINE);
    }
  }
  if(!is.null(results$MIRNA))
  {
    if(nrow(results$MIRNA) == 0)
    {cat("Under the Constraint of the ");cat(pre);cat(". No Micro RNA can be detected in Extracellular Circulating.\n");}
    cat("\n");
    if(nrow(results$MIRNA) > 0)
    {
      cat("Under the Constraint of the ");cat(pre);cat(".\n");cat(nrow(results$MIRNA));cat(" Micro RNA can be detected in Extracellular Circulating.\n");
      #print(results$MIRNA);
    }
  }
  cat("\n\n");
  if(kf == 1)
  {
    .C( "r_ddmarker", file = as.character(c(file, "*", "*")), iB = as.integer(nrow(results$BLOOD)), iU = as.integer(nrow(results$URINE)), iM = as.integer(nrow(results$MIRNA)), ktype = as.integer(ktype), kgorp = as.integer(kgorp));
    if(ktype == 1)
    {
      write.table(results$BLOOD, paste(file, "_BLOOD.txt", sep = "\t", quote = F));
      write.table(results$URINE, paste(file, "_URINE.txt", sep = "\t", quote = F));
    }
    if(ktype == 2)
    {
      write.table(results$BLOOD, paste(file, "_BLOOD.txt", sep = "\t", quote = F));
    }
    if(ktype == 3)
    {
      write.table(results$URINE, paste(file, "_URINE.txt", sep = "\t", quote = F));
    }
    if(ktype == 4)
    {
      write.table(results$BLOOD, paste(file, "_BLOOD.txt", sep = "\t", quote = F));
      write.table(results$URINE, paste(file, "_URINE.txt", sep = "\t", quote = F));
      write.table(results$MIRNA, paste(file, "_MIRNA.txt", sep = "\t", quote = F));
    }
  }
  if(Path)
  {
    temp = c();
    tmr = "FALSE";
    if(pvalue)
    {
      pvalue = as.numeric(as.matrix(pvalue));
    }
    temp = tall2;
    if(nrow(temp) > 0)
    {
      temp = temp[,2];
      temp = cbind(temp, 0.8);
      temp = cbind(temp, 0.8);
      temp = rbind(temp, tall[,c(4,6,7)]);
    }else
    {
      temp = tall[,c(4,6,7)];
    }
    if( (kgorp == 1) || (kgorp == 4) )
    {
      tmr = "TRUE";
    }
    PATHWAY = MMC(gene = temp, IDType = "Entrez", PV = PV, out = out, mirna = tmr, pvalue = pvalue);
    results$PATHWAY = PATHWAY;
    return (results);
  }
  return (results);
}
