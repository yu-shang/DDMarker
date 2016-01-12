#include "stdio.h"
#include "time.h"

int init_ddmarker( char **rfile, int *sB, int *sU, int *sM, int *k1, int *k2 )
{
  FILE *P;
  int i = 0, j = 0;
  time_t ti;
  time(&ti);
  P = fopen(rfile[0], "w");
  fprintf(P, "%s\nThe results of the prediction is writing to the file : %s .\nDDMarker will predict the", ctime(&ti), rfile[0]);
  if(*k2 == 1)  fprintf(P, " all the symbols in the deductive system, ");
  if(*k2 == 2)  fprintf(P, " the GENE ");
  if(*k2 == 3)  fprintf(P, " the PROTEIN ");
  if(*k2 == 4)  fprintf(P, " the MiRNA ");
  if(*k2 == 5)  fprintf(P, " the ENTREZ ");
  if(*k2 == 6)  fprintf(P, " the ISOFORM ");
  fprintf(P, "whether can be detected");
  if(*k1 == 1)  fprintf(P, " both in blood and urine .\n");
  if(*k1 == 2)  fprintf(P, " in blood .\n");
  if(*k1 == 3)  fprintf(P, " in urine .\n");
  if(*k1 == 4)  fprintf(P, " in all extracellular circulating .\n");
  fprintf(P, "\n");
  if(*sB == 0)
  {
    fprintf(P, "Under the Constraint of the ");
    if(*k2 == 1)  fprintf(P, " Mix ");
    if(*k2 == 2)  fprintf(P, " GENE ");
    if(*k2 == 3)  fprintf(P, " PROTEIN ");
    if(*k2 == 4)  fprintf(P, " MiRNA ");
    if(*k2 == 5)  fprintf(P, " ENTREZ ");
    if(*k2 == 6)  fprintf(P, " ISOFORM ");
    fprintf(P, ". Nothing can be detected in BLOOD. \n");
  }
  else if(*sB > 0)
  {
    fprintf(P, "Under the Constraint of the ");
    if(*k2 == 1)  fprintf(P, " Mix ");
    if(*k2 == 2)  fprintf(P, " GENE ");
    if(*k2 == 3)  fprintf(P, " PROTEIN ");
    if(*k2 == 4)  fprintf(P, " MiRNA ");
    if(*k2 == 5)  fprintf(P, " ENTREZ ");
    if(*k2 == 6)  fprintf(P, " ISOFORM ");
    fprintf(P, ".\n%d can be detected in BLOOD.\nGENE\tENTREZ\tISOFORM\tPROTEIN\t\n", *sB);
  }
  if(*sU == 0)
  {
    fprintf(P, "Under the Constraint of the ");
    if(*k2 == 1)  fprintf(P, " Mix ");
    if(*k2 == 2)  fprintf(P, " GENE ");
    if(*k2 == 3)  fprintf(P, " PROTEIN ");
    if(*k2 == 4)  fprintf(P, " MiRNA ");
    if(*k2 == 5)  fprintf(P, " ENTREZ ");
    if(*k2 == 6)  fprintf(P, " ISOFORM ");
    fprintf(P, ". Nothing can be detected in URINE. \n");
  }
  else if(*sU > 0)
  {
    fprintf(P, "Under the Constraint of the ");
    if(*k2 == 1)  fprintf(P, " Mix ");
    if(*k2 == 2)  fprintf(P, " GENE ");
    if(*k2 == 3)  fprintf(P, " PROTEIN ");
    if(*k2 == 4)  fprintf(P, " MiRNA ");
    if(*k2 == 5)  fprintf(P, " ENTREZ ");
    if(*k2 == 6)  fprintf(P, " ISOFORM ");
    fprintf(P, ".\n%d can be detected in URINE.\nGENE\tENTREZ\tISOFORM\tPROTEIN\t\n", *sU);
  }
  if(*sM == 0)
  {
    fprintf(P, "Under the Constraint of the ");
    if(*k2 == 1)  fprintf(P, " Mix ");
    if(*k2 == 2)  fprintf(P, " GENE ");
    if(*k2 == 3)  fprintf(P, " PROTEIN ");
    if(*k2 == 4)  fprintf(P, " MiRNA ");
    if(*k2 == 5)  fprintf(P, " ENTREZ ");
    if(*k2 == 6)  fprintf(P, " ISOFORM ");
    fprintf(P, ". No Micro RNA can be detected in Extracellular Circulating. \n");
  }
  else if(*sM > 0)
  {
    fprintf(P, "Under the Constraint of the ");
    if(*k2 == 1)  fprintf(P, " Mix ");
    if(*k2 == 2)  fprintf(P, " GENE ");
    if(*k2 == 3)  fprintf(P, " PROTEIN ");
    if(*k2 == 4)  fprintf(P, " MiRNA ");
    if(*k2 == 5)  fprintf(P, " ENTREZ ");
    if(*k2 == 6)  fprintf(P, " ISOFORM ");
    fprintf(P, ".\n%d Micro RNA can be detected in Extracellular Circulating.\nCLASS\tECTYPE\tJournal\t\n", *sM);
  }
  fclose(P);
  return 1;
}
