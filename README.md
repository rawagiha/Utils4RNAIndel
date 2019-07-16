## RNAIndel maintenance

#### refGene2refCodingExon
Make a coding exon database in BED format
```
refGene2refCodingExon.py -r refGene.txt
```

#### extractExonicIndelFromVcf
Make an exonic indel database in VCF format (mostly to reduce the database size)
```
extractExonicIndelFromVcf.py -v INPUT_VCF -d DATA_DIR -o OUTPUT_VCF
```

#### annotateDbSnpIndelsWithNonCancerAf
Annotate dbSNP database with gnomAD non-cancer allele frequency
```
annotateDbSnpIndelsWithNonCancerAf.py -d DBSNP -g GNOMAD -f FASTA -o OUTPUT_VCF
```

#### NCBI_CDD.tar
Written by Dr. Jian Wang (St Jude Research Hosp). Preparing Conserved Domain Database in JSON.
```
tar xvzf NCBI_CDD.tar
cd NCBI_CDD
./NCBI_CDD.sh
```
 
