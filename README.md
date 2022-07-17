# RBPed
RBPed is a bioinformatics pipeline to identify differential RNA editing events between eCLIP-seq and RNA-seq

## Authors
Xiaolin Hu, Qin Zou, Li Yao and Xuerui Yang

## 1. python dependency
RBPed was compromised by samtools, REDItools V1, and R. Among those softwares, REDItools V1 depends on python2.7, besides, it depends on pysam 0.7.6
* you can install with code: python -m pip install pysam==0.7.6

Also, it depends on some bionformatics packages
* samtools

## 2. reference file
the ference genome (hg19) downloaded from GENCODE website. Known editing events downloaded from REDIportal database.
thos data can be downloaded within this pipeline.

## 3. Run
* ### A. downloading mapping files (bam) from ENCODE websites.
'''
wget http://srv00.recas.ba.infn.it/webshare/ATLAS/donwload/TABLE1_hg19.txt.gz
gunzip TABLE1_hg19.txt.gz
grep -v 'Region' TABLE1_hg19.txt |awk '{print $1,$2,$5}' | sed 's/ /\t/g' | sed 's/chr//g' | sort -k1,1 -k2,2n  > REDIportals.forREDItools.txt

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_40/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz
gunzip GRCh37.primary_assembly.genome.fa.gz
sed 's/chr//g' GRCh37.primary_assembly.genome.fa > hg19.ref.fa
'''
* ### B. sort, index bam files
'''
samtools  sort --threads 16  -o XXX.sorted.bam XXX.eCLIP.bam
samtools index XXX.sorted.ba
'''
* ### C. detection RNA editing events with REDItools.
'''
python2.7 REDItoolKnown.py  -i XXX.sorted.bam -f hg19.ref.fa -l REDIportals.forREDItools.txt -t 16 -c 10 -T 6-0  -p -e -d -u -m20  -v 0 -n 0.0 -o outputDir
'''
* ### D. merge RNA editing events in different replicate for different fractional RNA-seq dataset
* ### D. merge RNA editing events in different replicate for specific RBP's eCLIP
* ### E. merge RNA editing events for specific eCLIP and RNA-seq
* ### F. using R to calculate delta editing levels and fisher test P-values.


## License
RBPed is licensed under the MIT license



