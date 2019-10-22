cd /gpfs/home/guosa/PRJNA577192/srr
fastq-dump --skip-technical --split-files -X 500000 --gzip SRR10290287

wget https://raw.githubusercontent.com/Shicheng-Guo/SmartMeth/master/bin/smartbismark.pl -O smartbismark.pl
perl smartbismark.pl --input config.txt --genome hg19 --server MCRI --queue shortq --submit submit


