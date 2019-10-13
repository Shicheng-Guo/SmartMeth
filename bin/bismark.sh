cd ~/hpc/methylation/clep/wgbs/test/data
zcat ../../SRR9888304_1.fastq.gz | head -n 500000 | gzip -c > SRR9888304_1.fastq.gz
zcat ../../SRR9888304_2.fastq.gz | head -n 500000 | gzip -c > SRR9888304_2.fastq.gz
wget https://raw.githubusercontent.com/Shicheng-Guo/SmartMeth/master/bin/smartbismark.pl -O smartbismark.pl
perl smartbismark.pl --input SraRunTable.txt --genome hg19 --server MCRI --queue shortq --submit submit

