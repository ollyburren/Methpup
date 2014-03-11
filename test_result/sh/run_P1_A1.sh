echo "check qc ..."
/home/xin/sources/fastqc_v0.10.0/FastQC/fastqc -o test_result/fastqc_raw test/P1_A1.1.fq
/home/xin/sources/fastqc_v0.10.0/FastQC/fastqc -o test_result/fastqc_raw test/P1_A1.2.fq

echo "trim sequences..."
java -jar /home/xin/sources/Trimmomatic-0.27/trimmomatic-0.27.jar PE -phred33 -trimlog test_result/trim/P1_A1.log  test/P1_A1.1.fq  test/P1_A1.2.fq test_result/trim/P1_A1.1.paired.fq test_result/trim/P1_A1.1.unpaired.fq test_result/trim/P1_A1.2.paired.fq test_result/trim/P1_A1.2.unpaired.fq HEADCROP:18 TRAILING:20 SLIDINGWINDOW:4:20

/home/xin/sources/cutadapt-1.2.1/bin/cutadapt -a GGTCATAGCTGTTTCCTG -O 5 test_result/trim/P1_A1.1.paired.fq > test_result/trim/P1_A1.1.paired.nolinker.fq

/home/xin/sources/cutadapt-1.2.1/bin/cutadapt -a ACTGGCCGTCGTTTTACA -O 4 test_result/trim/P1_A1.2.paired.fq > test_result/trim/P1_A1.2.paired.nolinker.fq

echo "stitch forward and reverse reads..."
flash -M 222 -m 10 -o P1_A1 -d test_result/stitch test_result/trim/P1_A1.1.paired.nolinker.fq test_result/trim/P1_A1.2.paired.nolinker.fq

echo "check qc after stitching..."
/home/xin/sources/fastqc_v0.10.0/FastQC/fastqc -o test_result/fastqc_stitch test_result/stitch/P1_A1.extendedFrags.fastq

echo "aligning..."
export BOWTIE2_INDEXES=test_result/ref
bowtie2 --phred33 --np 0 --no-unal --local --no-head -x myref -U test_result/stitch/P1_A1.extendedFrags.fastq -S test_result/sam/P1_A1.sam

echo "pick up CpG sites..."
indel.py test_result/sam/P1_A1.sam test_result/sam_indel/P1_A1.site
methsite.py test_result/total_CpG_sites.txt test_result/sam_indel/P1_A1.site test_result/meth_site/P1_A1.site

