# usage
#python svf.py feature -b /data/xwu/bam/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam -c /data/xwu/gsd/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz  -v /home/xwu/opsv/NA12878high.vcf.gz -o NA12878high_feature.txt -g hg19 -f /home/xwu/reference/hs37d5/hs37d5.fa
#python svf.py feature -b /data/xwu/bam/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam -c test_train_set.bed  -v /home/xwu/opsv/NA12878high.vcf.gz -o test_NA12878train_feature.txt -g hg19 -f /home/xwu/reference/hs37d5/hs37d5.fa

sample='NA12878'
bam='/data/xwu/bam/NA12878.mapped.ILLUMINA.bwa.CEU.high_coverage_pcr_free.20130906.bam'
snv='/home/xwu/opsv/NA12878high.vcf.gz'
sv='test_train_set.bed'
call='test_test_set.bed'

feature='NA12878/feature.txt'
gtype='hg19'
gfile='/home/xwu/reference/hs37d5/hs37d5.fa'
pop='CEU'
echo choose module [feature train evaluate assess]
read x
case $x in 
	feature) python Sven.py feature -s $sample -b $bam -v $snv  -c $sv;;

	train) python Sven.py train -i $feature;;

	evaluate) python Sven.py evaluate -i $feature;;

	assess) python Sven.py assess -s $sample -b $bam -v $snv -c $call;;

	*) echo module name error!
esac