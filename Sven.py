#!/home/xwu/bin/python
#-*- coding:utf-8 -*-
import argparse,sys,os
from Sample import *
from multiprocessing import Pool
from Model import *

def feature(args):
	mytime('feature construction start')
	sample=args.sample
	bam=args.bam
	snv=args.snv
	sv=args.sv
	if not os.path.exists(sample):
		os.mkdir(sample)
	feature=sample+'/feature.txt' if args.feature_file==None else sample+'/'+args.feature_file
	gtype=args.genome_type
	gfile=args.genome_file
	pop=args.population

	Sample(sample,bam,snv,sv,feature,gtype,gfile,pop)

	mytime('feature construction end')



def train(args):
	mytime('train start')

	feature=args.feature_file
	suffix=args.model_suffix
	index=feature.rfind('/')
	assessment=feature[:index+1]+'train_assessment.txt' if args.assessment_file==None else args.assessment_file

	Model(feature,suffix,assessment)

	mytime('train end')




def evaluate(args):
	mytime('evaluate start')

	feature=args.feature_file
	suffix=args.model_suffix
	index=feature.rfind('/')
	assessment=feature[:index+1]+'assessment.txt' if args.assessment_file==None else args.assessment_file
	filteration=feature[:index+1]+'filteration.txt' if args.filteration_file==None else args.filteration_file

	Model(feature,suffix,assessment,filteration)

	mytime('evaluate end')




def assess(args):
	mytime('assess start')
	sample=args.sample
	bam=args.bam
	snv=args.snv
	sv=args.sv
	if not os.path.exists(sample):
		os.mkdir(sample)
	feature=sample+'/feature.txt' if args.feature_file==None else sample+'/'+args.feature_file
	gtype=args.genome_type
	gfile=args.genome_file
	pop=args.population

	Sample(sample,bam,snv,sv,feature,gtype,gfile,pop)

	suffix=args.model_suffix
	index=feature.rfind('/')
	assessment=feature[:index+1]+'assessment.txt' if args.assessment_file==None else args.assessment_file
	filteration=feature[:index+1]+'filteration.txt' if args.filteration_file==None else args.filteration_file

	Model(feature,suffix,assessment,filteration)

	mytime('assess end')



my_version='1.0.1'

svenfig='''
	 __________ 
	|   ____   |
	|  |    |__|
	|  |           __       __     ______     __________
	|  |_______   |  |     |  |   / ____ \   /   ____   \ 
	|_______   |  |  |     |  |  / /____\ \  |  |    |  |
	 __     |  |  \  \     /  /  | _______|  |  |    |  |
	|  |    |  |   \  \   /  /   | |     _   |  |    |  |
	|  |____|  |    \  \_/  /    \ \____/ |  |  |    |  |
	|__________|     \_____/      \______/   |__|    |__|

	Structural Variation Evaluation and Filteration Tool
	
'''
my_usage=svenfig+'''
1.python Sven.py feature *
2.python Sven.py train **
3.python Sven.py test ***
4.python Sven.py assess ****
'''
my_description=svenfig+'''
[S]tructural [v]ariation [e]valuation and filteratio[n] tool(Sven) is a toolkit
for scoring and filtering SVs in structural variation detection sets produced 
by various SV detection tools.
Sven extracts useful information from mapping data, snv data and reference genome
to construct several features, train a deletion model and a duplication model,
assess SVs from structural variation detection set and finally produces a set of
SVs with high-confidence.

Version: 	v%s
Author: 	Xiaoliang Wu
Contact: 	xlwu@hit.edu.cn
Github: 	https://github.com/xlwuHIT/Sven 
'''%(my_version)

my_epilog='''
Sven includes 4 modules as follows:

1.feature module:
$Sven feature [-h] -s SAMPLE -b BAM -v SNV -c SV [-e FEATURE_FILE]
              [-t {hg18,hg19,hg38}] [-f GENOME_FILE] [-p POPULATION]

2.train module:
$Sven train [-h] -i FEATURE_FILE [-m MODEL_SUFFIX] [-a ASSESSMENT_FILE]

3.evaluate module:
$Sven evaluate [-h] -i FEATURE_FILE [-m MODEL_SUFFIX] 
			   [-a ASSESSMENT_FILE] [-l FILTER_FILE]

4.assess module: 
$Sven assess [-h] -s SAMPLE -b BAM -v SNV -c SV [-e FEATURE_FILE]
             [-t {hg18,hg19,hg38}] [-f GENOME_FILE] [-p POPULATION]
             [-m MODEL_SUFFIX] [-a ASSESSMENT_FILE] [-l FILTER_FILE]
	
	<Input "Sven [module name] -h" for more details>
	
'''








parser=argparse.ArgumentParser(prog='Sven',description=my_description,epilog=my_epilog,formatter_class=argparse.RawDescriptionHelpFormatter)
sub_parsers=parser.add_subparsers()
parser.add_argument('--version', '-v', 
		action = 'version', 
		version = '%(prog)s {version}'.format(version=my_version))

#feature
feature_parser=sub_parsers.add_parser('feature',help='construct features of one sample')
feature_parser.add_argument('-s','--sample',required=True,help='sample name corresponding bam and snv file')
feature_parser.add_argument('-b','--bam',required=True,help='mapping file(*.bam)')
feature_parser.add_argument('-v','--snv',required=True,help='snv file(*.vcf.gz) with index file(*.vcf.gz.tbi)')
feature_parser.add_argument('-c','--sv',required=True,help='structural variation calling file(*.vcf.gz or *.bed)')
feature_parser.add_argument('-e','--feature_file',help='output feature file')
feature_parser.add_argument('-t','--genome_type',required=False,default='hg19',choices=['hg18','hg19','hg38'],help='version of reference genome file(hg18 or hg19 or hg38)')
feature_parser.add_argument('-f','--genome_file',required=False,default='hs37d5.fa',help='reference genome file(*.fa)')
feature_parser.add_argument('-p','--population',required=False,default='CEU',help='sample population group code')
feature_parser.set_defaults(function=feature)

#train
train_parser=sub_parsers.add_parser('train',help='train and dump deletion and duplication model')
train_parser.add_argument('-i','--feature_file',required=True,help='produced feature file(*.txt)')
train_parser.add_argument('-m','--model_suffix',required=False,default='model.pkl',help='deletion and duplication model suffix for dumping model')
train_parser.add_argument('-a','--assessment_file',help='assessment file for train set')
train_parser.set_defaults(function=train)

#evaluate
evaluate_parser=sub_parsers.add_parser('evaluate',help='evaluate and filter SVs according to their features')
evaluate_parser.add_argument('-i','--feature_file',required=True,help='produced feature file(*.txt)')
evaluate_parser.add_argument('-m','--model_suffix',required=False,default='model.pkl',help='deletion and duplication model suffix for loading model')
evaluate_parser.add_argument('-a','--assessment_file',help='assessment file for structural variation detection set')
evaluate_parser.add_argument('-l','--filteration_file',help='final filteration file for structural variation detection set')
evaluate_parser.set_defaults(function=evaluate)


#assess
assess_parser=sub_parsers.add_parser('assess',help='sequential execution of feature and evaluate process with pre-trained models')
assess_parser.add_argument('-s','--sample',required=True,help='sample name corresponding bam and snv file')
assess_parser.add_argument('-b','--bam',required=True,help='mapping file(*.bam)')
assess_parser.add_argument('-v','--snv',required=True,help='snv file(*.vcf.gz) with index file(*.vcf.gz.tbi)')
assess_parser.add_argument('-c','--sv',required=True,help='structural variation calling file(*.vcf.gz or *.bed)')
assess_parser.add_argument('-e','--feature_file',help='output feature file')
assess_parser.add_argument('-t','--genome_type',required=False,default='hg19',choices=['hg18','hg19','hg38'],help='version of reference genome file(hg18 or hg19 or hg38)')
assess_parser.add_argument('-f','--genome_file',required=False,default='hs37d5.fa',help='reference genome file(*.fa)')
assess_parser.add_argument('-p','--population',required=False,default='CEU',help='sample population group code')
assess_parser.add_argument('-m','--model_suffix',required=False,default='model.pkl',help='deletion and duplication model suffix for loading model')
assess_parser.add_argument('-a','--assessment_file',help='assessment file for structural variation detection set')
assess_parser.add_argument('-l','--filteration_file',help='final filteration file for structural variation detection set')
assess_parser.set_defaults(function=assess)


args=parser.parse_args()
args.function(args)
