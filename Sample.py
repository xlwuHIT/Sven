from Bam import *
from Snv import *
from SV import *
import datetime

class Sample():
	def __init__(self,sample_name,bam_file,vcf_file,call_file,output_file,gt='hg19',gf='/home/xwu/reference/hs37d5/hs37d5.fa',population='CEU'):
		
		self._name=sample_name
		self._gt=gt
		self._gf=gf
		self._pop=population

		print ('ref sv')
		self._sv_ref=SV('config/'+gt+'.bed',gt=gt)
		#print self._sv_ref._svdict
		print ('bam')
		self._bam=Bam(bam_file)
		print ('snv')
		self._snv=Snv(vcf_file,sample_name)

		#print 'gc'
		# self._gc=self.getGC('PCRFREE-DOC')
		# #print self._gc

		self._ref={}
		print ('SL')
		self._ref['isize'],self._ref['length']=self._bam.SL2(self._sv_ref)#ci & length
		#print 'rpkm1'
		#self._ref['rpkm']=self._bam.rpkm(self._sv_ref)# rpkm
		#print 'coverage1'
		#self._ref['coverage']=self._bam.coverage(self._sv_ref)# coverage
		print ('snvcov1')
		self._ref['snvcov']=self._snv.snv2(self._sv_ref)# snv_cov nsnv ratio nhet

		print ('sample sv')
		self._sv_sample=SV(call_file,sample_name,gt,gf,sample_flag=True)
		self._sample={}
		print ('rpkm2')
		self._sample['rpkm']=self._bam.rpkm2(self._sv_sample)
		print ('discordant')
		self._sample['discordant']=self._bam.discordant2(self._sv_sample,median(self._ref['isize'].values()))
		print ('snvcov2')
		self._sample['snvcov']=self._snv.snv2(self._sv_sample)
		print ('tagcnv')
		self._sample['tagcnv']=self.tagCnvPercent2()
		

		if len(self._sv_sample._svs[0].fields)<7:
			print ('write test')
			self.write_features_test2(output_file)
		else:
			print ('write train')
			self.write_features_train2(output_file)


	def getGC(self,choice):
		gc={}
		with open('config/GCfactor.txt','r') as fp:
			for l in fp:
				l=l.strip().split('\t')
				if l[0]==choice:
					gc[int(l[1])]=float(l[2])
		return gc

	def tagCnvPercent(self):
		tagcnv={}
		taglist=[]
		with open('config/tagSNP_'+self._gt+'.txt') as f:
			next(f)
			for l in f:
				l=l.strip().split('\t')
				chr=l[1][3:]
				start=int(l[2])
				end=int(l[3])
				pos=int(l[4])
				r2=float(l[5])
				pop=l[6]
				if pop==self._pop:
					snvs=self._snv._fp.fetch(self._snv._chr_prefix+chr,pos-1,pos)
					for snv in snvs:
						sample_genotype=snv.genotype(self._name)
						if sample_genotype.is_het:
							taglist.append((chr,start,end,r2))
		
		for i in self._sv_sample._svdict:
			c=self._sv_sample._svs[int(i)][0]
			s=int(self._sv_sample._svs[int(i)][1])
			e=int(self._sv_sample._svs[int(i)][2])
			r2=0
			for j in taglist:
				if c==j[0]:
					if e<j[1] or s>j[2]:
						continue
					else:
						overlap=min(e,j[2])-max(s,j[1])
						ratio=overlap*1.0/(e-s)
						if ratio>0.5 and j[3]>r2:
							r2=j[3]
			tagcnv[i]=r2
		return tagcnv
		#print tagcnv

	def tagCnvPercent2(self):
		tagcnv={}
		taglist=[]
		with open('config/tagSNP_'+self._gt+'.txt') as f:
			next(f)
			for l in f:
				l=l.strip().split('\t')
				chr=l[1][3:]
				start=int(l[2])
				end=int(l[3])
				pos=int(l[4])
				r2=float(l[5])
				pop=l[6]
				if pop==self._pop:
					snvs=self._snv._fp.fetch(self._snv._chr_prefix+chr,pos-1,pos)
					for snv in snvs:
						sample_genotype=snv.genotype(self._name)
						if sample_genotype.is_het:
							taglist.append((chr,start,end,r2))
		
		for target in self._sv_sample._svs:
			c=target[0]
			s=int(target[1])
			e=int(target[2])
			r2=0
			for j in taglist:
				if c==j[0]:
					if e<j[1] or s>j[2]:
						continue
					else:
						overlap=min(e,j[2])-max(s,j[1])
						ratio=overlap*1.0/(e-s)
						if ratio>0.5 and j[3]>r2:
							r2=j[3]
			tagcnv[target[-1]]=r2
		return tagcnv


	def write_features_test(self,output_file):
		#features=['coverage','GCcoverage','ndiscordant','nsplit','nconcordant','nsnvs','snv_coverage','nhets','allele_ratio']
		op=open(output_file,'w')
		op.write('chromosome\tstart\tend\tlength\tgenotype\trpkm\tndiscordant\tnsplit\tnconcordant\tsnv_coverage\tnsnvs\tallele_ratio\tnhets\tRLCR\ttagCNV\n')
		for i in self._sv_sample._svdict:
			s=self._sv_sample._svs[int(i)][0]+'\t'+self._sv_sample._svs[int(i)][1]+'\t'+self._sv_sample._svs[int(i)][2]
			s+='\t'+str(int(self._sv_sample._svs[int(i)][2])-int(self._sv_sample._svs[int(i)][1]))+'\t'+self._sv_sample._svs[int(i)][3]
			#s+='\t'+str(self._sample['rpkm'][i])+'\t'+str(self._sample['rpkm'][i]*self._gc[self._sv_sample._GC[i]])
			s+='\t'+str(self._sample['rpkm'][i])
			for j in range(3):
				s+='\t'+str(self._sample['discordant'][i][j])
			for j in range(4):
				s+='\t'+str(self._sample['snvcov'][i][j])
			s+='\t'+str(self._sv_sample._RLCR[i])
			s+='\t'+str(self._sample['tagcnv'][i])
			op.write(s+'\n')
		op.close()

	def write_features_train(self,output_file):
		#features=['coverage','GCcoverage','ndiscordant','nsplit','nconcordant','nsnvs','snv_coverage','nhets','allele_ratio']
		op=open(output_file,'w')
		op.write('chromosome\tstart\tend\tlength\tgenotype\tsample\tcategory\trpkm\tndiscordant\tnsplit\tnconcordant\tsnv_coverage\tnsnvs\tallele_ratio\tnhets\tRLCR\ttagCNV\n')
		for i in self._sv_sample._svdict:
			s=self._sv_sample._svs[int(i)][0]+'\t'+self._sv_sample._svs[int(i)][1]+'\t'+self._sv_sample._svs[int(i)][2]
			s+='\t'+str(int(self._sv_sample._svs[int(i)][2])-int(self._sv_sample._svs[int(i)][1]))+'\t'+self._sv_sample._svs[int(i)][3]+'\t'+self._sv_sample._svs[int(i)][4]+'\t'+self._sv_sample._svs[int(i)][5]
			#s+='\t'+str(self._sample['rpkm'][i])+'\t'+str(self._sample['rpkm'][i]*self._gc[self._sv_sample._GC[i]])
			s+='\t'+str(self._sample['rpkm'][i])
			for j in range(3):
				s+='\t'+str(self._sample['discordant'][i][j])
			for j in range(4):
				s+='\t'+str(self._sample['snvcov'][i][j])
			s+='\t'+str(self._sv_sample._RLCR[i])
			s+='\t'+str(self._sample['tagcnv'][i])
			op.write(s+'\n')
		op.close()

	def write_features_test2(self,output_file):
		#features=['coverage','GCcoverage','ndiscordant','nsplit','nconcordant','nsnvs','snv_coverage','nhets','allele_ratio']
		op=open(output_file,'w')
		op.write('chromosome\tstart\tend\tlength\tgenotype\trpkm\tndiscordant\tnsplit\tnconcordant\tsnv_coverage\tnsnvs\tallele_ratio\tnhets\tRLCR\ttagCNV\n')
		for target in self._sv_sample._svs:
			i=target[-1]
			s=self._sv_sample._svs[int(i)][0]+'\t'+self._sv_sample._svs[int(i)][1]+'\t'+self._sv_sample._svs[int(i)][2]
			s+='\t'+str(int(self._sv_sample._svs[int(i)][2])-int(self._sv_sample._svs[int(i)][1]))+'\t'+self._sv_sample._svs[int(i)][3]
			#s+='\t'+str(self._sample['rpkm'][i])+'\t'+str(self._sample['rpkm'][i]*self._gc[self._sv_sample._GC[i]])
			s+='\t'+str(self._sample['rpkm'][i])
			for j in range(3):
				s+='\t'+str(self._sample['discordant'][i][j])
			for j in range(4):
				s+='\t'+str(self._sample['snvcov'][i][j])
			s+='\t'+str(self._sv_sample._RLCR[i])
			s+='\t'+str(self._sample['tagcnv'][i])
			op.write(s+'\n')
		op.close()

	def write_features_train2(self,output_file):
		#features=['coverage','GCcoverage','ndiscordant','nsplit','nconcordant','nsnvs','snv_coverage','nhets','allele_ratio']
		op=open(output_file,'w')
		op.write('chromosome\tstart\tend\tlength\tgenotype\tsample\tcategory\trpkm\tndiscordant\tnsplit\tnconcordant\tsnv_coverage\tnsnvs\tallele_ratio\tnhets\tRLCR\ttagCNV\n')
		for target in self._sv_sample._svs:
			i=target[-1]
			s=self._sv_sample._svs[int(i)][0]+'\t'+self._sv_sample._svs[int(i)][1]+'\t'+self._sv_sample._svs[int(i)][2]
			s+='\t'+str(int(self._sv_sample._svs[int(i)][2])-int(self._sv_sample._svs[int(i)][1]))+'\t'+self._sv_sample._svs[int(i)][3]+'\t'+self._sv_sample._svs[int(i)][4]+'\t'+self._sv_sample._svs[int(i)][5]
			#s+='\t'+str(self._sample['rpkm'][i])+'\t'+str(self._sample['rpkm'][i]*self._gc[self._sv_sample._GC[i]])
			s+='\t'+str(self._sample['rpkm'][i])
			for j in range(3):
				s+='\t'+str(self._sample['discordant'][i][j])
			for j in range(4):
				s+='\t'+str(self._sample['snvcov'][i][j])
			s+='\t'+str(self._sv_sample._RLCR[i])
			s+='\t'+str(self._sample['tagcnv'][i])
			op.write(s+'\n')
		op.close()



def mytime(content,logFile='Sven_log.txt'):
	fp=open(logFile,'a')
	now=datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	fp.write(now+'\t'+content+'\n')
	fp.close()
