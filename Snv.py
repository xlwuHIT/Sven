import vcf
import numpy as np

class Snv():
	def __init__(self, vcf_file,sample_name):
		self._name = sample_name
		self._fp=vcf.Reader(filename=vcf_file)
		self._format,self._chr_prefix=self.vcfFormat()


	def snv(self,sv):
		snv_dict={}
		for i in sv._svdict:
			deps=[]
			ratio=[]
			for target in sv._svdict[i]:
				# if target[0]!='14':
				# 	continue
				#print target
				snvs=self._fp.fetch(self._chr_prefix+target[0],int(target[1]),int(target[2]))
				for snv in snvs:
					try:
						sample_genotype=snv.genotype(self._name)
						if sample_genotype.gt_type==None:
							continue
						#print sample_genotype
						deps.append(sample_genotype['DP'])
						if sample_genotype.is_het:
							if self._format=='AD':
								a0=sample_genotype['AD'][-2]
								a1=sample_genotype['AD'][-1]
							else:
								a0=sample_genotype['RO']
								a1=sample_genotype['AO']
							#print sample_genotype
							if a0==0:
								ratio.append(float('inf'))
							else:
								ratio.append(a1*1.0/a0)
					except:
						pass
			snv_dict[i]=(median(deps),len(deps),median(ratio),len(ratio))
		return snv_dict

	def snv2(self,sv):
		snv_dict={}
		for target in sv._svs:
			i=target[-1]
			deps=[]
			ratio=[]
			snvs=self._fp.fetch(self._chr_prefix+target[0],int(target[1]),int(target[2]))
			for snv in snvs:
				try:
					sample_genotype=snv.genotype(self._name)
					if sample_genotype.gt_type==None:
						continue
					#print sample_genotype
					deps.append(sample_genotype['DP'])
					if sample_genotype.is_het:
						if self._format=='AD':
							a0=sample_genotype['AD'][-2]
							a1=sample_genotype['AD'][-1]
						else:
							a0=sample_genotype['RO']
							a1=sample_genotype['AO']
						#print sample_genotype
						if a0==0:
							ratio.append(float('inf'))
						else:
							ratio.append(a1*1.0/a0)
				except:
					pass
			snv_dict[i]=(median(deps),len(deps),median(ratio),len(ratio))
		return snv_dict


	def vcfFormat(self):
		tmp=next(self._fp)
		vcf_format=tmp.FORMAT
		ref_name=tmp.CHROM
		if 'AD' in vcf_format:
			vcf_format='AD'
		elif 'DPR' in vcf_format:
			vcf_format='DPR'
		else:
			print ('unknown vcf format',vcf_file)
			exit(1)
		chr_prefix='chr' if 'chr' in ref_name else ''
		return vcf_format,chr_prefix

def median(l,choice=0):
	if len(l)==0:
		if choice==0:
			return 0
		else:
			return float('nan')
	else:
		return np.nanmedian(list(l))
