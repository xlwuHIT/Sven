#!/home/xwu/bin/python
#-*- coding:utf-8 -*-

from sklearn.linear_model import LogisticRegression
from sklearn.externals import joblib
from sklearn.model_selection import train_test_split
import numpy as np
import os

class Model():
	def __init__(self,filename,suffix,assessment,filteration=None,coefs_file='coefs.txt',):
		self._filename=filename
		self._suffix=suffix
		self._assessment=assessment
		self._filteration=filteration
		self._coefs_file=coefs_file

		data=self.loadFeature()
		if data[0]=='train':
			self.Train(data)
		else:
			self.Test(data)


	def Train(self,data):
		del_model = LogisticRegression(C=100000.0, penalty='l2', tol=0.0001,solver='liblinear')
		dup_model = LogisticRegression(C=100000.0, penalty='l2', tol=0.0001,solver='liblinear')
		del_model.fit(data[1][0],data[1][1])
		dup_model.fit(data[2][0],data[2][1])
		print ('The training score of deletion group is',del_model.score(data[1][0],data[1][1]))
		print ('The training score of duplication group is',dup_model.score(data[2][0],data[2][1]))
		joblib.dump(del_model,'del_'+self._suffix)
		joblib.dump(dup_model,'dup_'+self._suffix)
		#coefs
		model_dict={'del_model':del_model,'dup_model':dup_model}
		f=open(self._coefs_file,'w')
		f.write('#\t'+'\t'.join(data[-1][0][7:])+'\tintercept\n')
		for m in ['del_model','dup_model']:
			f.write(m+'\t')
			#scale=abs(model_dict[m].intercept_)
			coefs=model_dict[m].coef_[0]
			intercept=model_dict[m].intercept_
			for i in range(len(coefs)):
				f.write(str(coefs[i])+'\t')
			f.write(str(intercept)+'\n')
		f.close()
		#new data file
		del_dpp=np.array([del_model.decision_function(data[1][0]),
				np.max(del_model.predict_proba(data[1][0]),1),
				del_model.predict(data[1][0])]).T
		dup_dpp=np.array([dup_model.decision_function(data[2][0]),
				np.max(dup_model.predict_proba(data[2][0]),1),
				dup_model.predict(data[2][0])]).T
		dpp=np.vstack((del_dpp,dup_dpp)).astype(str)
		titled_dpp=np.vstack((np.array([['z-value','p-value','prediction']]),dpp))
		assessment_data=np.hstack((data[-1],titled_dpp))
		np.savetxt(self._assessment,assessment_data,fmt='%s',delimiter='\t')


	
	def Test(self,data):
		if os.path.exists('del_'+self._suffix) and os.path.exists('dup_'+self._suffix):
			del_model=joblib.load('del_'+self._suffix)
			dup_model=joblib.load('dup_'+self._suffix)
			#new data file
			del_dpp=np.array([del_model.decision_function(data[1]),
					np.max(del_model.predict_proba(data[1]),1),
					del_model.predict(data[1])]).T
			dup_dpp=np.array([dup_model.decision_function(data[2]),
					np.max(dup_model.predict_proba(data[2]),1),
					dup_model.predict(data[2])]).T
			dpp=np.vstack((del_dpp,dup_dpp)).astype(str)
			titled_dpp=np.vstack((np.array([['z-value','p-value','prediction']]),dpp))
			assessment_data=np.hstack((data[-1],titled_dpp))
			np.savetxt(self._assessment,assessment_data,fmt='%s',delimiter='\t')
			
			col_name=assessment_data[0]
			targets=assessment_data[1:]
			filteration_targets=np.array([l for l in targets if float(l[17])==1])
			filteration_data=np.vstack((col_name,filteration_targets))[:,(0,1,2,3,4,16)]
			np.savetxt(self._filteration,filteration_data,fmt='%s',delimiter='\t')

		else:
			print ('model is missing!')

	def loadFeature(self):
		data=np.loadtxt(self._filename,skiprows=0,dtype=str,delimiter='\t',)
		col_name=data[0]
		targets=data[1:]
		del_group=np.array([l for l in targets if l[4]=='DEL'])
		dup_group=np.array([l for l in targets if l[4]=='DUP'])
		group=np.vstack((del_group,dup_group))
		group=np.vstack((col_name,group))
		if 'category' in col_name:
			del_y=del_group[:,6].astype(float)
			del_x=del_group[:,7:].astype(float)
			dup_y=dup_group[:,6].astype(float)
			dup_x=dup_group[:,7:].astype(float)
			return 'train',(del_x,del_y),(dup_x,dup_y),group
		else:
			del_x=del_group[:,5:].astype(float)
			dup_x=dup_group[:,5:].astype(float)
			return 'test',del_x,dup_x,group












