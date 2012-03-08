#!/frodo/shared/epd/bin/python

import os
import sys
from ConfigParser import SafeConfigParser

analysisM = {'all':['nki_pipeline.py','alff_pipeline.py','rsfc_pipeline.py'], 'basic':['nki_pipeline.py'],
	'nuisance':['nuisance_pipeline.py'],'alff':['alff_pipeline.py'], 'rsfc': ['rsfc_pipeline.py'], 'nuisance+alff' : ['nuisance_pipeline.py','alff_pipeline.py'], 
	'nuisance+alff+rsfc' : ['nuisance_pipeline.py', 'alff_pipeline.py', 'rsfc_pipeline.py'] , 'nuisance+rsfc' :['nuisance_pipeline.py','rsfc_pipeline.py'], 
	'rsfc+alff': ['rsfc_pipeline.py','alff_pipeline.py'], 'basic+alff' : ['nki_pipeline.py','alff_pipeline.py'], 'basic+rsfc' : ['nki_pipeline.py','rsfc_pipeline.py']}

def getValues(value):

	values = (value.strip(' ')).split(',')

	return values
	

def generateDefinition(idx, parsermap):


	cwd = os.getcwd()
	fpath = os.path.join(cwd,'dir_setup_%d.ini' %idx)
	file = open(fpath, 'w')

	print >>file,'[setup]'
	for key, value in parsermap.items():

		key = key.strip(' ')
		if not (key == 'strategy' or key == 'analysis'):

			if idx < len(value):
				value[idx] = value[idx].strip(' ')
				print >>file, key+ ":"+ value[idx]
			else:
				value[len(value) - 1] = value[len(value) - 1].strip(' ')
				print >>file,key+ ":"+ value[ len(value) - 1 ]
			  
	file.close()

	return fpath

def getPipelineName(idx, parsermap):

	pipelines = []

	value = parsermap['analysis']

	if idx < len(value):
		key = value[idx]

	else:
		key = value[len(value) - 1]

	if key in analysisM.keys():
		pipelines = analysisM[key]

	if pipelines == []:

		print 'invalid pipeline option: ' + key
		sys.exit()

	return pipelines 

def runAnalysis(parsermap):

	import subprocess
	strategies = parsermap['strategy']

	
	for i in range(0,len(strategies)):

		file = generateDefinition(i, parsermap)

		strat = strategies[i]
		pipelines = getPipelineName(i,parsermap)

		for pipeline in pipelines:
			cmd = './%s %s %s %s' %(pipeline, sys.argv[2], file, strat)
			
			print cmd
			subprocess.Popen(cmd, shell=True)
			

def readDefinition():

	parsermap = {}
	parser = SafeConfigParser()
	parser.read(sys.argv[1])

	for section in parser.sections():
		for variable , value in parser.items(section):

			parsermap[variable] = getValues(value)
			#print variable + ' -> ' + str(parsermap[variable])

	runAnalysis(parsermap)

def main():

	readDefinition()

if __name__ == "__main__":


	sys.exit(main())
