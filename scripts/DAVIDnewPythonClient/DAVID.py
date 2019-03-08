#import ssl
#ssl._create_default_https_context = ssl._create_unverified_context

import sys, getopt
# reload(sys)
# sys.setdefaultencoding('utf8')
sys.path.append('../')
import os
import logging
import traceback as tb
import suds.metrics as metrics
from tests import *
from suds import *
from suds.client import Client
from datetime import datetime
import csv

#
# input and output files.
#
#https://www.tutorialspoint.com/python/python_command_line_arguments.html
inputfile = ''
outputDIR = '../'
argv = sys.argv[1:]
try:
  opts, args = getopt.getopt(argv,"hi:o:d:t:",["iFILE=","oDIR=","DAVID=","idType="])
except getopt.GetoptError:
  print 'python2 chartReport.py -i <inputfile> -o <outputDIR> -d <DAVID version> -t <Identifier>'
  sys.exit(2)
for opt, arg in opts:
  if opt == '-h':
     print 'python2 chartReport.py -i <inputfile> -o <outputDIR> -d <DAVID version> -t <Identifier>'
     sys.exit()
  elif opt == '--help':
     print 'python2 chartReport.py -i <inputfile> -o <outputDIR> -d <DAVID version> -t <Identifier>'
     sys.exit()
  elif opt in ("-i", "--iFILE"):
     inputfile = arg
  elif opt in ("-o", "--oDIR"):
    if os.path.isdir(arg):  
      print("--oDIR is a directory")
      outputDIR = arg  
    else :  
      print("--oDIR is a normal file -> directory is used") 
      outputDIR = os.path.dirname(arg)
      print 'Output Dir is :', outputDIR
    
  elif opt in ("-d", "--DAVID"):
     DAVID = arg     
  elif opt in ("-t", "--idType"):
     idT = arg    
     # print client.service.getConversionTypes()
### call 
# python2 DAVID.py -i geneList.csv -o "./testRES" -d DAVID68 -t "ENSEMBL_GENE_ID"

# /home/adsvy/GitHubRepo/SnakeWF_HIF/viper/scripts/DAVIDnewPythonClient
# python2 DAVID.py -i geneList.csv -o /home/adsvy/GitHubRepo/SnakeWF_HIF/viper/scripts/DAVIDnewPythonClient/TEST/a.txt -d DAVID68 -t "ENSEMBL_GENE_ID"
# python2 DAVID.py -i /home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/salmonAlignment/estcount_HSQ-vs-NSQsi_edegR_Res_sig_0.05_MYlog2FC_1.csv -o /home/adsvy/GitHubRepo/SnakeWF_HIF/viper/scripts/DAVIDnewPythonClient/TEST/ -d DAVID68 -t "ENSEMBL_GENE_ID"

print 'Input file is :', inputfile
print 'Output Dir is :', outputDIR
print 'DAVID version is :', DAVID    
print 'idType is :', idT    

# print 'Number of arguments:', len(sys.argv), 'arguments.'
# print 'Argument List:', str(sys.argv)


# #
# # Configuring Logging
# #
# errors = 0
# # setup_logging()
# #logging.getLogger('suds.client').setLevel(logging.DEBUG)

# #create logger
# logger = logging.getLogger('suds.client')
# logger.setLevel(logging.DEBUG)

# # create console handler and set level to debug
# ch = logging.StreamHandler()
# ch.setLevel(logging.DEBUG)

# # create formatter
# formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# # add formatter to ch
# ch.setFormatter(formatter)

# # add ch to logger
# logger.addHandler(ch)

# # 'application' code
# # logger.debug('debug message')
# # logger.info('info message')
# logger.warn('warn message')
# logger.error('error message')
# logger.critical('critical message')


#
# create a service client using the wsdl.
#
if DAVID == "DAVID68":
    print("DAVID 6.8")
    url = 'https://david.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl' ### DAVID 6.8
    client = Client(url)
    client.wsdl.services[0].setlocation('https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/') ### DAVID 6.8
elif DAVID == "DAVID67":
    print("DAVID 6.7")
    url = 'https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl' ### DAVID 6.7 (with ILLUMINA_ID)
    client = Client(url)
    client.wsdl.services[0].setlocation('https://david-d.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap11Endpoint/') ### DAVID 6.7 (with ILLUMINA_ID)
else:
    print '--DAVID/-d must be "DAVID68" or "DAVID67"'
    sys.exit(2)

print 'url=%s' % url

#
# print the service (introspection)
#
#print client

### TODO add commandline parameter for user mail 

#authenticate user email  
userMail = 'claus.weinholdt@informatik.uni-halle.de'
if client.service.authenticate(userMail):
  print "authenticated user: ", userMail
else :
  sys.exit(2)

#check authenticate user email in browser
#https://david.ncifcrf.gov/webservice/services/DAVIDWebService/authenticate?args0=claus.weinholdt@informatik.uni-halle.de


# print client.service.getAllListNames()
# print client.service.getAllPopulationNames()
# print client.service.getAllAnnotationCategoryNames()
# print client.service.getConversionTypes()

#
# read gene list csv
#

# inputfile = "/home/adsvy/GitHubRepo/SnakeWF_HIF/results/deg/edegR/hg38/PE/salmonAlignment/estcount_NSQ-vs-NSQsi_edegR_Res_sig_0.05_MYlog2FC_1.csv"
# inputfile = "/home/adsvy/GitHubRepo/SnakeWF_HIF/viper/scripts/DAVIDnewPythonClient/geneList.csv"
ResCSV = False
with open(inputfile, 'rb') as csvfile:
  CSVreader = csv.reader(csvfile, delimiter=';', quotechar='"')
  tmp_list = list()
  for row in CSVreader:
    if len(row) > 1 :
      # print ''.join(str(row[1]))
      tmp_list.append(row[1])         #second column with geneID !!!
      ResCSV = True
    else:     
      # print ','.join(row)
      tmp_list.append(''.join(row))

if(ResCSV):
  inputIds = inputIds =  ','.join(tmp_list[1:])
else:
  inputIds =  ','.join(tmp_list)

# print inputIds
print 'Number of genes:', len(tmp_list)

# inputIds = '1112_g_at,1331_s_at,1355_g_at,1372_at,1391_s_at,1403_s_at,1419_g_at,1575_at,1645_at,1786_at,1855_at,1890_at,1901_s_at,1910_s_at,1937_at,1974_s_at,1983_at,2090_i_at,31506_s_at,31512_at,31525_s_at,31576_at,31621_s_at,31687_f_at,31715_at,31793_at,31987_at,32010_at,32073_at,32084_at,32148_at,32163_f_at,32250_at,32279_at,32407_f_at,32413_at,32418_at,32439_at,32469_at,32680_at,32717_at,33027_at,33077_at,33080_s_at,33246_at,33284_at,33293_at,33371_s_at,33516_at,33530_at,33684_at,33685_at,33922_at,33963_at,33979_at,34012_at,34233_i_at,34249_at,34436_at,34453_at,34467_g_at,34529_at,34539_at,34546_at,34577_at,34606_s_at,34618_at,34623_at,34629_at,34636_at,34703_f_at,34720_at,34902_at,34972_s_at,35038_at,35069_at,35090_g_at,35091_at,35121_at,35169_at,35213_at,35367_at,35373_at,35439_at,35566_f_at,35595_at,35648_at,35896_at,35903_at,35915_at,35956_s_at,35996_at,36234_at,36317_at,36328_at,36378_at,36421_at,36436_at,36479_at,36696_at,36703_at,36713_at,36766_at,37061_at,37096_at,37097_at,37105_at,37166_at,37172_at,37408_at,37454_at,37711_at,37814_g_at,37898_r_at,37905_r_at,37953_s_at,37954_at,37968_at,37983_at,38103_at,38128_at,38201_at,38229_at,38236_at,38482_at,38508_s_at,38604_at,38646_s_at,38674_at,38691_s_at,38816_at,38926_at,38945_at,38948_at,39094_at,39187_at,39198_s_at,39469_s_at,39511_at,39698_at,39908_at,40058_s_at,40089_at,40186_at,40271_at,40294_at,40317_at,40350_at,40553_at,40735_at,40790_at,40959_at,41113_at,41280_r_at,41489_at,41703_r_at,606_at,679_at,822_s_at,919_at,936_s_at,966_at'
# inputIds = 'ENSG00000100644,ENSG00000071967,ENSG00000176171,ENSG00000140263,ENSG00000235162,ENSG00000162433,ENSG00000143847,ENSG00000152256,ENSG00000164176,ENSG00000119471,ENSG00000099256,ENSG00000179958,ENSG00000156017,ENSG00000104419,ENSG00000074410,ENSG00000126261,ENSG00000104765,ENSG00000159399,ENSG00000107159,ENSG00000156531,ENSG00000147526,ENSG00000135074,ENSG00000159167,ENSG00000164332,ENSG00000171388,ENSG00000111674,ENSG00000213625,ENSG00000169756,ENSG00000173786,ENSG00000003989,ENSG00000123728,ENSG00000148110,ENSG00000046604,ENSG00000170379,ENSG00000178695,ENSG00000163516,ENSG00000167601,ENSG00000247095,ENSG00000106484,ENSG00000164168,ENSG00000172292,ENSG00000151806,ENSG00000101236,ENSG00000122884,ENSG00000102144,ENSG00000077713,ENSG00000185787,ENSG00000111845,ENSG00000003436,ENSG00000136167,ENSG00000166295,ENSG00000185515,ENSG00000136261,ENSG00000024526,ENSG00000002834,ENSG00000117862,ENSG00000066084,ENSG00000113658,ENSG00000104369,ENSG00000135622,ENSG00000131941,ENSG00000104635,ENSG00000145604,ENSG00000141425,ENSG00000134333,ENSG00000177917,ENSG00000114268,ENSG00000180011,ENSG00000142156,ENSG00000151233,ENSG00000139117,ENSG00000164603,ENSG00000140743,ENSG00000151882,ENSG00000170961,ENSG00000141401,ENSG00000196950,ENSG00000259479,ENSG00000197860,ENSG00000181458,ENSG00000225205,ENSG00000170949,ENSG00000171204,ENSG00000178462,ENSG00000105877,ENSG00000109107,ENSG00000139508,ENSG00000146674,ENSG00000147852,ENSG00000141562,ENSG00000061656,ENSG00000253882,ENSG00000159860,ENSG00000113083,ENSG00000141552,ENSG00000249992,ENSG00000064225,ENSG00000123570,ENSG00000198729,ENSG00000077585,ENSG00000227617,ENSG00000139832,ENSG00000160233,ENSG00000120910,ENSG00000118804,ENSG00000225206,ENSG00000186352,ENSG00000174804,ENSG00000065534,ENSG00000107738,ENSG00000114698,ENSG00000149506,ENSG00000255737,ENSG00000147036,ENSG00000272887,ENSG00000123572,ENSG00000225791,ENSG00000197358,ENSG00000143590,ENSG00000163898,ENSG00000173890,ENSG00000164283,ENSG00000272870,ENSG00000116711,ENSG00000242193,ENSG00000133195,ENSG00000183018,ENSG00000175764,ENSG00000047346,ENSG00000172638,ENSG00000015568,ENSG00000128805,ENSG00000106031,ENSG00000282339,ENSG00000135678,ENSG00000272079,ENSG00000142623,ENSG00000100311,ENSG00000167619,ENSG00000146205,ENSG00000231298,ENSG00000131094,ENSG00000278266,ENSG00000277701,ENSG00000139636,ENSG00000255641,ENSG00000204682,ENSG00000150051,ENSG00000226415,ENSG00000249786,ENSG00000259330,ENSG00000231890,ENSG00000253161,ENSG00000044524,ENSG00000130303,ENSG00000106302,ENSG00000254204,ENSG00000167131,ENSG00000100167,ENSG00000197632,ENSG00000167772,ENSG00000164932,ENSG00000164849,ENSG00000154310,ENSG00000274049,ENSG00000179869,ENSG00000277027,ENSG00000254665,ENSG00000213453,ENSG00000143320,ENSG00000261534,ENSG00000267034,ENSG00000280800,ENSG00000130829,ENSG00000253741,ENSG00000202198,ENSG00000280614,ENSG00000281181,ENSG00000189149,ENSG00000079385,ENSG00000088726,ENSG00000171408,ENSG00000172967,ENSG00000254806,ENSG00000144834,ENSG00000258940,ENSG00000263535,ENSG00000188305'

idType = str(idT) #'ENSEMBL_GENE_ID' #'AFFYMETRIX_3PRIME_IVT_ID'
listName = 'make_up'
listType = 0
print 'percent of found genes:',client.service.addList(inputIds, idType, listName, listType)

# print client.service.getDefaultCategoryNames()
# setCategories
categoryString = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE,GOTERM_BP_DIRECT,GOTERM_CC_DIRECT,GOTERM_MF_DIRECT'))
print categoryString

# As of October 2016 there has been a new GO category added (GO Direct), which provides GO mappings directly annotated by the source database (no parent terms included). I'm not sure why DAVID selects this as the default but it has very low specificity.
# BP_1 to BP_5 are the levels of GO terms with increasing specificity (1 is least,5 is most). BP_ALL is all of the levels combined..
# GO FAT filters out very broad GO terms based on a measured specificity of each term, so it's like filtered BP_ALL.
# New GO category (GO FAT) filters out very broad GO terms based on a measured specificity of each term (not level-specificity)
# New GO category (GO Direct) provides GO mappings directly annotated by the source database (no parent terms included)



#
# getChartReport
#
print "run getChartReport"

threshold = 0.1 #set to 1, to get the full result -- || -- EASE Score default 0.1, a modified Fisher Exact P-Value https://david-d.ncifcrf.gov/helps/functional_annotation.html#fisher
count = 2   # The threshold of minimum gene counts belonging to an annotation term. It has to be equal or greater than 0. Default is 2. In short, you do not trust the term only having one gene involved
chartReport = client.service.getChartReport(threshold,count)
chartRow = len(chartReport)
print 'Total chart records:',chartRow

#parse and print chartReport
resF = outputDIR+"/"+DAVID+"_"+'chartReport.txt'
with open(resF, 'w') as fOut:
  fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\taFDR\trFDR\tFisherExact\n')
  for simpleChartRecord in chartReport:
            categoryName = simpleChartRecord.categoryName
            termName = simpleChartRecord.termName
            listHits = simpleChartRecord.listHits
            percent = simpleChartRecord.percent
            ease = simpleChartRecord.ease
            Genes = simpleChartRecord.geneIds
            listTotals = simpleChartRecord.listTotals
            popHits = simpleChartRecord.popHits
            popTotals = simpleChartRecord.popTotals
            foldEnrichment = simpleChartRecord.foldEnrichment
            bonferroni = simpleChartRecord.bonferroni
            benjamini = simpleChartRecord.benjamini
            aFDR = simpleChartRecord.afdr
            rFDR = simpleChartRecord.rfdr
            Fisher = simpleChartRecord.fisher
            rowList = [categoryName,termName,str(listHits),str(percent),str(ease),Genes,str(listTotals),str(popHits),str(popTotals),str(foldEnrichment),str(bonferroni),str(benjamini),str(aFDR),str(rFDR),str(Fisher)]
            fOut.write('\t'.join(rowList).encode('utf-8').strip()+'\n')
print 'write file:', resF, 'finished!'

threshold = 1 #set to 1, to get the full result -- || -- EASE Score default 0.1, a modified Fisher Exact P-Value https://david-d.ncifcrf.gov/helps/functional_annotation.html#fisher
count = 2   # The threshold of minimum gene counts belonging to an annotation term. It has to be equal or greater than 0. Default is 2. In short, you do not trust the term only having one gene involved
chartReport = client.service.getChartReport(threshold,count)
chartRow = len(chartReport)
print 'Total chart records:',chartRow

#parse and print chartReport
resFt = outputDIR+"/"+DAVID+"_"+'chartReport_T1.txt'
with open(resFt, 'w') as fOut:
	fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\taFDR\trFDR\tFisherExact\n')
	for simpleChartRecord in chartReport:
            categoryName = simpleChartRecord.categoryName
            termName = simpleChartRecord.termName
            listHits = simpleChartRecord.listHits
            percent = simpleChartRecord.percent
            ease = simpleChartRecord.ease
            Genes = simpleChartRecord.geneIds
            listTotals = simpleChartRecord.listTotals
            popHits = simpleChartRecord.popHits
            popTotals = simpleChartRecord.popTotals
            foldEnrichment = simpleChartRecord.foldEnrichment
            bonferroni = simpleChartRecord.bonferroni
            benjamini = simpleChartRecord.benjamini
            aFDR = simpleChartRecord.afdr
            rFDR = simpleChartRecord.rfdr
            Fisher = simpleChartRecord.fisher
            rowList = [categoryName,termName,str(listHits),str(percent),str(ease),Genes,str(listTotals),str(popHits),str(popTotals),str(foldEnrichment),str(bonferroni),str(benjamini),str(aFDR),str(rFDR),str(Fisher)]
            fOut.write('\t'.join(rowList).encode('utf-8').strip()+'\n')
print 'write file:', resFt, 'finished!'

# print chartReport


#
# getGeneClusterReport
#
print "run getGeneClusterReport"


if DAVID == "DAVID68":
  ### Classification Stringency => Medium
  overlap = 3 # Similarity Term Overlap
  linkage = 0.5 # Similarity Threshold
  initialSeed = 3 # Initial Group Membership
  finalSeed = 3 # Final Group Membership
  kappa = 50 #Multiple Linkage Threshold

elif DAVID == "DAVID67":
  ### David 6.7
  # ### Classification Stringency => Low
  # overlap = 4 # Similarity Term Overlap
  # linkage = 0.3 # Similarity Threshold
  # initialSeed = 3 # Initial Group Membership
  # finalSeed = 3 # Final Group Membership
  # kappa = 50 #Multiple Linkage Threshold

  ### Classification Stringency => Medium
  overlap = 4 # Similarity Term Overlap
  linkage = 0.35 # Similarity Threshold
  initialSeed = 4 # Initial Group Membership
  finalSeed = 4 # Final Group Membership
  kappa = 50 #Multiple Linkage Threshold

  # ### Classification Stringency => high
  # overlap = 4 # Similarity Term Overlap
  # linkage = 0.4 # Similarity Threshold
  # initialSeed = 5 # Initial Group Membership
  # finalSeed = 5 # Final Group Membership
  # kappa = 50 #Multiple Linkage Threshold


geneClusterReport = client.service.getGeneClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)

totalClusters = len(geneClusterReport)
print 'Total groups:',totalClusters
resF2 = outputDIR+"/"+DAVID+"_"+'geneClusterReport.txt'

print 'type:', type(geneClusterReport[0])
ttt = str(type(geneClusterReport[0]))
if "NoneType" in ttt:
# if str(type(geneClusterReport[0])) == 'NoneType':
  print 'NoneType'
  with open(resF2, 'w') as fOut:
    fOut.write('NULL\n')
  # sys.exit(2)
else:
  # resF2 = 'list1.geneClusterReport.txt'
  # print geneClusterReport
  with open(resF2, 'w') as fOut:
      for simpleGeneClusterRecord in geneClusterReport:
          #EnrichmentScore = simpleGeneClusterRecord.score
          fOut.write( simpleGeneClusterRecord.name + '\tEnrichmentScore: ' + str(simpleGeneClusterRecord.score) + '\n')
          fOut.write(idType+'\tGene Name\n')
          for listRecord in simpleGeneClusterRecord.listRecords:
              gene_id = ','.join(listRecord.values)
              fOut.write(gene_id + '\t' + listRecord.name + '\n')   
              
  print 'write file:', resF2, 'finished!'
# sys.exit(2)

#
# getTermClusteringReport
#
print "run getTermClusteringReport"


overlap=3
initialSeed = 3
finalSeed = 3
linkage = 0.5
kappa = 50
termClusteringReport = client.service.getTermClusterReport(overlap, initialSeed, finalSeed, linkage, kappa)

#parse and print report
totalRows = len(termClusteringReport)
print 'Total clusters:',totalRows
# resF3 = 'list1.termClusteringReport.txt'
resF3 = outputDIR+"/"+DAVID+"_"+'termClusteringReport.txt'

print 'type:', type(termClusteringReport[0])
ttt = str(type(termClusteringReport[0]))
if "NoneType" in ttt:
  print 'NoneType'
  with open(resF3, 'w') as fOut:
    fOut.write('NULL\n')
  # sys.exit(2)
else:
  with open(resF3, 'w') as fOut:
    i = 0
    for simpleTermClusterRecord in termClusteringReport:
      i = i+1
      EnrichmentScore = simpleTermClusterRecord.score
      fOut.write('Annotation Cluster '+str(i) + '\tEnrichmentScore:'+str(EnrichmentScore)+'\n')
      # fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\tFDR\n')
      fOut.write('Category\tTerm\tCount\t%\tPvalue\tGenes\tList Total\tPop Hits\tPop Total\tFold Enrichment\tBonferroni\tBenjamini\taFDR\trFDR\tFisherExact\n')
      for simpleChartRecord in  simpleTermClusterRecord.simpleChartRecords:
        categoryName = simpleChartRecord.categoryName
        termName = simpleChartRecord.termName
        listHits = simpleChartRecord.listHits
        percent = simpleChartRecord.percent
        ease = simpleChartRecord.ease
        Genes = simpleChartRecord.geneIds
        listTotals = simpleChartRecord.listTotals
        popHits = simpleChartRecord.popHits
        popTotals = simpleChartRecord.popTotals
        foldEnrichment = simpleChartRecord.foldEnrichment
        bonferroni = simpleChartRecord.bonferroni
        benjamini = simpleChartRecord.benjamini
        # FDR = simpleChartRecord.afdr
        # rowList = [categoryName,termName,str(listHits),str(percent),str(ease),Genes,str(listTotals),str(popHits),str(popTotals),str(foldEnrichment),str(bonferroni),str(benjamini),str(FDR)]
        aFDR = simpleChartRecord.afdr
        rFDR = simpleChartRecord.rfdr
        Fisher = simpleChartRecord.fisher
        rowList = [categoryName,termName,str(listHits),str(percent),str(ease),Genes,str(listTotals),str(popHits),str(popTotals),str(foldEnrichment),str(bonferroni),str(benjamini),str(aFDR),str(rFDR),str(Fisher)]
        fOut.write('\t'.join(rowList)+'\n')
  print 'write file:', resF3, 'finished!'



#
# getTableReport
#
print "run getTableReport"

# #print client.service.getDefaultCategoryNames()
# categoryString = str(client.service.setCategories('BBID,BIOCARTA,COG_ONTOLOGY,GOTERM_BP_FAT,GOTERM_CC_FAT,GOTERM_MF_FAT,INTERPRO,KEGG_PATHWAY,OMIM_DISEASE,PIR_SUPERFAMILY,SMART,SP_PIR_KEYWORDS,UP_SEQ_FEATURE'))
categories = categoryString.split(',')

tableReport = client.service.getTableReport()
tableRow = len(tableReport)
print 'Total table records:',tableRow
# resF4 = 'list1.tableReport.txt'
resF4 = outputDIR+"/"+DAVID+"_"+'tableReport.txt'

print 'type:', type(tableReport[0])
ttt = str(type(tableReport[0]))
if "NoneType" in ttt:
  print 'NoneType'
  with open(resF4, 'w') as fOut:
    fOut.write('NULL\n')
  # sys.exit(2)
else:
  with open(resF4, 'w') as fOut:
    categoryConcat = '\t'.join(categories);
    fOut.write('ID\tGene Name\tSpecies\t'+categoryConcat) 
    for tableRecord in tableReport:
      name = tableRecord.name
      species = tableRecord.species
      for arrayString in tableRecord.values:
        gene_id = ','.join(x for x in arrayString.array)
        rowList = [gene_id,name,species]
        fOut.write('\n'+'\t'.join(rowList))
      for annotationRecord in tableRecord.annotationRecords:
        default_value = ''
        category_dict = dict.fromkeys(categories,default_value)     
        termsConcat = '';
      for term in annotationRecord.terms:
        termString = term.split("$")[1]
        termList = [termString,termsConcat]
        termsConcat = ','.join(termList)    
        category_dict[str(annotationRecord.category)] = termsConcat;
      for key in category_dict:
        fOut.write('\t'+category_dict[key])

  print 'write file:', resF4, 'finished!'



