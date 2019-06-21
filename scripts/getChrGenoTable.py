import sys
from collections import defaultdict
import json
import csv
from flaskext.mysql import MySQL
sys.path.insert(0, "../")
from app import app

mysql = MySQL(app)
db = mysql.connect()

#Usage: python getChrGenoTable.py [chromosome_name] [output_path]

def getSNPData(chromosome,outpath):
	db = mysql.connect()
	cur = db.cursor()
	sql = "SELECT SNPArray.SNPArrayName,Specie.SpecieName,Cultivar.CultivarName,Accession.Accession,SNP.SNPid,SNP.Chromosome,SNP.Reference,SNP.Strand,SNP.Coordination,Accession_has_SNP.Genotype FROM mydb.SNP LEFT join mydb.SNPArray ON SNP.SNPArray_idSNPArray=SNPArray.idSNPArray LEFT join mydb.Accession_has_SNP ON SNP.SNPTableid = Accession_has_SNP.SNP_SNPTableid LEFT join mydb.Accession ON Accession_has_SNP.Accession_idAccession = Accession.idAccession LEFT join mydb.Cultivar ON Accession.Cultivar_idCultivar = Cultivar.idCultivar left join mydb.Specie on Cultivar.Specie_idSpecie = Specie.idSpecie WHERE SNP.Chromosome LIKE"
	sql_cmd = sql + "'" + chromosome + "%'" 
	cur.execute(sql_cmd)
	cultivars = set()
	snpdict = defaultdict(list)
	snpinfodict = defaultdict(list)
	fields = cur.fetchall()
	counter = 0
	for field in fields:
		if field[2] is not None:
			unid = field[2].replace(" ","_") + "_" + field[3]
			cultivars.add(unid)
			SNPid = field[4]
			if SNPid in snpdict:
				snpdict[SNPid].append((unid,field[9]))
			else:
				snpdict[SNPid] = [(unid,field[9])]
			if SNPid in snpinfodict:
				pass	
			else:
				snpinfodict[SNPid] = [(field[8],field[7],field[6],field[5],field[0])] 
#get all unids and add them to each snpdict entry with "?" value
#iterate over SNPid in dict and append missing unid's
	all_cult_of_SNP = set()
	sorted_cult = sorted(list(cultivars))
	for SNPID in snpdict: #for each key in dict
		for SNPcultivars in snpdict[SNPID]:
			all_cult_of_SNP.add(SNPcultivars[0])

		missing_cult = list(cultivars - all_cult_of_SNP)
		all_cult_of_SNP = set()
		for cultivar in missing_cult:
			snpdict[SNPID].append((cultivar,"?"))
	out = outpath + '{}.csv'.format(_chr) 
	header = sorted_cult
	columnlist = ['SNPid','Array','Chr','Ref_SNP','Strand','SNP_Position_On_Chr']
	columnlist.reverse()
	for column in columnlist:
        	header.insert(0,column)
	with open(out, 'w+') as csvhead:
		writer = csv.writer(csvhead, dialect='excel', quoting=csv.QUOTE_NONE)
		writer.writerow(header)

	for SNPID in snpdict: 
		header = sorted_cult
		columnlist = ['SNPid','Array','Chr','Ref_SNP','Strand','SNP_Position_On_Chr']
		columnlist.reverse()
		for column in columnlist:
			header.insert(0,column)

		genolist = list()
		cultlist = list()
		for cult in sorted_cult:
			for SNPcultivars in snpdict[SNPID]:
				if cult == SNPcultivars[0]:
					genolist.append(SNPcultivars[1])
		for value in snpinfodict[SNPID][0][0:]:
			genolist.insert(0,value)
		genolist.insert(0,SNPID)

		with open(out, 'a') as csvfile:
			writer = csv.writer(csvfile, dialect='excel', quoting=csv.QUOTE_NONE)
			writer.writerow(genolist)
	print("Successfully printed genotype table to {}".format(out))
	cur.close()
	db.close()

_chr = sys.argv[1]
outpath = sys.argv[2]
getSNPData(_chr, outpath)
