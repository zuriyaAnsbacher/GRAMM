#!/usr/bin/env python
import re
import typing
import requests
import sys
import re
from pyard import ARD # NMDP allele codes


# converts CT validation dataset into GL String format to load into graph database

# global XX code GL hash
XX_hash = {}
dna2ser = {}
seroGL = {}


def expand_XX(ac, XX_hash):
	# uses hashing to limit calls to MAC service 
	if (ac in XX_hash):
		gl = XX_hash[ac]
		# print ("Hash Lookup")
	else:
		url = "https://hml.nmdp.org/mac/api/decode?typing="
		response = requests.get(url + ac)
		gl = response.text
		XX_hash[ac] = gl
		# print ("Service Call")

	gl_ars = ard.redux_gl(gl,'lg')
	return(gl_ars)

def expand_AC(ac):
	# print ("AC " + ac)
	gl_ars = ard.redux_gl(ac,'lg')
	return(gl_ars)

def loc_gl (id, loc, typ1, typ2, dna, XX_hash):
	if (typ2 == ""):
		typ2 = typ1
	if ((typ1 == "") or (typ2 == "") or (typ1 == "NEW") or (typ2 == "NEW") or (typ1 is None) or (typ2 is None) or (typ1 is " ") or (typ2 is " ")):
		# print ("SKIP")
		return "SKIP"

	# handle ID 8541096 when DNA = 1 and DQB1*01:XX - invalid allele code
	if ((id == "8541096") and (loc == "DQB1")):
		dna = "0"
		typ1 = "1"
		typ2 = "1"
	if ((id == "9862798") and (loc == "DQB1")):
		dna = "0"
		typ1 = "3" # typing at other locus makes it difficult to switch to serology in the general case	
		typ2 = "1"


	# print (loc + " TYP: " + typ1 + " " + typ2 + " DNA: " + dna)	
	
	# serology
	if (dna == "0"):
		if (loc == "DRB1"):
			loc = "DR"
		if (loc == "DQB1"):
			loc = "DQ"
		if ((typ1 == "X") and (loc != "C")):
		   typ1 = "NULL"
		if ((typ2 == "X") and (loc != "C")):
		   typ2 = "NULL"
		loctyp1 = loc + typ1
		loctyp2 = loc + typ2
		# print ("SERO " + loctyp1 + "" + loctyp2)
		al1 = seroGL[loctyp1]
		al2 = seroGL[loctyp2]
		gl = '+'.join([al1, al2])
		return gl

	# dna typing
	elif (dna == "1"):
		typ1_fields = typ1.split(':')
		subtype1 = typ1_fields[1]
		subtype1_trunc = subtype1[:2] # only take first two characters because of expression chars
		typ2_fields = typ2.split(':')
		subtype2 = typ2_fields[1]
		subtype2_trunc = subtype2[:2] # only take first two characters because of expression chars
		loctyp1 = '*'.join([loc,typ1])
		loctyp2 = '*'.join([loc,typ2])
		if (subtype1 == "XX"):
			al1 = expand_XX(loctyp1,XX_hash)
			# print ("loctyp1 " + loctyp1 + " XX")
		elif (subtype1_trunc.isnumeric()):
			al1 = loctyp1
			# print ("loctyp1 " + loctyp1 + " HR")
		else:
			al1 = expand_AC(loctyp1)
			# print ("loctyp1 " + loctyp1 + " AC")
		if (subtype2 == "XX"):
			al2 = expand_XX(loctyp2,XX_hash)
			# print ("loctyp2 " + loctyp2 + " XX")
		elif (subtype2_trunc.isnumeric()):
			al2 = loctyp2
			# print ("loctyp2 " + loctyp2 + " HR")
		else:
			al2 = expand_AC(loctyp2)
			# print ("loctyp2 " + loctyp2 + " AC")
		# print (loctyp1 + " " + loctyp2 + " " + al1 + " " + al2)
		gl = '+'.join([al1, al2])
		return gl

	else:
		return "SKIP"



##############################################################################
# Function: loadSeroWMDAMap - load WMDA serology
##############################################################################
def loadSeroWMDAMap (ser_map_wmda,ser_gl_wmda,XX_hash):

	wmda_sero_filename = "rel_ser_dna.txt"
	wmda_sero_file = open(wmda_sero_filename, 'r')

	family = {}
	for line in wmda_sero_file:
		(loc, sero, locstar, glstring) = line.split("\t")

		if ((loc != "A") and (loc != "C") and (loc != "B") and (loc != "DR") and (loc != "DQB1")):
			continue


		loc_ser = loc

		if (loc == "DR"):
			loc = "DRB1"
		if (loc_ser == "DQB1"):
			loc_ser = "DQ"

		# skip 3-digit and 4-digit serologic nomenclature
		if (len(sero) == 4):
			continue
		if (len(sero) == 3):
			continue
		if (sero == "0"):
			sero = "NULL"
		if (sero is None):
			sero = "NULL"
		if (sero == "?"):
			sero = "NULL" # C serology
		sero_null = 0 # flag if serology is null
		if (sero == "NULL"):
			sero_null = 1
		sero = loc_ser + sero

		glstring = glstring.strip() # remove line ending
		# print(glstring)
		glstring_array = glstring.split('/')
		# glstring_array[-1] = glstring_array[-1].strip()

		# print(glstring_array)

		ars_dupe_check = {}
		alleles_sero_ars = []
		for allele in glstring_array:
			# print(allele)
			hla_fields = allele.split(':')
			family = hla_fields[0]
			protein = hla_fields[1]
			allele = family + ":" + protein
			if len(hla_fields) == 3:
				synony = hla_fields[2]
				if (re.search('N',synony)):
					allele = allele + "N"
				if (re.search('Q',synony)):
					allele = allele + "Q"
				if (re.search('L',synony)):
					allele = allele + "L"
			if len(hla_fields) == 4:
				noncoding = hla_fields[3]
				if (re.search('N',noncoding)):
						allele = allele + "N"
				if (re.search('Q',noncoding)):
					allele = allele + "Q"
				if (re.search('L',noncoding)):
					allele = allele + "L"

			loctyp = loc + "*" + allele
			ser_map_wmda[loctyp] = sero

			allele_ars = ard.redux_gl(loctyp,'lg')
			# print(loctyp + " " + allele_ars)
			if (sero_null == 1):
				# skip "g" alleles in list of nulls when serology is null
				if (re.search('g',allele_ars)):
					continue
			if (allele_ars in ars_dupe_check):
				alleles_sero_ars.append(allele_ars)
				ars_dupe_check[allele_ars] = 1;

		glstring_ars = '/'.join(alleles_sero_ars);
		if sero in ser_gl_wmda:
			ser_gl_wmda[sero] = glstring_ars + "/" + ser_gl_wmda[sero]
		else:
			ser_gl_wmda[sero] = glstring_ars

	# Missing 3-digit and 4-digit serotypes from WMDA file
	ser_gl_wmda["A2403"] = "A*24:03"
	ser_gl_wmda["DR103"] = "DRB1*01:03"
	ser_gl_wmda["B5102"] = "B*51:02"
	ser_gl_wmda["B4005"] = "B*40:05"
	ser_gl_wmda["A203"] = "A*02:03"
	ser_gl_wmda["B2708"] = "B*27:08"

	# Add C blank serology as "X"
	# TODO - confirm that HapLogic matching algorithm does not consider nulls?
	ser_gl_wmda["CX"] = "/".join([expand_XX("C*12:XX",XX_hash),expand_XX("C*14:XX",XX_hash),expand_XX("C*15:XX",XX_hash),expand_XX("C*16:XX",XX_hash),expand_XX("C*17:XX",XX_hash),expand_XX("C*18:XX",XX_hash)])



	wmda_sero_file.close()
	return (0)


### Main 


# initialize ARD object
ard = ARD()

loadSeroWMDAMap(dna2ser,seroGL,XX_hash)

ct_filename = sys.argv[1]
gl_filename = sys.argv[2]
subject_type = sys.argv[3]

# ct_filename = "./data/recip_dpre.in.txt"
# don_filename = ./data/don.gl.txt"
# pat_filename = ./data/pat.gl.txt"

# print (ct_filename)
# print (gl_filename)
# print (subject_type)


ct_file = open(ct_filename, 'r')
gl_file = open(gl_filename, 'w')

dupeid_checker = {}

for line in ct_file:
	(recip_id, recip_phen_seq, recip_detail_race, recip_ethnicity, recip_broad_race, 
		recip_A_dna, recip_A_1, recip_A_2, 
		recip_B_dna, recip_B_1, recip_B_2, 
		recip_DRB1_dna, recip_DRB1_1, recip_DRB1_2,
		recip_C_dna, recip_C_1, recip_C_2,
		recip_DQB1_dna, recip_DQB1_1, recip_DQB1_2,
		recip_DRB3_dna, recip_DRB3_1, recip_DRB3_2, 
		recip_DRB4_dna, recip_DRB4_1, recip_DRB4_2,
		recip_DRB5_dna, recip_DRB5_1, recip_DRB5_2, 
		donor_id,  donor_detail_race, donor_ethnicity, donor_broad_race, 
		donor_A_dna, donor_A_1, donor_A_2, donor_A_pdtl, 
		donor_B_dna, donor_B_1, donor_B_2, donor_B_pdtl, 
		donor_DRB1_dna, donor_DRB1_1, donor_DRB1_2, donor_DRB1_pdtl,
		donor_C_dna, donor_C_1, donor_C_2, donor_C_pdtl, 
		donor_DQB1_dna, donor_DQB1_1, donor_DQB1_2, donor_DQB1_pdtl,
		donor_DRB3_dna, donor_DRB3_1, donor_DRB3_2, donor_DRB3_pdtl, 
		donor_DRB4_dna, donor_DRB4_1, donor_DRB4_2, donor_DRB4_pdtl, 
		donor_DRB5_dna, donor_DRB5_1, donor_DRB5_2, donor_DRB5_pdtl) = line.split('\t')

	if (subject_type == "R"):
		if recip_id in dupeid_checker:
			continue
		else:
			dupeid_checker[recip_id] = 1
			gl_A = loc_gl(recip_id,"A",recip_A_1,recip_A_2,recip_A_dna,XX_hash)
			gl_C = loc_gl(recip_id,"C",recip_C_1,recip_C_2,recip_C_dna,XX_hash)
			gl_B = loc_gl(recip_id,"B",recip_B_1,recip_B_2,recip_B_dna,XX_hash)
			gl_DRB1 = loc_gl(recip_id,"DRB1",recip_DRB1_1,recip_DRB1_2,recip_DRB1_dna,XX_hash)
			gl_DQB1 = loc_gl(recip_id,"DQB1",recip_DQB1_1,recip_DQB1_2,recip_DQB1_dna,XX_hash)
	else:
		if donor_id in dupeid_checker:
			continue
		else:
			dupeid_checker[donor_id] = 1
			gl_A = loc_gl(donor_id,"A",donor_A_1,donor_A_2,donor_A_dna,XX_hash)
			gl_C = loc_gl(donor_id,"C",donor_C_1,donor_C_2,donor_C_dna,XX_hash)
			gl_B = loc_gl(donor_id,"B",donor_B_1,donor_B_2,donor_B_dna,XX_hash)
			gl_DRB1 = loc_gl(donor_id,"DRB1",donor_DRB1_1,donor_DRB1_2,donor_DRB1_dna,XX_hash)
			gl_DQB1 = loc_gl(donor_id,"DQB1",donor_DQB1_1,donor_DQB1_2,donor_DQB1_dna,XX_hash)

	gl_array = list()
	if (gl_A != "SKIP"):
		gl_array.append(gl_A)
	if (gl_C != "SKIP"):
		gl_array.append(gl_C)
	if (gl_B != "SKIP"):
		gl_array.append(gl_B)
	if (gl_DRB1 != "SKIP"):
		gl_array.append(gl_DRB1)
	if (gl_DQB1 != "SKIP"):
		gl_array.append(gl_DQB1)

	# print (recip_id + " " + subject_type)
	# print (recip_A_1 + " " + recip_A_2)
	# print (gl_A + " " + gl_C + " " + gl_B + " " + gl_DRB1 + " " + gl_DQB1)
	# print (gl_array)
	gl_string = '^'.join(gl_array)
	gl_string.replace("HLA-","")

	if (subject_type == "R"):
		gl_file.write('%'.join([recip_id, gl_string]) + '\n')
	else:
		gl_file.write('%'.join([donor_id, gl_string]) + '\n')

ct_file.close()
gl_file.close()