"""
Author: kzhu
Date: 09/25/2017
Generate N random patient mutation profiles
"""
import random
import sys
import getopt
import copy
import numpy as np
import math
import multiprocessing

#colors = ["DEL", "AMP", "SNV", "EXPROUT"]
"""
NODE_MAP: map each gene in the network to the connected component
"""
node_map = dict()
"""
GENES: gene names; PATIENTS: patient names
"""
genes = []
patients = []
"""
MUTATIONS: binary relation {(gene/supernode, color)}
"""
mutations = []
mutations_super = []
"""
Statistics
"""
supernode_sizes = []
maxavg5_supernode_sizes = []
color_counts = []
mutations_color_counts = []

def addNodeGroup(s0, s1, group_id, cc_id):
	mutations_super[patients.index(s0)][0].append(group_id)
	mutations_super[patients.index(s0)][1].append(cc_id)			
	mutations_super[patients.index(s0)][2].append([s1])

def readProfile(profile_name):
	fp = open(profile_name, 'r')
	for line in fp:
        	line = line.rstrip()
		s = line.split()
		if s[1] in genes and s[0] not in patients:
			patients.append(s[0])
	fp.close()
	for i in range (0, len(patients)):
		mutations.append([[], []])
		mutations_color_counts.append(dict())
		mutations_super.append([[], [], []])
		supernode_sizes.append([])
		color_counts.append(dict())
	fp = open(profile_name, 'r')
	gid = 0
	for line in fp:
        	line = line.rstrip()
		s = line.split()
		if s[1] in genes:
			pidx = patients.index(s[0])
			if s[1] in mutations_color_counts[pidx].keys():
				mutations_color_counts[pidx][s[1]] += 1
			else:
				mutations_color_counts[pidx][s[1]] = 1
			mutations[pidx][0].append(s[1])
			mutations[pidx][1].append(s[2])
			if s[2] in color_counts[pidx].keys():
				color_counts[pidx][s[2]] += 1
			else:
				color_counts[pidx][s[2]] = 1
			"""
			Group genes on the same chromosome arm as supernodes
			"""
			if node_map[s[1]] not in mutations_super[pidx][1]:
				addNodeGroup(s[0], s[1], gid, node_map[s[1]])
				gid += 1
			else:
				gidx = mutations_super[pidx][1].index(node_map[s[1]])
				if s[1] not in mutations_super[pidx][2][gidx]:
					mutations_super[pidx][2][gidx].append(s[1])
	fp.close()
	supernode_stat = [0, 0, 0]
	for i in range (0, len(patients)):
		for j in range (0, len(mutations_super[i][0])):
			supernode_sizes[i].append(len(mutations_super[i][2][j]))
		maxavg5_supernode_sizes.append(sum(sorted(supernode_sizes[i])[::-1][:5]) / 5.0)
		supernode_stat[0] += max(supernode_sizes[i])
		supernode_stat[1] += len(supernode_sizes[i])
		supernode_stat[2] += len(mutations[i][0])
	print "AVG MAX SUPERNODE SIZE"
	print supernode_stat[0] * 1.0 / len(patients)
	print "AVG NUM SUPERNODES"
	print supernode_stat[1] * 1.0 / len(patients)
	print "AVG NUM NODES PER SUPERNODE"
	print supernode_stat[2] * 1.0 / supernode_stat[1]

def readNetworkPartitionProfile(partition_file_name):
	fp = open(partition_file_name, 'r')
	for line in fp:
       		line = line.rstrip()
		s = line.split()
		genes.append(s[0])
	fp.close()
	num_cc = 0
	max_cc_idx = -1
	fp = open(partition_file_name, 'r')
	for line in fp:
       		line = line.rstrip()
		s = line.split()
		node_map[s[0]] = int(s[2])
		if int(s[2]) > max_cc_idx:
			max_cc_idx = int(s[2])
		num_cc += 1.0 / int(s[1])
	fp.close()
	if int(round(num_cc)) != max_cc_idx + 1:
		raise Exception("Failed to shuffle the colors.")

"""
Random shuffle strategy: preserving the num nodes with i = 1, 2, 3, 4 colors
Used for p-value lower bound estimation
"""
def duplRemShuffle(gene_color_count, color_list):
	gene_list_length = len(color_list) 
	gene_list_ = gene_color_count.keys()
	random.shuffle(gene_list_)
	multicolor_count = [v for v in gene_color_count.values() if v > 1]
	multicolor_count.sort(reverse = True)
	random.shuffle(color_list)
	mui = [[], []]
	for v in multicolor_count:
		vc_set = set([])
		for i in range (0, len(color_list)):
			vc_set.add(color_list[i])
			if len(vc_set) >= v:
				break
		vc_set = list(vc_set)
		for j in range (0, v):
			mui[0].append(gene_list_[0])
			mui[1].append(vc_set[j])
			color_list.remove(vc_set[j])
		gene_list_ = gene_list_[1:]
		random.shuffle(color_list)
	mui[0] += gene_list_
	mui[1] += color_list
	if len(mui[0]) != gene_list_length or len(mui[1]) != gene_list_length:
		raise Exception("Failed to shuffle the colors.")
	return mui

"""
Return TRUE if there is a gene(node) with i >= 2 colors for given patient
FALSE otherwise
"""
def duplicationTest(gene_list, color_list):
	if len(set(gene_list)) == len(gene_list):
		return False
	return True

"""
Shuffle the supernodes randomly for given patient
"""
def shuffleSupernodes(mus, sns):
	for i in range (0, len(patients)):
		if len(mus[i][0]) != len(sns[i]) or len(mus[i][1]) != len(sns[i]) or \
			len(mus[i][2]) != len(sns[i]):
			raise Exception("Failed to sort the supernodes.")
		muidx = [j for j in range (0, len(sns[i]))]
		random.shuffle(muidx)
		mus[i][0] = [mus[i][0][j] for j in muidx]
		mus[i][1] = [mus[i][1][j] for j in muidx]
		mus[i][2] = [mus[i][2][j] for j in muidx]
	return mus

def randomAssignPriorityColors(musi, cri, sigmai):
	mui = [[], []]
	epsilon = 0.05
	priority_colors = [0, 0]
	total = sum(cri.values())		
	for c in cri.keys():
		if c == "AMP":
			priority_colors[0] += cri[c] * (1.0 + epsilon) / total
		if c == "DEL":
			priority_colors[1] += cri[c] * (1.0 + epsilon) / total
	while True:
		cri_ = copy.deepcopy(cri)
		mui = [[], []]
		for j in range (0, len(musi[0])):
			r = random.random()
			length = len(musi[2][j])
			if r < priority_colors[0]:
				for k in range (0, length):
					mui[0].append(musi[2][j][k])
					mui[1].append("AMP")
				cri_["AMP"] -= length
			elif r < priority_colors[0] + priority_colors[1]:
				for k in range (0, length):
					mui[0].append(musi[2][j][k])
					mui[1].append("DEL")
				cri_["DEL"] -= length
		valid = True
		if ("AMP" in cri.keys() and abs(cri_["AMP"]) > sigmai) or ("DEL" in cri.keys() and \
			abs(cri_["DEL"]) > sigmai):
			valid = False
		if valid:
			return mui, cri_
	return mui, cri

def greedyAssignPriorityColors(musi, cri):
	mui = [[], []]
	priority_colors = [c for c in cri.keys() if c == "AMP" or c == "DEL"]
	j = 0
	for c in priority_colors:
		while True:
			if j >= len(musi[0]):
				break
			length = len(musi[2][j])
			if cri[c] - length > 0:
				for k in range (0, length):
					mui[0].append(musi[2][j][k])
					mui[1].append(c)
				cri[c] -= length
				j += 1
			else:
				r = random.random()
				if cri[c] - length == 0:
					r = 0
				if r < cri[c] * 1.0 / length:
					for k in range (0, length):
						mui[0].append(musi[2][j][k])
						mui[1].append(c)
					cri[c] -= length
					j += 1
				break
	return mui, cri

"""
Random shuffle strategy: for random order of supernodes, first assign color AMP and DEL
			then break down supernodes and assign the other colors
Used for p-value upper bound estimation
"""
def greedyShuffleSupernodeColors(mus, cr, mode__):
	mu = []
	for i in range (0, len(patients)):
		"""
		First assign AMP and DEL
		"""
		mui = []
		if mode__ == 2:
			mui, cr[i] = randomAssignPriorityColors(mus[i], cr[i], maxavg5_supernode_sizes[i])
		else:
			mui, cr[i] = greedyAssignPriorityColors(mus[i], cr[i])
		"""
		Assign multi-colored nodes
		"""
		multicolor_count = [v for v in mutations_color_counts[i].values() if v > 2]
		mtg = []
		if len(mui[0]) >= len(multicolor_count):
			mtg = random.sample(mui[0], len(multicolor_count))
		else:
			mtg = mui[0]
		for gene in mtg:
			mui[0] += [gene, gene]
			mui[1] += ["SNV", "EXPROUT"]
			cr[i]["SNV"] -= 1
			cr[i]["EXPROUT"] -= 1
		"""
		Redistribute and then assign the remaining colors
		"""
		remaining_colors = [c for c in cr[i].keys() if cr[i][c] > 0 and (c == "SNV" or c == "EXPROUT")]
		remaining_genes = [gene for gene in mutations_color_counts[i].keys() if gene not in mui[0]]
		random.shuffle(remaining_genes)
		total = sum(cr[i].values()) # The number of remaining colored nodes
		if total < len(remaining_genes):
			raise Exception("Failed to shuffle the colors.")
		if len(remaining_colors) == 1:
			mtg_candidates = [gene for gene in mui[0] if gene not in mtg]
			if total - len(remaining_genes) > len(mtg_candidates):
				print [remaining_colors[0]], total - len(remaining_genes), len(mtg_candidates)
				remaining_colors = ["SNV", "EXPROUT"]
				if "SNV" not in cr[i].keys():
					cr[i]["SNV"] = 0
			else:
				mtg = random.sample(mtg_candidates, total - len(remaining_genes))
				mui[0] += remaining_genes
				mui[1] += [remaining_colors[0] for j in range (0, len(remaining_genes))]
				for gene in mtg:
					mui[0].append(gene)
					mui[1].append(remaining_colors[0])
				cr[i][remaining_colors[0]] -= total
		if len(remaining_colors) > 1:
			s = total * 1.0 / sum([cr[i][c] for c in remaining_colors])
			nexprout = [gene for gene in mutations_color_counts[i].keys() if gene not in \
				[mui[0][j] for j in range (0, len(mui[0])) if mui[1][j] == "EXPROUT"]]
			acr_exprout = min(int(math.floor(cr[i]["EXPROUT"] * s)), len(nexprout))
			acr_snv = total - acr_exprout
			cr[i]["SNV"] -= acr_snv
			cr[i]["EXPROUT"] -= acr_exprout
			single_colors = random.sample(["SNV" for j in range (0, acr_snv)] + \
					["EXPROUT" for j in range (0, acr_exprout)], len(remaining_genes))
			mui[0] += remaining_genes
			mui[1] += single_colors
			acr_snv -= single_colors.count("SNV")
			acr_exprout -= single_colors.count("EXPROUT")
			nsnv = [gene for gene in mutations_color_counts[i].keys() if gene not in \
				[mui[0][j] for j in range (0, len(mui[0])) if mui[1][j] == "SNV"]]
			nexprout = [gene for gene in mutations_color_counts[i].keys() if gene not in \
				[mui[0][j] for j in range (0, len(mui[0])) if mui[1][j] == "EXPROUT"]]
			mui[0] += (random.sample(nsnv, acr_snv) + random.sample(nexprout, acr_exprout))
			mui[1] += (["SNV" for j in range (0, acr_snv)] + ["EXPROUT" for j in range (0, acr_exprout)])
		mu.append(mui)
	return mu, cr

def shuffleProfile(N, b):
	for n in range (0, N):
		mutations_ = copy.deepcopy(mutations)
		for i in range (0, len(patients)):
			"""
			if duplicationTest(mutations_[i][0], mutations_[i][1]):
				print i
			"""
			if not duplicationTest(mutations_[i][0], mutations_[i][1]):
				random.shuffle(mutations_[i][1])
			else:
				mutations_[i] = duplRemShuffle(mutations_color_counts[i], mutations_[i][1])
		fp = open("sim" + str(n + b) + ".tsv", 'w')
		for i in range (0, len(patients)):
			for j in range (0, len(mutations_[i][0])):
				fp.write('%s\t%s\t%s\n' %(patients[i], mutations_[i][0][j], mutations_[i][1][j]))
		fp.close()

def shuffleSuperProfile(N, b, mode_):
	n = 0
	while n < N:
		mutations_super_ = copy.deepcopy(mutations_super)
		"""
		Random order of supernodes
		"""
		mutations_super_ = shuffleSupernodes(mutations_super_, supernode_sizes)
		color_remaining = copy.deepcopy(color_counts)
		try:
			mutations_, color_error = greedyShuffleSupernodeColors(mutations_super_, \
						color_remaining, mode_)
		except Exception:
			continue
		fp = open("sim_" + str(n + b) + ".tsv", 'w')
		for i in range (0, len(patients)):
			for j in range (0, len(mutations_[i][0])):
				fp.write('%s\t%s\t%s\n' %(patients[i], mutations_[i][0][j], mutations_[i][1][j]))
		fp.close()
		fp = open("color_difference" + str(n + b) + ".tsv", 'w')
		for i in range (0, len(patients)):
			for key, value in color_error[i].items():
				fp.write('%s\t%d\t' %(key, value))
			fp.write('\n')
		fp.close()
		n += 1

if __name__ == "__main__":
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'i:p:n:t:luU')
	except getopt.error:
		print "Illegal Input."
		exit()

	filename = ""
	partition_filename = ""
	N = 1000
	"""
	mode = 0: p-value lower bound; mode = 1, 2: p-value upper bound; 
	"""
	mode = 0
	num_threads = 10
	for (op, opVal) in optlist:
		if op == "-i":
			filename = opVal
		if op == "-p":
			partition_filename = opVal
		if op == "-n":
			N = int(opVal)
		if op == "-t":
			num_threads = int(opVal)
		if op == "-l":
			mode = 0
		if op == "-u":
			mode = 1
		if op == "-U":
			mode = 2
	if N % num_threads != 0:
		print "Illegal Input."
		exit()
	readNetworkPartitionProfile(partition_filename)
	readProfile(filename)
	jobs = []
	N = N / num_threads
	if mode == 0:
		for i in range(0, num_threads):
			p = multiprocessing.Process(target = shuffleProfile, args = (N, i * N))
			jobs.append(p)
			p.start()
	elif mode == 1:
		for i in range(0, num_threads):
			p = multiprocessing.Process(target = shuffleSuperProfile, args = (N, i * N, mode))
			jobs.append(p)
			p.start()
	else:
		for i in range(0, num_threads):
			p = multiprocessing.Process(target = shuffleSuperProfile, args = (N, i * N, mode))
			jobs.append(p)
			p.start()
