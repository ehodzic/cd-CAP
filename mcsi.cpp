#include <cstdio>
#include <cstring>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <math.h>
#include <queue>
#include <algorithm>
#include <vector>
#include <utility>
using namespace std;

typedef unsigned long long llu;

typedef std::pair<int, int> colored_node;
typedef std::set<colored_node> colored_node_set;

/* Fancy Output Functions. */
void printHeader(const char * text) {
	int n = strlen(text);
	fprintf(stderr, "\n");
	for (int i = 0; i < n + 4; i++)
		fprintf(stderr, "*");
	fprintf(stderr, "\n* %s *\n", text);
	for (int i = 0; i < n + 4; i++)
		fprintf(stderr, "*");
	fprintf(stderr, "\n");
}

/* Global variables. */
struct stringPairHash {
	size_t operator()(const pair<string, string> & Q) const {
		size_t h1 = hash<string>()(Q.first);
		size_t h2 = hash<string>()(Q.second);
		return h1 ^ (h2 << 1);
	}
};

struct colored_node_set_hash {
	size_t operator()(const colored_node_set &candidate) const {
		size_t hashvalue = 0;
		for (auto node : candidate)
			hashvalue += node.first;
		return hashvalue;
	}
};

struct Entry {
	string * names;
	unordered_map<string, int> indices;
} samples, genes, alterations;

struct PatientBitmask {
	int maxSize = 0;
	int size = 0;
	int len = 0;
	llu * bits = NULL;

	PatientBitmask() {
	}

	PatientBitmask(int maxSize) :
			maxSize(maxSize) {
		size = 0;
		len = ceil(maxSize / 64.0);
		bits = new llu[len];
		memset(bits, 0, sizeof(llu) * len);
	}

	PatientBitmask(PatientBitmask *Q) {
		maxSize = Q->maxSize;
		size = Q->size;
		len = Q->len;
		bits = new llu[len];
		for (int i = 0; i < len; i++)
			bits[i] = Q->bits[i];
	}

	~PatientBitmask() { 
		delete[] bits;
	}

	void setBit(int pos, bool val) {
		int idx = pos / 64;
		if (idx >= len) {
			fprintf(stderr,
					"< Error > Cannot assign bit because bitmask is too small. pos: %d | len: %d | idx: %d | maxSize: %d\n",
					pos, len, idx, maxSize);
			exit(0);
		}
		int bitIdx = pos % 64;
		bool oldVal = getBit(pos);
		if (oldVal ^ val) {	// the bit is about to get changed
			bits[idx] ^= llu(1) << bitIdx;
			if (val)
				size++;
			else
				size--;
		}
	}

	bool getBit(int pos) const {
		int idx = pos / 64;
		if (idx >= len) {
			fprintf(stderr,
					"< Error > Cannot return bit because bitmask is too small.\n");
			exit(0);
		}
		int bitIdx = pos % 64;
		return bits[idx] & (llu(1) << bitIdx);
	}

	int getSize() const {
		return size;
	}

	int getPositionOfFirstSetBit() const {
		if (size == 0) {
			fprintf(stderr,
					"< Error > Cannot return position of first set bit because bitmask is empty.\n");
			exit(0);
		}
		for (int i = 0; i < len; i++) {
			if (bits[i])
				return i * 64 + __builtin_ctzll(bits[i]);
		}
	}

	void mergeBitmask(PatientBitmask *Q) {
		size = 0;
		for (int i = 0; i < len; i++) {
			bits[i] = bits[i] & Q->bits[i];
			size += __builtin_popcountll(bits[i]);
		}
	}
};

struct Graph {
	int V, E;
	int * NSize;
	int ** N;
	string * nodeNames;
	unordered_map<int, int> * idxInNeighbourList;// idxInNeighbourList[a][b] = c means that "In node a's neighbor list, node b is at position c."
	unordered_map<int, int> * incomingEdges;// incomingEdges[a][b] = c means that "There is an incoming edge into a from b's edge list at position c"
	unordered_map<string, int> nodeIndices;
	int numCC;		// number of connected components.
	int * ccIndex;	// connected component index for each node.
	int * ccSize;	// size of each connected component
} G;

unordered_map<int, unsigned int> * geneAlterations;	// geneAlterations[i][j] = c means that "gene i has color c in patient j". The colors are bitmasks (so supporting max 32 different alterations).

void findConnectedComponents() {
	G.numCC = 0;
	G.ccIndex = new int[G.V];
	memset(G.ccIndex, -1, sizeof(G.ccIndex[0]) * G.V);
	int * nodeStack = new int[G.V];
	int nodeStackSize = 0;
	for (int i = 0; i < G.V; i++) {
		if (G.ccIndex[i] == -1) {
			G.ccIndex[i] = G.numCC++;
			nodeStack[nodeStackSize++] = i;
			while (nodeStackSize) {
				int node = nodeStack[--nodeStackSize];
				for (int j = 0; j < G.NSize[node]; j++) {
					int const & neighbor = G.N[node][j];
					if (G.ccIndex[neighbor] == -1) {
						G.ccIndex[neighbor] = G.ccIndex[node];
						nodeStack[nodeStackSize++] = neighbor;
					}
				}
			}
		}
	}
	G.ccSize = new int[G.numCC];
	memset(G.ccSize, 0, sizeof(G.ccSize[0]) * G.numCC);
	for (int i = 0; i < G.V; i++)
		G.ccSize[G.ccIndex[i]]++;
	delete[] nodeStack;
	fprintf(stderr, "\tInput network contains %d connected components.\n",
			G.numCC);
}

/* Reads the input -n parameter as a collection of undirected edges (pairs of node names, separated by whitespace) and stores the information in global Graph object G. */
void readUndirectedNetwork(const char * filename) {
	fprintf(stderr, "Reading the network... ");
	int timerStart = clock();
	unordered_set<pair<string, string>, stringPairHash> uniqueEdges;
	char u[1000], v[1000];
	FILE * fin = fopen(filename, "r");
	G.V = 0;
	while (fscanf(fin, "%s%s", u, v) == 2) {
		char * a = u, *b = v;
		if (strcmp(a, b) > 0)
			swap(a, b);
		else if (strcmp(a, b) == 0)
			continue;
		pair < string, string > e = make_pair(string(a), string(b));
		if (uniqueEdges.find(e) == uniqueEdges.end()) {	// Previously NOT seen undirected edge
			uniqueEdges.insert(e);
			if (G.nodeIndices.find(e.first) == G.nodeIndices.end()) {
				G.nodeIndices[e.first] = G.V;
				G.V++;
			}
			if (G.nodeIndices.find(e.second) == G.nodeIndices.end()) {
				G.nodeIndices[e.second] = G.V;
				G.V++;
			}
		}
	}
	fclose(fin);
	G.E = uniqueEdges.size() * 2;
	G.NSize = new int[G.V];
	G.N = new int *[G.V];
	G.nodeNames = new string[G.V];
	memset(G.NSize, 0, sizeof(G.NSize[0]) * G.V);
	for (auto e : uniqueEdges) {
		int idx1 = G.nodeIndices[e.first];
		if (G.nodeNames[idx1].size() == 0)
			G.nodeNames[idx1] = e.first;

		int idx2 = G.nodeIndices[e.second];
		if (G.nodeNames[idx2].size() == 0)
			G.nodeNames[idx2] = e.second;

		G.NSize[idx1]++;
		G.NSize[idx2]++;
	}
	for (int i = 0; i < G.V; i++) {
		G.N[i] = new int[G.NSize[i]];
	}

	int * tempNSize = new int[G.V];
	memset(tempNSize, 0, sizeof(tempNSize[0]) * G.V);
	G.idxInNeighbourList = new unordered_map<int, int> [G.V];
	G.incomingEdges = new unordered_map<int, int> [G.V];
	for (auto e : uniqueEdges) {
		int idx1 = G.nodeIndices[e.first];
		int idx2 = G.nodeIndices[e.second];
		G.N[idx1][tempNSize[idx1]] = idx2;
		G.idxInNeighbourList[idx1][idx2] = tempNSize[idx1];
		G.incomingEdges[idx2][idx1] = tempNSize[idx1];
		G.N[idx2][tempNSize[idx2]] = idx1;
		G.idxInNeighbourList[idx2][idx1] = tempNSize[idx2];
		G.incomingEdges[idx1][idx2] = tempNSize[idx2];
		tempNSize[idx1]++;
		tempNSize[idx2]++;
	}
	delete[] tempNSize;
	fprintf(stderr, "done. (%.2lf seconds)\n",
			double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf(stderr,
			"\tInput network contains %d nodes and %d undirected edges.\n", G.V,
			G.E / 2);
}

/* Reads the input -l parameter as a collection of "sample gene alterationType" triples, separated by whitespace. */
void readAlterationProfiles(const char * filename, int init_) {
	fprintf(stderr, "Reading the alteration profiles... ");
	int timerStart = clock();
	char sample[1000], gene[1000], alterationType[1000];
	FILE * fin = NULL;
	if (!(fin = fopen(filename, "r"))) {
		fprintf(stderr,
				"\n< Error > Cannot open file '%s'. Please make sure the file exists.\n",
				filename);
		exit(0);
	}
	while (fscanf(fin, "%s%s%s", sample, gene, alterationType) == 3) {
		if (G.nodeIndices.find(string(gene)) == G.nodeIndices.end())
			continue;
		if (samples.indices.find(string(sample)) == samples.indices.end()) {
			int idx = samples.indices.size();
			samples.indices[string(sample)] = idx;
		}
		if (genes.indices.find(string(gene)) == genes.indices.end()) {
			int idx = genes.indices.size();
			genes.indices[string(gene)] = idx;
		}
		if (alterations.indices.find(string(alterationType))
				== alterations.indices.end()) {
			int idx = alterations.indices.size();
			alterations.indices[string(alterationType)] = idx;
		}
	}
	samples.names = new string[samples.indices.size()];
	genes.names = new string[genes.indices.size()];
	for (auto it : samples.indices)
		samples.names[it.second] = it.first;
	for (auto it : genes.indices)
		genes.names[it.second] = it.first;
	if (init_ == 0) {
		alterations.names = new string[alterations.indices.size()];
		for (auto it : alterations.indices)
			alterations.names[it.second] = it.first;
	}
	rewind(fin);
	geneAlterations = new unordered_map<int, unsigned int> [G.V];
	while (fscanf(fin, "%s%s%s", sample, gene, alterationType) == 3) {
		if (G.nodeIndices.find(string(gene)) == G.nodeIndices.end())
			continue;
		int geneIndex = G.nodeIndices[string(gene)];
		int sampleIndex = samples.indices[sample];
		int alterationIndex = alterations.indices[alterationType];
		if (geneAlterations[geneIndex].find(sampleIndex)
				== geneAlterations[geneIndex].end())
			geneAlterations[geneIndex][sampleIndex] = 0;
		int bitMask = (1 << alterationIndex);
		geneAlterations[geneIndex][sampleIndex] |= bitMask;
	}
	fclose(fin);
	fprintf(stderr, "done. (%.2lf seconds)\n",
			double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf(stderr,
			"\tThere are %lu samples, with a total of %lu genes, harboring %lu different alterations.\n",
			samples.indices.size(), genes.indices.size(),
			alterations.indices.size());
}

bool colorful_test(colored_node_set &candidate, int colorful_option) {
	int color_count[] = { 0, 0, 0, 0, 0 };
	switch (colorful_option) {
	case 1: /* Colorful Subnetworks. */{
		std::set<int> color_set;
		for (auto node : candidate) {
			color_set.insert(node.second);
			if (color_set.size() >= 2)
				return true;
		}
		return false;
	}
	case 2: /* At most 2 non-expression outliers. */{
		for (auto node : candidate)
			color_count[node.second - 1]++;
		for (int i = 0; i < 4; i++)
			if (i != alterations.indices["EXPROUT"])
				color_count[4] += color_count[i];
		if (color_count[4] > 0 && color_count[4] <= 2)
			return true;
		return false;
	}
	case 3: /* Copy num only subnetworks. */{
		for (auto node : candidate)
			color_count[node.second - 1]++;
		for (int i = 0; i < 4; i++)
			if (i != alterations.indices["AMP"])
				color_count[4] += color_count[i];
		if (color_count[4] == 0)
			return true;
		return false;
	}
	case 4: /* Expression outlier only subnetworks. */{
		for (auto node : candidate)
			color_count[node.second - 1]++;
		for (int i = 0; i < 4; i++)
			if (i != alterations.indices["EXPROUT"])
				color_count[4] += color_count[i];
		if (color_count[4] == 0)
			return true;
		return false;
	}
	default:
		break;
	}
	return false;
}

/*
 The main function that preprocesses the data and finds the maximum conserved subnetworks.
 colorful_option = 0: no additional constraints;
 colorful_option = 1: maximum colorful conserved subnetworks;
 colorful_option = 2: maximum colorful conserved subnetworks with at most 2 non-expression outliers;
 colorful_option = 3: maximum copy num only subnetworks;
 colorful_option = 4: maximum expression outlier only subnetworks;
 */
int runMcsiSolver(int S, int t, double errorRate, int colorful_option,
		bool output_option) {
	std::vector<std::unordered_map<colored_node_set, int, colored_node_set_hash>> subNetwork;
	std::vector<std::unordered_map<int, PatientBitmask*>> patientProfile;
	std::unordered_map<colored_node_set, int, colored_node_set_hash> subNetwork_0;
	subNetwork.push_back(subNetwork_0);
	std::unordered_map<int, PatientBitmask*> patientProfile_0;
	patientProfile.push_back(patientProfile_0);

	/* STAGE 1: Initialization of the maximum subnetwork discovery process with the single-node networks of all colored nodes. */
	int color_count[G.V][alterations.indices.size() + 1] = { 0 };
	for (int j = 0; j < G.V; j++)
		for (size_t i = 0; i < samples.indices.size(); i++)
			for (size_t a = 0; a < alterations.indices.size(); a++)
				if ((geneAlterations[j][i] & (1 << a)) != 0)
					color_count[j][a + 1]++;
	int sum = 0;
	for (int j = 0; j < G.V; j++) {
		for (size_t i = 1; i <= alterations.indices.size(); i++)
			if (color_count[j][i] >= t) {
				color_count[j][i] = 1;
				color_count[j][0] = 1;
				if (color_count[j][0] == 1)
					sum++;
			} else
				color_count[j][i] = 0;
	}
	fprintf(stderr,
			"There are %d nodes where at least %d patients are mutated.\n", sum,
			t);

	/* STAGE 2: Initialization of the maximum subnetwork discovery process with the colored edges. */
	int valid = 0;
	for (int j = 0; j < G.V; j++)
		if (color_count[j][0] == 1)
			for (auto index : G.idxInNeighbourList[j])
				if (color_count[index.first][0] == 1 && index.first > j)
					for (size_t k = 1; k <= alterations.indices.size(); k++)
						for (size_t l = 1; l <= alterations.indices.size(); l++)
							if (color_count[j][k] == 1
									&& color_count[index.first][l] == 1) {
								PatientBitmask *profile_e_k_l =
										new PatientBitmask(
												samples.indices.size());
								colored_node_set edge;
								colored_node node_1(j, (int) k);
								edge.insert(node_1);
								colored_node node_2(index.first, (int) l);
								edge.insert(node_2);
								sum = 0;
								for (size_t m = 0; m < samples.indices.size();
										m++)
									if ((geneAlterations[j][m] & (1 << (k - 1)))
											&& (geneAlterations[index.first][m]
													& (1 << (l - 1)))) {
										sum++;
										profile_e_k_l->setBit(m, 1);
									}
								if (sum >= t) {
									auto result = subNetwork[0].insert(
											std::pair<colored_node_set, int>(
													edge, valid));
									if (result.second) {
										patientProfile[0][valid] =
												profile_e_k_l;
										valid++;
									}
									else
										delete profile_e_k_l;
								}
								else
									delete profile_e_k_l;
							}
	fprintf(stderr,
			"There are %lu edges where at least %d patients are mutated at each node.\n",
			subNetwork[0].size(), t);

	/* STAGE 3: Incremental identification of subnetworks with n nodes by extending already identified subnetworks with n - 1 nodes. */
	size_t network_size = 1;
	if (colorful_option == 0) {
		while (true) {
			std::unordered_map<colored_node_set, int, colored_node_set_hash> nextsubNetwork;
			std::unordered_map<int, PatientBitmask*> nextpatientProfile;
			valid = 0;
			for (auto s1 : subNetwork[network_size - 1]) {
				for (auto s2 : subNetwork[0]) {
					colored_node_set candidate = s1.first;
					candidate.insert(s2.first.begin(), s2.first.end());
					if (candidate.size() == network_size + 2) {
						PatientBitmask *profile_candidate = new PatientBitmask(
								patientProfile[network_size - 1][s1.second]);
						profile_candidate->mergeBitmask(
								patientProfile[0][s2.second]);
						if (profile_candidate->getSize() >= t) {
							auto result = nextsubNetwork.insert(
									std::pair<colored_node_set, int>(candidate,
											valid));
							if (result.second) {
								nextpatientProfile[valid] = profile_candidate;
								valid++;
							}
							else
								delete profile_candidate;
						}
						else
							delete profile_candidate;
					}
				}
			}
			/* Found the largest subnetwork. */
			if (nextsubNetwork.empty()) {
				fprintf(stderr, "The maximum subnetwork size is %lu.\n",
						network_size + 1);
				break;
			}
			/* Try next level. */
			fprintf(stderr,
					"There are %lu subnetworks of size %lu, where at least %d patients are mutated at each node.\n",
					nextsubNetwork.size(), network_size + 2, t);
			subNetwork.push_back(nextsubNetwork);
			patientProfile.push_back(nextpatientProfile);
			network_size++;
			/* Network size exceeds limit. */
			if (network_size + 1 >= S)
				break;
		}
	} else {
		while (true) {
			std::unordered_map<colored_node_set, int, colored_node_set_hash> nextsubNetwork;
			std::unordered_map<int, PatientBitmask*> nextpatientProfile;
			valid = 0;
			for (auto s1 : subNetwork[network_size - 1]) {
				for (auto s2 : subNetwork[0]) {
					colored_node_set candidate = s1.first;
					candidate.insert(s2.first.begin(), s2.first.end());
					if (candidate.size() == network_size + 2) {
						PatientBitmask *profile_candidate = new PatientBitmask(
								patientProfile[network_size - 1][s1.second]);
						profile_candidate->mergeBitmask(
								patientProfile[0][s2.second]);
						if (colorful_test(candidate, colorful_option)
								&& (profile_candidate->getSize() >= t)) {
							auto result = nextsubNetwork.insert(
									std::pair<colored_node_set, int>(candidate,
											valid));
							if (result.second) {
								nextpatientProfile[valid] = profile_candidate;
								valid++;
							}
							else
								delete profile_candidate;
						}
						else
							delete profile_candidate;
					}
				}
			}
			/* Found the largest subnetwork. */
			if (nextsubNetwork.empty()) {
				fprintf(stderr,
						"The maximum colorful subnetwork size is %lu.\n",
						network_size + 1);
				break;
			}
			/* Try next level. */
			fprintf(stderr,
					"There are %lu subnetworks of size %lu, where at least %d patients are mutated at each node.\n",
					nextsubNetwork.size(), network_size + 2, t);
			subNetwork.push_back(nextsubNetwork);
			patientProfile.push_back(nextpatientProfile);
			network_size++;
			/* Network size exceeds limit. */
			if (network_size + 1 >= S)
				break;
		}
	}

	/* STAGE 5: Extending the maximum subnetworks to include samples in which there aren't all exact matches, but there aren't colour conflicts either. */
	if (network_size > 5 && errorRate >= 1.0 / (network_size + 1)) {
		llu numSubnetworksExtended = 0;
		llu numSamplesAdded = 0;
		fprintf(stderr, "Extending subnetworks to include %d%% errors...\n",
				int(errorRate * 100));
		for (auto subgraph : subNetwork[subNetwork.size() - 1]) {
			bool extended = false;
			for (size_t i = 0; i < samples.indices.size(); i++) {
				if (!patientProfile[subNetwork.size() - 1][subgraph.second]->getBit(
						i)) {
					int agree = 0;
					for (auto node : subgraph.first)
						if (geneAlterations[node.first][i]
								& (1 << (node.second - 1)))
							agree++;
						else if (geneAlterations[node.first][i]) // Colors conflict
							continue;
					if (agree * 1.0 / (network_size + 1) >= 1 - errorRate) {
						numSamplesAdded++;
						extended = true;
					}
				}
			}
			if (extended)
				numSubnetworksExtended++;
		}
		fprintf(stderr, "%llu subnetworks have been extended.\n",
				numSubnetworksExtended);
		fprintf(stderr, "Average number of samples added is %.1lf\n",
				double(numSamplesAdded) / numSubnetworksExtended);
	}

	/* Output solution(s). Distinct from MCSC, The MCSI solver output file(s) to the current directory. */
	if (output_option) {
		string outFile;
		if (colorful_option == 0 || colorful_option >= 3)
			outFile = "subnetworks_t=" + std::to_string(t) + ".tsv";
		else
			outFile = "colorful_subnetworks_t=" + std::to_string(t) + ".tsv";
		FILE * fout = fopen(outFile.c_str(), "w");
		size_t i = 0;
		fprintf(fout, "Solution\tNodes\tColors\tSampleID\n");
		for (auto subgraph : subNetwork[subNetwork.size() - 1]) {
			fprintf(fout, "Solution_%lu\t", i + 1);
			size_t k = 0;
			for (auto node : subgraph.first) {
				if (k < subgraph.first.size() - 1)
					fprintf(fout, "%s:", G.nodeNames[node.first].c_str());
				else
					fprintf(fout, "%s\t", G.nodeNames[node.first].c_str());
				k++;
			}
			k = 0;
			for (auto node : subgraph.first) {
				if (k < subgraph.first.size() - 1)
					fprintf(fout, "%s:",
							alterations.names[node.second - 1].c_str());
				else
					fprintf(fout, "%s\t",
							alterations.names[node.second - 1].c_str());
				k++;
			}
			k = 0;
			for (size_t j = 0; j < samples.indices.size(); j++)
				if (patientProfile[subNetwork.size() - 1][subgraph.second]->getBit(
						j) == 1) {
					if (k
							< patientProfile[subNetwork.size() - 1][subgraph.second]->getSize()
									- 1)
						fprintf(fout, "%s:", samples.names[j].c_str());
					else
						fprintf(fout, "%s\n", samples.names[j].c_str());
					k++;
				}
			i++;
		}
		fclose(fout);
	}
	/* Deconstruct patient alteration profiles. */
	for (size_t i = 0; i < patientProfile.size(); i++) {
		for (auto profile_map : patientProfile[i]) {
			delete profile_map.second;
			//profile_map.second = NULL;
		}
		patientProfile[i].clear();
	}
	return network_size + 1;
}

/* Reset alteration profiles. */
void clearAlterationProfiles() {
	delete[] samples.names;
	delete[] genes.names;
	//delete []alterations.names;
	samples.indices.clear();
	genes.indices.clear();
	//alterations.indices.clear();
	delete[] geneAlterations;
}

int main(int argc, char * argv[]) {
	printHeader("MSCI SOLVER");
	/* INPUT CHECK */
	if (argc <= 1) {
		fprintf(stderr,
				"./mcsi -p [for p value simulation; optional] -n [network] -l [alteration profiles] -r [color options in subnetwork] -s [maximum subnetwork size] -t [minimum subgraph recurrence] -e [error; optional]\n\n");
		return 0;
	}
	char consoleFlags[] = { 'n', 'l', 's', 't', 'e', 'r', 'p', 0 };
	bool optional[200] = { };
	optional['e'] = true;
	optional['p'] = true;
	unordered_map<char, string> consoleParameters;
	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-' && argv[i][1] && i + 1 < argc
				&& argv[i + 1][0] != '-') {
			consoleParameters[argv[i][1]] = string(argv[i + 1]);
			i++;
		}
		/* Special option p: run */
		if (argv[i][0] == '-' && argv[i][1] == 'p' && i + 1 < argc) {
			consoleParameters[argv[i][1]] = string(0);
			i++;
		}
	}
	for (char * ptrFlag = consoleFlags; *ptrFlag; ptrFlag++) {
		if (!optional[*ptrFlag] && !consoleParameters.count(*ptrFlag)) {
			fprintf(stderr,
					"\n< Error > Missing value for parameter '%c'. Exiting program.\n",
					*ptrFlag);
			exit(0);
		}
	}

	int maxSubnetworkSize;
	int minSubnetworkRecurrence;
	double errorRate = 0;
	int colorfulMode = 0;
	sscanf(consoleParameters['s'].c_str(), "%d", &maxSubnetworkSize);
	sscanf(consoleParameters['t'].c_str(), "%d", &minSubnetworkRecurrence);
	sscanf(consoleParameters['r'].c_str(), "%d", &colorfulMode);
	if (consoleParameters.count('e')) {
		sscanf(consoleParameters['e'].c_str(), "%lf", &errorRate);
		fprintf(stderr, "Error rate set to %.2lf\n", errorRate);
	}
	printHeader("Reading Input");
	readUndirectedNetwork(consoleParameters['n'].c_str());
	findConnectedComponents();
	if (!consoleParameters.count('p')) {
		readAlterationProfiles(consoleParameters['l'].c_str(), 0);
		printHeader("Solving the problem");
		runMcsiSolver(maxSubnetworkSize, minSubnetworkRecurrence, errorRate,
				colorfulMode, true);
	} else {
		for (int depth = minSubnetworkRecurrence;
				depth <= minSubnetworkRecurrence + 9; depth++) {
			std::vector<int> network_sizes;
			for (int i = 0; i < 1000; i++) {
				size_t file_separator = consoleParameters['l'].find('.');
				string input_profiles;
				if (file_separator != string::npos)
					input_profiles = consoleParameters['l'].substr(0,
							file_separator) + std::to_string(i)
							+ consoleParameters['l'].substr(file_separator);
				else
					exit(0);
				readAlterationProfiles(input_profiles.c_str(), i);
				int res = runMcsiSolver(maxSubnetworkSize, depth, errorRate,
						colorfulMode, false);
				network_sizes.push_back(res);
				clearAlterationProfiles();
			}
			string outFile;
			switch (colorfulMode) {
			case 0:
				outFile = "subnetwork_sizes_t=" + std::to_string(depth)
						+ ".txt";
				break;
			case 1:
				outFile = "colorful_subnetwork_sizes_t=" + std::to_string(depth)
						+ ".txt";
				break;
			default:
				outFile = "special_subnetwork_sizes_t=" + std::to_string(depth)
						+ ".txt";
				break;
			}
			FILE *fout = fopen(outFile.c_str(), "w");
			for (int i = 0; i < 1000; i++)
				fprintf(fout, "%d\n", network_sizes[i]);
			fclose(fout);
		}
	}
	for (int i = 0; i < G.V; i++) {
		delete[] G.N[i];
	}
	delete[] G.N;
	delete[] G.NSize;
	//G.NSize = NULL;
	delete[] G.nodeNames;
	delete[] G.ccIndex;
	delete[] G.ccSize;
	delete[] G.idxInNeighbourList;
	delete[] G.incomingEdges;
	delete[] samples.names;
	delete[] genes.names;
	delete[] alterations.names;
	delete[] geneAlterations;
	return 0;
}

