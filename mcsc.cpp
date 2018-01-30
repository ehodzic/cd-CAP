/*
	KNOWN BUG: none so far
*/

#include <cstdio>
#include <cstring>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <ilcplex/ilocplex.h>
using namespace std;

typedef unsigned long long llu;
typedef long long ll;

// Fancy Output Functions
void printHeader(const char * text) {
	int n = strlen(text);
	fprintf( stderr, "\n" );
	for (int i = 0; i < n + 4; i++) fprintf(stderr, "*");
	fprintf(stderr, "\n* %s *\n", text);
	for (int i = 0; i < n + 4; i++) fprintf(stderr, "*");
	fprintf(stderr, "\n");
}

// Global variables
struct stringPairHash {
	size_t operator() (const pair<string, string> & Q) const {
		size_t h1 = hash<string>() (Q.first);
		size_t h2 = hash<string>() (Q.second);
		return h1 ^ (h2 << 1);
	}
};

struct Entry {
	string * names;
	unordered_map<string, int> indices;
} samples, genes, alterations;

struct Bitmask {
	int maxSize;
	int size;
	int len;
	llu * bits;

	Bitmask(int maxSize) : maxSize(maxSize) {
		size = 0;
		len = ceil(maxSize / 64.0);
		bits = new llu [len];
		memset(bits, 0, sizeof(llu)*len);
	}

	Bitmask(const Bitmask & Q) {
		maxSize = Q.maxSize;
		size = Q.size;
		len = Q.len;
		bits = new llu [len];
		for (int i = 0; i < len; i++)
			bits[i] = Q.bits[i];
	}

	Bitmask(Bitmask * Q) {
		maxSize = Q -> maxSize;
		size = Q -> size;
		len = Q -> len;
		bits = new llu [len];
		for (int i = 0; i < len; i++)
			bits[i] = Q -> bits[i];
	}

	~Bitmask() {
		if (bits)
			delete bits;
	}

	bool operator== (const Bitmask & Q) const {
		if (size != Q.size) return false;
		int minLen = (len < Q.len) ? len : Q.len;
		for (int i = 0; i < minLen; i++) if (bits[i] != Q.bits[i]) return false;
		for (int i = minLen; i < len; i++) if (bits[i]) return false;
		for (int i = minLen; i < Q.len; i++) if (Q.bits[i]) return false;
		return true;
	}

	void copylluBitmask(llu x) {
		len = 1;
		bits[0] = x;
		size = __builtin_popcountll(x);
	}

	void setBit(int pos, bool val) {
		int idx = pos / 64;
		if (idx >= len) {
			fprintf(stderr, "< Error > Cannot assign bit because bitmask is too small. pos: %d | len: %d | idx: %d | maxSize: %d\n", pos, len, idx, maxSize);
			exit(0);
		}
		int bitIdx = pos % 64;
		bool oldVal = getBit(pos);
		if (oldVal ^ val) {	// the bit is about to get changed
			bits[idx] ^= llu(1) << bitIdx;
			if (val) size++;
			else size--;
		}
	}

	bool getBit(int pos) const {
		int idx = pos / 64;
		if (idx >= len) {
			fprintf(stderr, "< Error > Cannot return bit because bitmask is too small.\n");
			exit(0);
		}
		int bitIdx = pos % 64;
		return bits[idx] & (llu(1) << bitIdx);
	}

	int getSize() const { return size; }

	void invert() {
		for (int i = 0; i < len; i++)
			bits[i] = ~bits[i];
	}

	int getPositionOfFirstSetBit() const {
		if (size == 0) {
			fprintf(stderr, "< Error > Cannot return position of first set bit because bitmask is empty.\n");
			exit(0);
		}
		for (int i = 0; i < len; i++) {
			if (bits[i]) return i*64 + __builtin_ctzll(bits[i]);
		}
	}

	int extractLowestOrderSetBitIndex() {
		if (size == 0) {
			fprintf(stderr, "< Error > Cannot extract first set bit because bitmask is empty.\n");
			exit(0);
		}
		for (int i = 0; i < len; i++) {
			if (bits[i]) {
				int pos = i*64 + __builtin_ctzll(bits[i]);
				setBit(pos, 0);
				return pos;
			}
		}
	}
};

struct BitmaskHasher {
	std::size_t operator()(const Bitmask & Q) const {
		using std::size_t;
		using std::hash;
		size_t res = 17;
		for (int i = 0; i < Q.len; i++) {
			res = res * 31 + hash< unsigned long long >()(Q.bits[i]);
		}
		return res;
	}
};

struct SubnetworkEntry {
	vector<int> nodes;
	vector<int> nodeColourIdx;
	Bitmask * samples;
	int seedSampleIdx;
	bool isValid;

	SubnetworkEntry() : samples(0), isValid(true) {}

	SubnetworkEntry(int totalNumSamples) : isValid(true) {
		samples = new Bitmask(totalNumSamples);
	}

	SubnetworkEntry(vector<int> & nodeV, int totalNumSamples) : isValid(true) {
		nodes = nodeV;
		samples = new Bitmask(totalNumSamples);
	}

	SubnetworkEntry(vector<int> & nodeV, Bitmask & sampleB) : isValid(true) {
		nodes = nodeV;
		samples = new Bitmask(sampleB);
	}

	SubnetworkEntry(const SubnetworkEntry & Q) : isValid(Q.isValid), seedSampleIdx(Q.seedSampleIdx) {
		nodes = Q.nodes;
		nodeColourIdx = Q.nodeColourIdx;
		samples = new Bitmask(Q.samples);
	}

	~SubnetworkEntry(){
		if (samples) {
			delete samples;
			samples = 0;
		}
	}

	void print(unordered_map<int, llu> * geneAlterations, Entry & sampleInfo, Entry & alterationInfo, int * chrArm, string * nodeNames, FILE * fout) {
		fprintf(fout, "Patients\t%d\n", this->numSamples());
		Bitmask tempmask(this->samples);
		while (tempmask.getSize()) {
			int sampleIdx = tempmask.extractLowestOrderSetBitIndex();
			int numColoured = 0;
			for (int nodeIdx : this->nodes) {
				if (geneAlterations[nodeIdx].count(sampleIdx))
					numColoured++;
			}
			fprintf(fout, " %s(%d)", sampleInfo.names[sampleIdx].c_str(), numColoured);
		}
		fprintf(fout, "\nGenes\t%d\n", this->nodes.size());
		vector<llu> nodeColourMasks = this->getNodeColourBitmaskVector(geneAlterations);
		for (int i = 0; i < this->nodes.size(); i++) {
			int nodeIdx = this->nodes[i];
			fprintf(fout, "%s\t", nodeNames[nodeIdx].c_str());
			Bitmask tempmask(64);
			tempmask.copylluBitmask(nodeColourMasks[i]);
			while (tempmask.getSize()) {
				int colourIdx = tempmask.extractLowestOrderSetBitIndex();
				fprintf(fout, "\t%s", alterationInfo.names[colourIdx].c_str());
			}
			if (chrArm[nodeIdx]) {
				int chrNum = chrArm[nodeIdx] / 2;
				bool isQ = chrArm[nodeIdx] % 2;
				if (chrNum < 23)
					fprintf(fout, "\tchr%d", chrNum);
				else
					fprintf(fout, "\tchr%c", chrNum == 23 ? 'X' : 'Y');
				fprintf(fout, "%c", isQ ? 'q' : 'p');
			}
			else
				fprintf(fout, "\t-");
			fprintf(fout, "\n");
		}
	}

	int numSamples() const { return this->samples->getSize(); }

	// Returns true if the bitmask changes (shrinks), false if nothing gets changed
	bool fixSamplesViaNode(int nodeIdx, llu nodeColourBitmask, unordered_map<int, llu> * geneAlterations) {
		Bitmask tempmask(this->samples);
		bool aSampleDiscarded = false;
		while (tempmask.getSize()) {
			int sampleIdx = tempmask.extractLowestOrderSetBitIndex();
			// If the node is not coloured in this sample, or it is but the colours have no intersection
			if ( !geneAlterations[nodeIdx].count(sampleIdx) || (geneAlterations[nodeIdx][sampleIdx] & nodeColourBitmask) == 0 ) {
				this->samples->setBit(sampleIdx, 0);	// Discard this sample
				aSampleDiscarded = true;
			}
		}
		return aSampleDiscarded;
	}

	void buildSamplesViaNode(int nodeIdx, llu nodeColourBitmask, unordered_map<int, llu> * geneAlterations, int numSamples) {
		for (int i = 0; i < numSamples; i++) {
			// If the node is coloured in this sample and the colours have an intersection
			if (geneAlterations[nodeIdx].count(i) && (geneAlterations[nodeIdx][i] & nodeColourBitmask))
				this->samples->setBit(i, 1);
		}
	}

	llu getNodeColourBitmask(int nodeIdx, unordered_map<int, llu> * geneAlterations) const {
		llu colourMask = 0;
		colourMask--;
		Bitmask tempmask(this->samples);
		while (tempmask.getSize()) {
			int sampleIdx = tempmask.extractLowestOrderSetBitIndex();
			if (geneAlterations[nodeIdx].count(sampleIdx))	// Required if the subnetwork got extended to this sample with error
				colourMask &= geneAlterations[nodeIdx][sampleIdx];
		}
		return colourMask;
	}

	vector<llu> getNodeColourBitmaskVector(unordered_map<int, llu> * geneAlterations) const {
		vector<llu> nodeColourMasks(this->nodes.size());
		for (int j = 0; j < this->nodes.size(); j++)
			nodeColourMasks[j] = getNodeColourBitmask(this->nodes[j], geneAlterations);
		return nodeColourMasks;
	}

	// Checks whether the subnetwork can be extended to the given sample with the given error rate; conditioned upon no node having conflicting colouring, and only colourless nodes being allowed.
	bool supportsSampleWithError(int sampleIdx, unordered_map<int, llu> * geneAlterations, double errorRate) const {
		int agree = 0;
		for (int j = 0; j < this->nodes.size(); j++) {
			int const & nodeIdx = this->nodes[j];
			llu const & nodeColourMask = llu(1) << this->nodeColourIdx[j];
			if (geneAlterations[nodeIdx].count(sampleIdx)) {
				if (geneAlterations[nodeIdx][sampleIdx] & nodeColourMask)
					agree++;
				else 	// Colours conflict
					return false;
			}
		}
		int colourless = this->nodes.size() - agree;
		return (double(colourless)/this->nodes.size() <= errorRate);
	}

	void extendSubnetworkWithError(unordered_map<int, llu> * geneAlterations, Entry & sampleInfo, double errorRate) {
		// vector<llu> nodeColourMasks = this->getNodeColourBitmaskVector(geneAlterations);
		for (int i = 0; i < sampleInfo.indices.size(); i++) {
			if (!(this->samples -> getBit(i)) && this->supportsSampleWithError(i, geneAlterations, errorRate)) {
				this->samples -> setBit(i, 1);
			}
		}
	}
};

struct Graph {
	int V, E;
	int * NSize;
	int ** N;
	string * nodeNames;
	int * chrArm;
	unordered_map<string, int> nodeIndices;
} G;

struct Subgraph {
	int V;
	int * nodeNames;
	unordered_map<int, int> nameIdx;
	int * degrees;
	int ** edges;
	pair<int, int> ** incomingEdges;

	Subgraph(){}
	~Subgraph(){
		delete nodeNames;
		delete degrees;
		for (int i = 0; i < V; i++) {
			delete [] edges[i];
			delete [] incomingEdges[i];
		}
		delete edges;
		delete incomingEdges;
	}
	Subgraph( Graph * G, int * nodes, int numNodes ) : V(numNodes) {
		nodeNames = new int [V];
		degrees = new int [V];
		edges = new int * [V];
		incomingEdges = new pair<int, int> * [V];
		bool * flagged = new bool[G -> V];
		memset(degrees, 0, V * sizeof(degrees[0]));
		memset(flagged, 0, (G -> V) * sizeof(flagged[0]));
		for (int i = 0; i < V; i++) {
			flagged[nodes[i]] = true;
			nodeNames[i] = nodes[i];
			nameIdx[ nodeNames[i] ] = i;
		}
		/******************************************************
		 * Flagged all nodes that belong to current subgraph. *
		 ******************************************************/
		for (int i = 0; i < V; i++) {
			int node = nodes[i];
			int NSize = G -> NSize[node];
			for (int j = 0; j < NSize; j++) {
				int neighbour = G -> N[node][j];
				if (flagged[neighbour]) {
					degrees[i]++;
				}
			}
			edges[i] = new int [ degrees[i] + 1 ];
			incomingEdges[i] = new pair<int, int> [ degrees[i] + 1 ];
			edges[i][ degrees[i] ] = 0;
			incomingEdges[i][ degrees[i] ] = make_pair(0, 0);
		}
		/*********************************************
		 * Extracted subgraph degrees for all nodes. *
		 *********************************************/
		for (int i = 0; i < V; i++) {
			int node = nodes[i];
			int NSize = G -> NSize[node];
			for (int j = 0; j < NSize; j++) {
				int neighbour = G -> N[node][j];
				if (flagged[neighbour]) {
					// node -> neighbour directed edge discovered
					// int nodeIdx = nameIdx[node];
					int neighbourIdx = nameIdx[neighbour];
					int & edgesIndex = edges[i][ degrees[i] ];
					edges[i][ edgesIndex ] = neighbourIdx;
					int incomingEdgesIndex = incomingEdges[neighbourIdx][ degrees[neighbourIdx] ].first;
					incomingEdges[neighbourIdx][ incomingEdgesIndex ] = make_pair(i, edgesIndex);
					edgesIndex++;
					incomingEdges[neighbourIdx][ degrees[neighbourIdx] ] = make_pair(incomingEdgesIndex + 1, 0);
				}
			}
		}
		/**************************************
		 * Build edge and incoming edge lists *
		 **************************************/
		delete flagged;
	}
};

unordered_map<int, llu> * geneAlterations;	// geneAlterations[ i ][ j ] = c means that "gene i has colour c in patient j". The colours are bitmasks (so supporting max 64 different alteration types).
vector< pair<int,int> > subnetworkSeeds;

/*
	Reads the input -n parameter as a collection of undirected edges (pairs of node names, separated by whitespace) and stores the information in global Graph object G.
*/
void readUndirectedNetwork(const char * filename) {
	fprintf(stderr, "Reading the network... ");
	int timerStart = clock();
	unordered_set<pair<string, string>, stringPairHash> uniqueEdges;
	char u[1000], v[1000];
	FILE * fin = NULL;
	if (!(fin = fopen(filename, "r"))) {
		fprintf(stderr, "\n< Error > Cannot open file '%s'. Please make sure the file exists.\n", filename);
		exit(0);
	}
	G.V = 0;
	while (fscanf(fin, "%s%s", u, v) == 2) {
		char * a = u, * b = v;
		if (strcmp(a, b) > 0) swap(a, b);
		else if (strcmp(a, b) == 0) continue;
		pair<string, string> e = make_pair(string(a), string(b));
		if (!uniqueEdges.count(e)) {	// Previously NOT seen undirected edge
			uniqueEdges.insert(e);
			if (!G.nodeIndices.count(e.first)) {
				G.nodeIndices[e.first] = G.V;
				G.V++;
			}
			if (!G.nodeIndices.count(e.second)) {
				G.nodeIndices[e.second] = G.V;
				G.V++;
			}
		}
	}
	fclose(fin);
	G.E = uniqueEdges.size() * 2;

	G.NSize = new int [G.V];
	G.N = new int * [G.V];
	G.nodeNames = new string[G.V];
	G.chrArm = new int [G.V];
	memset(G.NSize, 0, sizeof(G.NSize[0]) * G.V);
	memset(G.chrArm, 0, sizeof(G.chrArm[0]) * G.V);

	for (auto e: uniqueEdges) {
		int idx1 = G.nodeIndices[e.first];
		if (G.nodeNames[idx1].size() == 0) G.nodeNames[idx1] = e.first;

		int idx2 = G.nodeIndices[e.second];
		if (G.nodeNames[idx2].size() == 0) G.nodeNames[idx2] = e.second;

		G.NSize[idx1]++;
		G.NSize[idx2]++;
	}

	for (int i = 0; i < G.V; i++) {
		G.N[i] = new int[ G.NSize[i] ];
	}

	int * tempNSize = new int [G.V];
	memset(tempNSize, 0, sizeof(tempNSize[0]) * G.V );

	for (auto e: uniqueEdges) {
		int idx1 = G.nodeIndices[e.first];
		int idx2 = G.nodeIndices[e.second];
		G.N[ idx1 ][ tempNSize[idx1] ] = idx2;
		G.N[ idx2 ][ tempNSize[idx2] ] = idx1;
		tempNSize[idx1]++;
		tempNSize[idx2]++;
	}
	delete [] tempNSize;
	fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf(stderr, "\tInput network contains %d nodes and %d undirected edges.\n", G.V, G.E / 2);
}

/*
	Reads the chromosome information for genes from the -c parameter
*/
void readChromosomeInfo(const char * filename) {
	fprintf(stderr, "Reading the chromosome information... ");
	int timerStart = clock();
	char gene[1000], chromosome[1000], karyotypeBand[1000];
	FILE * fin = fopen(filename, "r");
	while (fscanf(fin, "%s%s%s", gene, chromosome, karyotypeBand) == 3) {
		if ( !G.nodeIndices.count(string(gene)) ) continue;
		if (chromosome[0]=='X') strcpy(chromosome, "23");
		else if (chromosome[0]=='Y') strcpy(chromosome, "24");
		int chrIdx;
		sscanf(chromosome, "%d", &chrIdx);
		int chrArm = (karyotypeBand[0] == 'p') ? 0 : 1;
		int geneIndex = G.nodeIndices[string(gene)];
		G.chrArm[geneIndex] = chrIdx * 2 + chrArm;
	}
	fclose(fin);
	fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
}

/*
	Reads information about which genes should be excluded from the set of coloured nodes (for the purpose of results analysis) from the -x parameter
*/
void readExcludeInfo(const char * filename) {
	fprintf(stderr, "Reading the excluded genes... ");
	int timerStart = clock();
	char gene[1000];
	FILE * fin = fopen(filename, "r");
	unordered_set<int> excluded;
	while (fscanf(fin, "%s", gene) == 1) {
		if (!G.nodeIndices.count(string(gene))) continue;
		int geneIndex = G.nodeIndices[string(gene)];
		excluded.insert(geneIndex);
		geneAlterations[geneIndex].clear();	// We are removing the colour of this node.
	}
	fclose(fin);
	int newSize = 0;
	for (llu i = 0; i < subnetworkSeeds.size(); i++) {
		int geneIndex = subnetworkSeeds[i].second;
		if (!excluded.count(geneIndex)) {
			subnetworkSeeds[newSize] = subnetworkSeeds[i];
			newSize++;
		}
	}
	fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf(stderr, "%d nodes and %d possible subnetwork seeds excluded.\n", excluded.size(), subnetworkSeeds.size() - newSize);
	subnetworkSeeds.resize(newSize);
}

/*
	Reads the input -l parameter as a collection of "sample gene alterationType" triples, separated by whitespace.
*/
void readAlterationProfiles(const char * filename) {
	fprintf(stderr, "Reading the alteration profiles... ");
	int timerStart = clock();
	char sample[1000], gene[1000], alterationType[1000];
	FILE * fin = NULL;
	if (!(fin = fopen(filename, "r"))) {
		fprintf(stderr, "\n< Error > Cannot open file '%s'. Please make sure the file exists.\n", filename);
		exit(0);
	}
	while (fscanf(fin, "%s%s%s", sample, gene, alterationType) == 3) {
		if (!G.nodeIndices.count(string(gene))) continue;
		if (!samples.indices.count(string(sample))) {
			int idx = samples.indices.size();
			samples.indices[string(sample)] = idx;
		}
		if (!genes.indices.count(string(gene))) {
			int idx = genes.indices.size();
			genes.indices[string(gene)] = idx;
		}
		if (!alterations.indices.count(string(alterationType))) {
			int idx = alterations.indices.size();
			alterations.indices[string(alterationType)] = idx;
		}
	}
	samples.names = new string[samples.indices.size()];
	genes.names = new string[genes.indices.size()];
	alterations.names = new string[ alterations.indices.size() ];
	for (auto it : samples.indices) {
		samples.names[it.second] = it.first;
	}
	for (auto it : genes.indices) {
		genes.names[it.second] = it.first;
	}
	for (auto it : alterations.indices) {
		alterations.names[it.second] = it.first;
	}
	rewind(fin);
	geneAlterations = new unordered_map<int, llu> [ G.V ];
	while (fscanf(fin, "%s%s%s", sample, gene, alterationType) == 3) {
		if (!genes.indices.count(string(gene))) continue;
		int geneIndex = G.nodeIndices[string(gene)];
		int sampleIndex = samples.indices[sample];
		// int subnetworkIndex = G.nodeIndices[string(gene)] + sampleIndex * G.V;
		subnetworkSeeds.push_back(make_pair(sampleIndex, geneIndex));
		int alterationIndex = alterations.indices[alterationType];
		if ( !geneAlterations[geneIndex].count(sampleIndex) )
			geneAlterations[geneIndex][sampleIndex] = 0;
		llu bitMask = (llu(1) << alterationIndex);
		bitMask |= geneAlterations[geneIndex][sampleIndex];
		geneAlterations[geneIndex][sampleIndex] = bitMask;
	}
	fclose(fin);
	fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf( stderr, "\tThere are %lu samples, with a total of %lu genes, harboring %lu different alterations.\n", samples.indices.size(), genes.indices.size(), alterations.indices.size() );
	fprintf(stderr, "\tThere are %lu possible subnetwork seeds.\n", subnetworkSeeds.size());
}

/*
	The main function that preprocesses the data and runs the CPLEX ILP solver.
*/
void runSolver(int S, int t, int K, double errorRate, const char * folderName, int threads, int seconds, int minColours) {
	/*
		STAGE 1: Identification of connected components among the coloured nodes in the PPI network.
		Purpose: Minimization of the flow network size for each possible seed.
	*/
	fprintf(stderr, "Finding coloured patient-specific connected components... ");
	int timerStart = clock();
	int ** node_CCIndex = new int * [samples.indices.size()];	// node_CCIndex[sampleIdx][nodeIdx] = CCIndex
	int * CC_count = new int [samples.indices.size()];	// CC_count[sampleIdx] = numberOfConnectedColouredComponents
	int * nodeStack = new int [G.V + 1];
	int * temp = new int [G.V + 1];
	Subgraph *** CC = new Subgraph ** [samples.indices.size()];
	for (int sampleIdx = 0; sampleIdx < samples.indices.size(); sampleIdx++) {
		node_CCIndex[sampleIdx] = new int [G.V];
		memset(node_CCIndex[sampleIdx], -1, sizeof(node_CCIndex[sampleIdx][0])*G.V);
		CC_count[sampleIdx] = 0;	// Number of connected components
		int & num_of_CCs = CC_count[sampleIdx];
		int nodeStackSize = 0;
		for (int j = 0; j < G.V; j++) {
			if (geneAlterations[j].count(sampleIdx) && node_CCIndex[sampleIdx][j] == -1) {	// node is coloured and not assigned to a connected component
				node_CCIndex[sampleIdx][j] = num_of_CCs++;
				nodeStack[nodeStackSize++] = j;
				while (nodeStackSize) {
					int node = nodeStack[--nodeStackSize];
					for (int j1 = 0; j1 < G.NSize[node]; j1++) {
						int const & neighbour = G.N[node][j1];
						if (geneAlterations[neighbour].count(sampleIdx) && node_CCIndex[sampleIdx][neighbour] == -1) {	// neighbour is coloured and not assigned to a connected component
							node_CCIndex[sampleIdx][neighbour] = node_CCIndex[sampleIdx][node];
							nodeStack[nodeStackSize++] = neighbour;
						}
					}
				}
			}
		}
		/***************************************************************************************************************************************
		 * Computed number of connected coloured components in the sample. (CC_count[sampleIdx])											   *
		 * For every node, assigned index of the connected coloured component it belongs to in the sample. (node_CCIndex[sampleIdx][nodeIdx])  *
		 ***************************************************************************************************************************************/
		int * CCSizes = temp;
		memset(CCSizes, 0, num_of_CCs * sizeof(CCSizes[0]));
		for (int j = 0; j < G.V; j++) if (node_CCIndex[sampleIdx][j] != -1) CCSizes[ node_CCIndex[sampleIdx][j] ]++;
		for (int i = 1; i < num_of_CCs; i++) CCSizes[i] += CCSizes[i - 1];
		int totalColouredNodes = CCSizes[num_of_CCs - 1];
		for (int j = 0; j < G.V; j++) if (node_CCIndex[sampleIdx][j] != -1) nodeStack[ --CCSizes[ node_CCIndex[sampleIdx][j] ] ] = j;
		/***************************************************************************************
		 * Partitioned nodes of the sample network based on coloured connected component index *
		 ***************************************************************************************/
		CC[sampleIdx] = new Subgraph * [num_of_CCs];
		for (int i = 0, CCIndex = 0; i < totalColouredNodes; i++, CCIndex++) {
			int j = i + 1;
			while (j < totalColouredNodes && node_CCIndex[sampleIdx][nodeStack[j]] == node_CCIndex[sampleIdx][nodeStack[i]]) j++;
			CC[sampleIdx][CCIndex] = new Subgraph(&G, nodeStack + i, j - i);
			i = j - 1;
		}
		/****************************************************************
		 * Constructed subgraphs based on connected coloured components *
		 ****************************************************************/
	}
	delete [] nodeStack;
	delete [] temp;
	fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);

	// // 
	// // STAGE 1.5: Checking the number of pairs of differently coloured nodes that occur in at least 't' patients
	// // 
	// fprintf(stderr, "Checking the number of pairs of differently coloured neighbour nodes that occur in at least %d patients.\n", t); timerStart = clock();
	// llu colourfulPairs = 0;
	// unordered_set< llu > matchedPairs;
	// for (int k = 0; k < subnetworkSeeds.size(); k++) {
	// 	fprintf(stderr, "\r%.1lf%%", 100*double(k + 1) / double(subnetworkSeeds.size()));
	// 	auto & seedInfo 	= subnetworkSeeds[k];
	// 	int sampleIdx		= seedInfo.first;
	// 	int nodeIdx			= seedInfo.second;
	// 	int CCIndex			= node_CCIndex[sampleIdx][nodeIdx];
	// 	int nodeInternalIdx	= CC[sampleIdx][CCIndex] -> nameIdx[nodeIdx];
	// 	int NSize			= CC[sampleIdx][CCIndex] -> degrees[nodeInternalIdx];
	// 	for (int nIdx = 0; nIdx < NSize; nIdx++) {	// Going through all the neighbours of the current node of the current subnetwork
	// 		int neighbourInternalIdx	= CC[sampleIdx][CCIndex] -> edges[nodeInternalIdx][nIdx];
	// 		int neighbourIdx			= CC[sampleIdx][CCIndex] -> nodeNames[neighbourInternalIdx];
	// 		llu a = nodeIdx;
	// 		llu b = neighbourIdx;
	// 		llu r = (a < b) ? a * G.V + b : b * G.V + a;
	// 		if (matchedPairs.count(r))
	// 			continue;
	// 		Bitmask nodeColours(64);
	// 		Bitmask neighbourColours(64);
	// 		nodeColours.copylluBitmask(geneAlterations[nodeIdx][sampleIdx]);
	// 		neighbourColours.copylluBitmask(geneAlterations[neighbourIdx][sampleIdx]);
	// 		while (nodeColours.getSize()) {
	// 			int nodeColourIndex = nodeColours.extractLowestOrderSetBitIndex();
	// 			llu nodeColourMask = llu(1) << nodeColourIndex;
	// 			Bitmask tempmask(neighbourColours);
	// 			while (tempmask.getSize()) {
	// 				int neighbourColourIndex = tempmask.extractLowestOrderSetBitIndex();
	// 				llu neighbourColourMask = llu(1) << neighbourColourIndex;
	// 				// We isolated a pair of colours
	// 				if (nodeColourMask ^ neighbourColourMask) {	// We have two different colours
	// 					llu sampleCount = 0;
	// 					for (int i = 0; i < samples.indices.size(); i++) {
	// 						if (geneAlterations[nodeIdx].count(i) && (geneAlterations[nodeIdx][i] & nodeColourMask) && geneAlterations[neighbourIdx].count(i) && (geneAlterations[neighbourIdx][i] & neighbourColourMask))	// Identical colour match
	// 							sampleCount++;
	// 					}
	// 					if (sampleCount >= t) 
	// 						matchedPairs.insert(r);
	// 				}
	// 			}
	// 		}
	// 	}
	// }
	// fprintf(stderr, "\r%llu colourful neighbour pairs found. (%.2lf seconds)\n", (llu) matchedPairs.size(), double(clock() - timerStart) / CLOCKS_PER_SEC);

	// Allocation
	const int numPatients = samples.indices.size();
	vector<SubnetworkEntry> * candidateSubnetworks = new vector<SubnetworkEntry> [S];

	//
	//	STAGE 2: Initialization of the candidate subnetwork discovery process with the single-node networks of all coloured nodes.
	//
	fprintf(stderr, "Constructing initial coloured single-node subnetworks...\n"); timerStart = clock();
	for (int k = 0, lastProg = 0; k < subnetworkSeeds.size(); k++) {
		int progress = 1000 * double(k + 1) / double(subnetworkSeeds.size());
		if (progress > lastProg) {
			fprintf(stderr, "\r%.1lf%%", double(progress)/10.0);
			lastProg = progress;
		}
		auto & seedInfo = subnetworkSeeds[k];
		int sampleIdx	= seedInfo.first;
		int nodeIdx		= seedInfo.second;
		Bitmask nodeColourMask(64);
		nodeColourMask.copylluBitmask(geneAlterations[nodeIdx][sampleIdx]);
		while (nodeColourMask.getSize()) {
			int colourIndex = nodeColourMask.extractLowestOrderSetBitIndex();
			// if (alterations.names[colourIndex] == "EXPROUT")
			// if (alterations.names[colourIndex] != "EXPROUT")	// Using only expression outlier seeds
			// if (alterations.names[colourIndex] == "EXPROUT" || alterations.names[colourIndex] == "AMP")
				// continue;	// Not using expression-outlier seeds.
			SubnetworkEntry newEntry(numPatients);
			newEntry.seedSampleIdx = sampleIdx;
			newEntry.isValid = true;
			newEntry.nodes.push_back(nodeIdx);
			newEntry.nodeColourIdx.push_back(colourIndex);
			llu requiredColourBitmask = llu(1) << colourIndex;
			newEntry.buildSamplesViaNode(nodeIdx, requiredColourBitmask, geneAlterations, samples.indices.size());
			if (newEntry.numSamples() >= t)
				candidateSubnetworks[0].push_back(newEntry);
		}
	}
	fprintf(stderr, "\rDone. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);

	//
	//	STAGE 3: Incremental identification of valid candidate subnetworks with n nodes by extending already identified subnetworks with n-1 nodes.
	//
	llu * visitedNode = new llu[G.V + 1];
	llu * nodesHash = new llu [G.V * S + 1];
	unordered_set<Bitmask, BitmaskHasher> * nodesHash_buckets = new unordered_set<Bitmask, BitmaskHasher> [G.V * S + 1];
	memset(visitedNode, 0, sizeof(visitedNode[0])*(G.V + 1));
	memset(nodesHash, 0, sizeof(nodesHash[0])*(G.V * S + 1));
	llu visitedIdx			= 0;
	llu hashIdx				= 0;
	llu totalNumSubgraphs	= candidateSubnetworks[0].size();
	llu numContained		= 0;
	timerStart = clock();
	for (int cycle = 1; cycle < S; cycle++) {
		fprintf(stderr, "\nCurrent number of subgraphs is %llu. Constructing all candidate subnetworks of size %d...\n", totalNumSubgraphs, cycle + 1);
		hashIdx++;	// We want unique networks for each network size
		for (int i = 0, lastProg = 0; i < candidateSubnetworks[cycle - 1].size(); i++) {
			int progress = 1000 * double(i + 1) / double(candidateSubnetworks[cycle - 1].size());
			if (progress > lastProg) {
				fprintf(stderr, "\r%.1lf%%", double(progress)/10.0);
				lastProg = progress;
			}
			SubnetworkEntry & subnetInfo = candidateSubnetworks[cycle - 1][i];
			visitedIdx++;	// New visited flag for each new subnetwork that is attempted to be extended
			llu subnetHash = 0;
			// Marking nodes visited and calculating the hash value of the subnetwork nodes
			for (int nodeIdx : subnetInfo.nodes) {
				visitedNode[nodeIdx] = visitedIdx;
				subnetHash += nodeIdx;
			}
			int & sampleIdx = subnetInfo.seedSampleIdx;
			// Exploring neighbours and constructing new subnetworks of size greater by 1
			for (int nodeIdx : subnetInfo.nodes) {	// We test neighbours of every node in the current subnetwork that we are seeking to extend
				int CCIndex			= node_CCIndex[sampleIdx][nodeIdx];
				int nodeInternalIdx	= CC[sampleIdx][CCIndex] -> nameIdx[nodeIdx];
				int NSize			= CC[sampleIdx][CCIndex] -> degrees[nodeInternalIdx];
				for (int nIdx = 0; nIdx < NSize; nIdx++) {	// Going through all the neighbours of the current node of the current subnetwork
					int neighbourInternalIdx	= CC[sampleIdx][CCIndex] -> edges[nodeInternalIdx][nIdx];
					int neighbourIdx			= CC[sampleIdx][CCIndex] -> nodeNames[neighbourInternalIdx];
					if (geneAlterations[neighbourIdx].count(sampleIdx) && visitedNode[neighbourIdx] < visitedIdx) {	// The node is actually coloured and we haven't tried it yet with the current subnetwork
						visitedNode[neighbourIdx] = visitedIdx;
						llu newHash = subnetHash + neighbourIdx;	// Sum of a combination of unique 'cycle + 1' numbers has to be unique itself.
						if (nodesHash[newHash] < hashIdx) {	// First time visiting this bucket. Clear it.
							nodesHash[newHash] = hashIdx;
							nodesHash_buckets[newHash].clear();
						}
						// We may have same subnetwork but with different patients involved due to different node colouring. That forces us to continue and compare against patient bitmasks for this hash bucket.
						// Calculate the patient bitmask with this node added
						Bitmask alterationBitmask(64);
						alterationBitmask.copylluBitmask(geneAlterations[neighbourIdx][sampleIdx]);
						while (alterationBitmask.getSize()) {	// Go through its colours.
							int alterationIndex = alterationBitmask.extractLowestOrderSetBitIndex();
							llu singleColourBitmask = llu(1) << alterationIndex;
							SubnetworkEntry newEntry(subnetInfo);
							bool lostSamples = newEntry.fixSamplesViaNode(neighbourIdx, singleColourBitmask, geneAlterations);
							if (newEntry.numSamples() >= t && !nodesHash_buckets[newHash].count(*newEntry.samples)) {	// Number of patients is still high enough and the subgraph is not a duplicate
								newEntry.isValid = true; // Need to do this because the base subnetwork's flag may have got marked as invalid in the 'if' below, before all neighbours got considered.
								newEntry.nodes.push_back(neighbourIdx);
								newEntry.nodeColourIdx.push_back(alterationIndex);
								candidateSubnetworks[cycle].push_back(newEntry);
								nodesHash_buckets[newHash].insert(*newEntry.samples);
								if (!lostSamples && subnetInfo.isValid) {	// The newly identified subnetwork is a sample-wise-lossless extension, making the base subnetwork redundant
									subnetInfo.isValid = false;
									numContained++;
								}
							}
						}
					}
				}
			}
		}
		totalNumSubgraphs += candidateSubnetworks[cycle].size();
	}
	delete [] nodesHash_buckets;
	delete nodesHash;
	delete visitedNode;
	fprintf(stderr, "\nConstructed all candidate subnetworks. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf(stderr, "Total amount of subgraphs of all sizes up to %d which are recurrent in at least %d patients is %llu.\n", S, t, totalNumSubgraphs);
	fprintf(stderr, "%llu subgraphs are contained in a subgraph of larger size with the same patients, and are discarded.\n", numContained);

	//
	//	STAGE 4: Filtering out candidate subnetworks based on their colour properties.
	//
	fprintf(stderr, "Filtering out homogenous (colour-wise) subnetworks.\n"); timerStart = clock();
	llu numImproperlyColoured = 0;
	llu numConsidered = 0;
	llu * colourCount = new llu [alterations.indices.size() + 1];
	memset(colourCount, 0, sizeof(colourCount[0]) * (1 + alterations.indices.size()));
	for (llu sizeIdx = 0, currentSubgraphCounter = 0; sizeIdx < S; sizeIdx++) {
		for (int i = 0, lastProg = 0; i < candidateSubnetworks[sizeIdx].size(); i++) {
			currentSubgraphCounter++;
			int progress = 1000 * double(currentSubgraphCounter) / double(totalNumSubgraphs);
			if (progress > lastProg) {
				fprintf(stderr, "\r%.1lf%%", double(progress)/10.0);
				lastProg = progress;
			}
			SubnetworkEntry & subnetInfo = candidateSubnetworks[sizeIdx][i];
			if (subnetInfo.isValid) {
				numConsidered++;
				// vector<llu> nodeColourMasks = subnetInfo.getNodeColourBitmaskVector(geneAlterations);
				llu expressionOutlierBitmask = llu(1) << alterations.indices["EXPROUT"];
				// int numNonOutlier = 0;
				llu subnetColourBitmask = 0;
				for (llu colourIdx : subnetInfo.nodeColourIdx) {
					// if ((expressionOutlierBitmask & nodeColourBitmask) == 0)
						// numNonOutlier++;
					subnetColourBitmask |= llu(1) << colourIdx;
				}
				int numColours = __builtin_popcountll(subnetColourBitmask);
				colourCount[numColours]++;
				// if (numNonOutlier < 2) {
				if (numColours < minColours) {
					subnetInfo.isValid = false;
					numImproperlyColoured++;
				}
			}
		}
	}
	fprintf(stderr, " (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf(stderr, "%llu subgraphs were considered.\n", numConsidered);
	for (int i = 1; i <= alterations.indices.size(); i++) {
		fprintf(stderr, "\t%llu subnetworks have %d colour%s among their nodes.\n", colourCount[i], i, i > 1 ? "s" : "");
	}
	if (minColours > 1) {
		fprintf(stderr, "%llu subgraphs do not have at least %d differently-coloured nodes, and are discarded.\n", minColours, numImproperlyColoured);
	}
	delete colourCount;

	//
	//	STAGE 4.5: constructing array of all valid subnetworks.
	//

	llu numProperSubgraphs = totalNumSubgraphs - numContained - numImproperlyColoured;
	SubnetworkEntry ** properSubgraphs = new SubnetworkEntry * [numProperSubgraphs];
	numProperSubgraphs = 0;
	for (int sizeIdx = 0; sizeIdx < S; sizeIdx++) {
		for (int i = 0; i < candidateSubnetworks[sizeIdx].size(); i++) {
			if (candidateSubnetworks[sizeIdx][i].isValid) {
				properSubgraphs[numProperSubgraphs] = & candidateSubnetworks[sizeIdx][i];
				numProperSubgraphs++;
			}
		}
	}
	fprintf(stderr, "%llu proper subgraphs are considered.\n", numProperSubgraphs);

	//
	//	STAGE 5: Extending candidate subnetworks to include samples in which there aren't all exact matches, but there aren't colour conflicts either. 
	//

	if (errorRate >= 1.0/S) {
		timerStart = clock();
		llu numSubnetworksExtended = 0;
		llu numSamplesAdded = 0;
		fprintf(stderr, "Extending subnetworks to include %d%% errors...\n", int(errorRate * 100));
		for (int i = 0, lastProg = 0; i < numProperSubgraphs; i++) {
			int progress = 1000 * double(i + 1) / double(numProperSubgraphs);
			if (progress > lastProg) {
				fprintf(stderr, "\r%.1lf%%", double(progress)/10.0);
				lastProg = progress;
			}
			SubnetworkEntry * subnetInfo = properSubgraphs[i];
			if (subnetInfo->isValid) {
				int numSamplesBefore = subnetInfo->numSamples();
				subnetInfo->extendSubnetworkWithError(geneAlterations, samples, errorRate);
				int numSamplesAfter = subnetInfo->numSamples();
				if (numSamplesAfter > numSamplesBefore) {
					numSubnetworksExtended++;
					numSamplesAdded += numSamplesAfter - numSamplesBefore;
				}
			}
		}
		fprintf(stderr, "\rDone. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
		fprintf(stderr, "%llu subnetworks have been extended.\n", numSubnetworksExtended);
		fprintf(stderr, "Average number of samples added is %.1lf\n", double(numSamplesAdded) / numSubnetworksExtended);
	}
	
	vector<int> * nodeCover = new vector<int> [G.V * samples.indices.size()];
	timerStart = clock();
	for (llu i = 0; i < numProperSubgraphs; i++) {
		SubnetworkEntry * subnetInfo = properSubgraphs[i];
		Bitmask tempmask(subnetInfo->samples);
		while (tempmask.getSize()) {
			int sampleIdx = tempmask.extractLowestOrderSetBitIndex();
			for (int nodeIdx : subnetInfo->nodes) {
				llu commonIdx = sampleIdx * G.V + nodeIdx;
				nodeCover[commonIdx].push_back(i);
			}
		}
	}
	fprintf(stderr, "\tCalculated node covers. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);

	IloEnv env;
	IloModel model(env);
	IloExpr objective(env);
	IloBoolVarArray X(env, numProperSubgraphs);
	fprintf(stderr, "\tConstructed X variables.\n");
	IloBoolVarArray C(env, subnetworkSeeds.size());
	fprintf(stderr, "\tConstructed C variables.\n");

	// Maximize the number of covered nodes
	timerStart = clock();
	for (llu i = 0; i < subnetworkSeeds.size(); i++) {
		int sampleIdx		= subnetworkSeeds[i].first;
		int nodeIdx			= subnetworkSeeds[i].second;
		llu commonIdx		= sampleIdx * G.V + nodeIdx;
		if (!nodeCover[commonIdx].empty())
			objective += C[i];
	}
	model.add( IloMaximize(env, objective) );
	fprintf(stderr, "\tConstructed the objective function. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);

	// No element can be covered without a set that contains it being picked
	llu uncoveredNodeCnt = 0;
	timerStart = clock();
	llu averageCoverage = 0;
	unordered_set<int> samplesWithNodesThatCanBeCovered;
	fprintf(stderr, "\tAdding constraints: \n");
	for (int k = 0, lastProg = 0; k < subnetworkSeeds.size(); k++) {
		int progress = 1000 * double(k + 1) / double(subnetworkSeeds.size());
		if (progress > lastProg) {
			fprintf(stderr, "\r\t%.1lf%%", double(progress)/10.0);
			lastProg = progress;
		}
		int node_sampleIdx		= subnetworkSeeds[k].first;
		int node_nodeIdx		= subnetworkSeeds[k].second;
		int node_commonIdx		= node_sampleIdx * G.V + node_nodeIdx;
		// llu subnetworkCnt = 0;
		if (nodeCover[node_commonIdx].size()) {
			samplesWithNodesThatCanBeCovered.insert(node_sampleIdx);
			IloExpr e(env);
			for (int subgraphIdx : nodeCover[node_commonIdx]) {
				e += X[subgraphIdx];
			}
			model.add(e >= C[k]);
			averageCoverage += nodeCover[node_commonIdx].size();
		}
		else
			uncoveredNodeCnt++;
	}
	fprintf(stderr, "\r\tAdded coverage constraints. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);

	// Number of sets that we can pick is at most K
	{
		timerStart = clock();
		IloExpr e(env);
		for (llu i = 0; i < numProperSubgraphs; i++) {
			e += X[i];
		}
		model.add(e <= K);
		fprintf(stderr, "\r\tAdded set number constraint. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); 
	}

	fprintf(stderr, "\t%llu node%s cannot be covered by even a single subnetwork that is recurrent in at least %d patients, and are removed from the model.\n", uncoveredNodeCnt, uncoveredNodeCnt != 1 ? "s" : "", t);
	fprintf(stderr, "\tThe remaining %llu nodes have on average %.1lf subgraphs that can cover them.\n", subnetworkSeeds.size() - uncoveredNodeCnt, averageCoverage / double(subnetworkSeeds.size() - uncoveredNodeCnt));
	fprintf(stderr, "\tThose nodes belong to a total of %d samples, out of %d.\n", samplesWithNodesThatCanBeCovered.size(), samples.indices.size());
	fprintf(stderr, "\tRunning ILP to find %d subnetworks that cover as many of the remaining nodes as possible.\n", K);
	
	try {
		char command[1000];
		// char fullFolder[1000];
		// sprintf(fullFolder, "./%s_s%d_t%d_k%d", folderName, S, t, K);
		string outFolder = string(folderName);
		// sprintf(command, "rm -f -r %s", outFolder.c_str());
		// system(command);
		// sprintf(command, "mkdir -p %s", outFolder.c_str());
		// system(command);
		string outModel						= outFolder + "/ilp_model.lp";
		string outSol						= outFolder + "/ilp_solution.txt";
		string outSubnetworks				= outFolder + "/subnetworksOverview.txt";
		string outDistributionSize			= outFolder + "/subnetwork_sizes.txt";
		string outDistributionRecurrence	= outFolder + "/subnetwork_recurrence.txt";
		string outSubnFolder				= outFolder + "/subnetworks";
		string outCandidateSamples			= outFolder + "/candidateSamples.txt";
		string outCoveredSamples			= outFolder + "/coveredSamples.txt";
		FILE * foutSamples = fopen(outCandidateSamples.c_str(), "w");
		for (int sampleIdx : samplesWithNodesThatCanBeCovered) {
			fprintf(foutSamples, "%s\n", samples.names[sampleIdx].c_str());
		}
		fclose(foutSamples);
		// exit(0);	// STOPS HERE	// STOPS HERE	// STOPS HERE	// STOPS HERE	// STOPS HERE	// STOPS HERE	// STOPS HERE	// STOPS HERE	// STOPS HERE
		sprintf(command, "mkdir -p %s", outSubnFolder.c_str());
		system(command);
		IloCplex cplex(model);
		cplex.exportModel(outModel.c_str());
		fprintf(stderr, "ILP model file written to '%s'.\n", outModel.c_str());
		cplex.setParam(IloCplex::IntParam::Threads, threads);
		// cplex.setParam(IloCplex::NumParam::TiLim, 172800);	// 2 days
		cplex.setParam(IloCplex::NumParam::TiLim, seconds);	// 10 hours
		cplex.solve();
		cplex.writeSolution(outSol.c_str());
		fprintf(stderr, "ILP solution file written to '%s'.\n", outSol.c_str());
		vector<int> subnetworksizes;
		vector<int> subnetworkrecurrence;
		FILE * fout = fopen(outSubnetworks.c_str(), "w");
		unordered_set<int> subnetworkNodes;
		vector<int> visited(G.V, 0);
		llu visitedIdx = 0;
		unordered_set<int> samplesWithCoveredNodes;
		for (llu i = 0, subnIdx = 0; i < numProperSubgraphs; i++) {
			if (cplex.getValue(X[i]) != 0) {
				SubnetworkEntry * subnetInfo = properSubgraphs[i];
				subnIdx++;
				fprintf(fout, "Subnetwork\t%llu\n", subnIdx);
				subnetInfo->print(geneAlterations, samples, alterations, G.chrArm, G.nodeNames, fout);
				fprintf(fout, "\n");
				char filename[1000];
				sprintf(filename, "%s/%llu.edges", outSubnFolder.c_str(), subnIdx);
				FILE * foutEdges = fopen(filename, "w");
				sprintf(filename, "%s/%llu.adj", outSubnFolder.c_str(), subnIdx);
				FILE * foutAdj = fopen(filename, "w");
				sprintf(filename, "%s/%llu.nodes", outSubnFolder.c_str(), subnIdx);
				FILE * foutNodes = fopen(filename, "w");
				sprintf(filename, "%s/%llu.samples", outSubnFolder.c_str(), subnIdx);
				FILE * foutSamples = fopen(filename, "w");
				subnetworkNodes.clear();
				for (int nodeIdx : subnetInfo -> nodes) {
					subnetworkNodes.insert(nodeIdx);
				}
				Bitmask tempmask(subnetInfo->samples);
				while (tempmask.getSize()) {
					int sampleIdx = tempmask.extractLowestOrderSetBitIndex();
					samplesWithCoveredNodes.insert(sampleIdx);
					fprintf(foutSamples, "%s\n", samples.names[sampleIdx].c_str());
				}
				for (int nodeIdx : subnetInfo -> nodes) {
					visitedIdx++;
					for (int edgeIdx = 0; edgeIdx < G.NSize[nodeIdx]; edgeIdx++) {
						int neighbour = G.N[nodeIdx][edgeIdx];
						if (subnetworkNodes.count(neighbour)) {
							fprintf(foutEdges, "%s %s\n", G.nodeNames[nodeIdx].c_str(), G.nodeNames[neighbour].c_str());
							visited[neighbour] = visitedIdx;
						}
					}
					for (int nodeIdx2 : subnetInfo->nodes) {
						if (nodeIdx != nodeIdx2 && visited[nodeIdx2] == visitedIdx)
							fprintf(foutAdj, "1 ");
						else
							fprintf(foutAdj, "0 ");
					}
					fprintf(foutAdj, "\n");
					fprintf(foutNodes, "%s\n", G.nodeNames[nodeIdx].c_str());
				}
				fclose(foutEdges);
				fclose(foutAdj);
				fclose(foutNodes);
				fclose(foutSamples);
				subnetworksizes.push_back(subnetInfo->nodes.size());
				subnetworkrecurrence.push_back(subnetInfo->numSamples());
			}
		}
		fclose(fout);
		fprintf(stderr, "Subnetwork information written to '%s'.\n", outSubnetworks.c_str());
		fprintf(stderr, "Out of %d samples with nodes that could be covered, %d samples support one of the chosen subnetworks.\n", samplesWithNodesThatCanBeCovered.size(), samplesWithCoveredNodes.size());
		foutSamples = fopen(outCoveredSamples.c_str(), "w");
		for (int sampleIdx : samplesWithCoveredNodes) {
			fprintf(foutSamples, "%s\n", samples.names[sampleIdx].c_str());
		}
		fclose(foutSamples);
		// sort(subnetworksizes.begin(), subnetworksizes.end());
		// sort(subnetworkrecurrence.begin(), subnetworkrecurrence.end());
		fout = fopen(outDistributionSize.c_str(), "w");
		for (int i = 0; i < subnetworksizes.size(); i++) fprintf(fout, "%d\n", subnetworksizes[i]);
		fclose(fout);
		fout = fopen(outDistributionRecurrence.c_str(), "w");
		for (int i = 0; i < subnetworkrecurrence.size(); i++) fprintf(fout, "%d\n", subnetworkrecurrence[i]);
		fclose(fout);/**/
	}
	catch (IloException &ex) {}

	for (int i = 0; i < samples.indices.size(); i++) {
		delete [] node_CCIndex[i];
		for (int j = 0; j < CC_count[i]; j++) delete CC[i][j];
		delete CC[i];
	}
	delete node_CCIndex;
	delete [] CC_count;
	delete CC;/**/
}

int main( int argc, char * argv[] ) {
	printHeader( "MCSC ILP" );
	// INPUT CHECK
	if (argc <= 1) {
		fprintf(stderr, "./mcsc -n [network] -l [alteration profiles] -c [chromosome information; optional] -r [min number of colours in subnetwork] -x [exclude genes; optional] -s [maximum subnetwork size] -t [minimum subgraph recurrence]  -k [number of subnetworks] -e [error; optional] -f [outputFolder] -d [threads] -h [time limit in seconds]\n\n");
		return 0;
	}
	char consoleFlags[] = {'n', 'l', 's', 't', 'f', 'k', 'c', 'x', 'e', 'd', 'h', 'r', 0};
	bool optional[200] = {};
	optional['c'] = true;
	optional['x'] = true;
	optional['e'] = true;
	unordered_map<char, string> consoleParameters;
	for (int i = 1; i < argc; i++) {
		if ( argv[i][0] == '-' && argv[i][1] && i + 1 < argc && argv[i + 1][0] != '-' ) {
			consoleParameters[ argv[i][1] ] = string( argv[i + 1] );
			i++;
		}
	}
	for (char * ptrFlag = consoleFlags; *ptrFlag; ptrFlag++) {
		if ( !optional[*ptrFlag] && !consoleParameters.count(*ptrFlag) ) {
			fprintf(stderr, "\n< Error > Missing value for parameter '%c'. Exiting program.\n", *ptrFlag);
			exit(0);
		}
	}
	int maxSubnetworkSize;
	int minSubnetworkRecurrence;
	int K;
	double errorRate = 0;
	int threads;
	int seconds;
	int minColours;
	char folderName[1000] = {};
	sscanf(consoleParameters['s'].c_str(), "%d", &maxSubnetworkSize);
	sscanf(consoleParameters['t'].c_str(), "%d", &minSubnetworkRecurrence);
	sscanf(consoleParameters['k'].c_str(), "%d", &K);
	sscanf(consoleParameters['d'].c_str(), "%d", &threads);
	sscanf(consoleParameters['h'].c_str(), "%d", &seconds);
	sscanf(consoleParameters['r'].c_str(), "%d", &minColours);
	if (consoleParameters.count('e')) {
		sscanf(consoleParameters['e'].c_str(), "%lf", &errorRate);
		fprintf(stderr, "Error rate set to %.2lf\n", errorRate);
	}
	sscanf(consoleParameters['f'].c_str(), "%s", folderName);
	// Create directory structure for the output
	char command[1000];
	char fullFolder[1000];
	if (consoleParameters.count('e'))
		sprintf(fullFolder, "./%s_s%d_t%d_k%d_e%.2lf_r%d", ("output/" + string(folderName)).c_str(), maxSubnetworkSize, minSubnetworkRecurrence, K, errorRate, minColours);
	else 
		sprintf(fullFolder, "./%s_s%d_t%d_k%d_r%d", ("output/" + string(folderName)).c_str(), maxSubnetworkSize, minSubnetworkRecurrence, K, minColours);
	string outFolder = string(fullFolder);
	sprintf(command, "rm -f -r %s", outFolder.c_str());
	system(command);
	sprintf(command, "mkdir -p %s", outFolder.c_str());
	system(command);
	string outLog = outFolder + "/run.log";
	// freopen(outLog.c_str(), "w", stderr);
	printHeader("Reading Input");
	readUndirectedNetwork( consoleParameters['n'].c_str() );
	if (consoleParameters.count('c'))
		readChromosomeInfo( consoleParameters['c'].c_str() );
	readAlterationProfiles( consoleParameters['l'].c_str() );
	if (consoleParameters.count('x'))
		readExcludeInfo( consoleParameters['x'].c_str() );
	printHeader("Solving the problem");
	runSolver(maxSubnetworkSize, minSubnetworkRecurrence, K, errorRate, outFolder.c_str(), threads, seconds, minColours);
	return 0;
}
