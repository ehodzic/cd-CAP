#include <cstdio>
#include <cstring>
#include <ctime>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>
#include <ilcplex/ilocplex.h>
#include <algorithm>
using namespace std;

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

struct Graph {
	int V, E;
	int * NSize;
	int ** N;
	string * nodeNames;
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

unordered_map<int, unsigned int> * geneAlterations;	// geneAlterations[ i ][ j ] = c means that "gene i has colour c in patient j". The colours are bitmasks (so supporting max 32 different alteration types).
vector< pair<int,int> > subnetworkSeeds;

/*
	Reads the input -n parameter as a collection of undirected edges (pairs of node names, separated by whitespace) and stores the information in global Graph object G.
*/
void readUndirectedNetwork(const char * filename) {
	fprintf(stderr, "Reading the network... ");
	int timerStart = clock();
	unordered_set<pair<string, string>, stringPairHash> uniqueEdges;
	char u[1000], v[1000];
	FILE * fin = fopen(filename, "r");
	G.V = 0;
	while (fscanf(fin, "%s%s", u, v) == 2) {
		char * a = u, * b = v;
		if (strcmp(a, b) > 0) swap(a, b);
		else if (strcmp(a, b) == 0) continue;
		pair<string, string> e = make_pair(string(a), string(b));
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

	G.NSize = new int [G.V];
	G.N = new int * [G.V];
	G.nodeNames = new string[G.V];
	memset(G.NSize, 0, sizeof(G.NSize[0]) * G.V);

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
	Reads the input -l parameter as a collection of "sample gene alterationType" triples, separated by whitespace.
*/
void readAlterationProfiles(const char * filename) {
	fprintf(stderr, "Reading the alteration profiles... ");
	int timerStart = clock();
	char sample[1000], gene[1000], alterationType[1000];
	FILE * fin = fopen(filename, "r");
	while (fscanf(fin, "%s%s%s", sample, gene, alterationType) == 3) {
		if ( G.nodeIndices.find( string(gene) ) == G.nodeIndices.end() ) continue;
		if ( samples.indices.find( string(sample) ) == samples.indices.end() ) {
			int idx = samples.indices.size();
			samples.indices[ string(sample) ] = idx;
		}
		if ( genes.indices.find( string(gene) ) == genes.indices.end() ) {
			int idx = genes.indices.size();
			genes.indices[ string(gene) ] = idx;
		}
		if ( alterations.indices.find( string(alterationType) ) == alterations.indices.end() ) {
			int idx = alterations.indices.size();
			alterations.indices[ string(alterationType) ] = idx;
		}
	}
	samples.names = new string[ samples.indices.size() ];
	genes.names = new string[ genes.indices.size() ];
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
	geneAlterations = new unordered_map<int, unsigned int> [ G.V ];
	while (fscanf(fin, "%s%s%s", sample, gene, alterationType) == 3) {
		if ( G.nodeIndices.find( string(gene) ) == G.nodeIndices.end() ) continue;
		int geneIndex = G.nodeIndices[string(gene)];
		int sampleIndex = samples.indices[sample];
		// int subnetworkIndex = G.nodeIndices[string(gene)] + sampleIndex * G.V;
		subnetworkSeeds.push_back(make_pair(sampleIndex, G.nodeIndices[string(gene)]));
		int alterationIndex = alterations.indices[alterationType];
		if ( geneAlterations[geneIndex].find(sampleIndex) == geneAlterations[geneIndex].end() ) geneAlterations[geneIndex][sampleIndex] = 0;
		int bitMask = (1 << alterationIndex);
		geneAlterations[geneIndex][sampleIndex] |= bitMask;
	}
	fclose(fin);
	fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC);
	fprintf( stderr, "\tThere are %d samples, with a total of %d genes, harboring %d different alterations.\n", samples.indices.size(), genes.indices.size(), alterations.indices.size() );
	fprintf(stderr, "\tThere are %d possible subnetwork seeds.\n", subnetworkSeeds.size());
}

void runILPSolver(int maxSubnetworkSize) {
	fprintf(stderr, "Finding coloured patient-specific connected components... ");
	int timerStart = clock();
	int * node_CCIndex = new int [G.V];	// node_CCIndex[nodeIdx] = CCIndex
	int CC_count = 0;	// Number of connected components
	int * nodeStack = new int [G.V + 1];
	Subgraph ** CC;

	memset(node_CCIndex, -1, sizeof(node_CCIndex[0])*G.V);
	int nodeStackSize = 0;
	for (int j = 0; j < G.V; j++) {
		if (node_CCIndex[j] == -1) {	// node is coloured and not assigned to a connected component
			node_CCIndex[j] = CC_count++;
			nodeStack[nodeStackSize++] = j;
			while (nodeStackSize) {
				int node = nodeStack[--nodeStackSize];
				for (int j1 = 0; j1 < G.NSize[node]; j1++) {
					int const & neighbour = G.N[node][j1];
					if (node_CCIndex[neighbour] == -1) {	// neighbour is coloured and not assigned to a connected component
						node_CCIndex[neighbour] = node_CCIndex[node];
						nodeStack[nodeStackSize++] = neighbour;
					}
				}
			}
		}
	}
	/***************************************************************************************************************************
	 * Computed number of connected coloured components in the graph. (CC_count)											   *
	 * For every node, assigned index of the connected coloured component it belongs to in the graph. (node_CCIndex[nodeIdx])  *
	 ***************************************************************************************************************************/
	int * CCSizes = new int [G.V + 1];
	memset(CCSizes, 0, CC_count * sizeof(CCSizes[0]));
	for (int j = 0; j < G.V; j++) if (node_CCIndex[j] != -1) CCSizes[ node_CCIndex[j] ]++;
	for (int i = 1; i < CC_count; i++) CCSizes[i] += CCSizes[i - 1];
	int totalNodes = CCSizes[CC_count - 1];
	for (int j = 0; j < G.V; j++) if (node_CCIndex[j] != -1) nodeStack[ --CCSizes[ node_CCIndex[j] ] ] = j;
	/***************************************************************************************
	 * Partitioned nodes of the sample network based on coloured connected component index *
	 ***************************************************************************************/
	CC = new Subgraph * [CC_count];
	for (int i = 0, CCIndex = 0; i < totalNodes; i++, CCIndex++) {
		int j = i + 1;
		while (j < totalNodes && node_CCIndex[nodeStack[j]] == node_CCIndex[nodeStack[i]]) j++;
		CC[CCIndex] = new Subgraph(&G, nodeStack + i, j - i);
		i = j - 1;
	}
	/****************************************************************
	 * Constructed subgraphs based on connected coloured components *
	 ****************************************************************/
	delete [] nodeStack;
	delete [] CCSizes;
	fprintf(stderr, "done. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); timerStart = clock();
	
	IloEnv env;
	IloModel model(env);
	IloExpr objective(env);
	// node j is part of the optimal subnetwork
	IloBoolVarArray X(env, G.V);
	fprintf(stderr, "\rX created. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); timerStart = clock();

	// patient i is chosen for the optimal subnetwork
	IloBoolVarArray P(env, samples.indices.size());
	fprintf(stderr, "\rP created. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); timerStart = clock();

	// node j is the seed of the optimal subnetwork
	IloBoolVarArray S(env, G.V);
	fprintf(stderr, "\rS created. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); timerStart = clock();

	// patients i and i1 are both chosen for the optimal subnetwork that node j is part of
	IloNumVarArray ** p = new IloNumVarArray * [G.V];
	for (int j = 0; j < G.V; j++) {
		p[j] = new IloNumVarArray [samples.indices.size()];
		for (int i = 0; i < samples.indices.size(); i++) {
			fprintf(stderr, "\rp: %.1lf%%", 100*double(j * samples.indices.size() + i + 1) / double(G.V * samples.indices.size()));
			// if ( (geneAlterations[j][i] & geneAlterations[nodeIdx][sampleIdx]) == 0 ) continue;	// Node isn't coloured in patient i
			p[j][i] = IloNumVarArray(env, samples.indices.size(), 0, 1);
		}
	}
	fprintf(stderr, "\rp created. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); timerStart = clock();

	IloNumVarArray sinkEdges(env, G.V);
	fprintf(stderr, "\rsinkEdges created. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); timerStart = clock();

	IloNumVarArray sourceEdges(env, G.V);
	fprintf(stderr, "\rsinkEdges created. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); timerStart = clock();

	IloNumVarArray * capacities = new IloNumVarArray [ CC_count ];
	for (int CCIndex = 0; CCIndex < CC_count; CCIndex++) {
		;
	}

	for (int k = 0; k < subnetworkSeeds.size(); k++) {
		fprintf(stderr, "\rCapacities: %.1lf%%", 100*double(k + 1) / double(subnetworkSeeds.size()));
		int sampleIdx = subnetworkSeeds[k].first;
		int nodeIdx = subnetworkSeeds[k].second;
		int CCIndex = node_CCIndex[sampleIdx][nodeIdx];
		if (CCIndex == -1) {
			fprintf(stderr, "< Error > Seed %d in patient %d has coloured connected component with index -1.\n", sampleIdx, nodeIdx);
			exit(0);
		}
		int CCSize = CC[sampleIdx][CCIndex] -> V;
		capacities[k] = new IloNumVarArray [CCSize];
		for (int j_sub = 0; j_sub < CCSize; j_sub++) {
			int degree = CC[sampleIdx][CCIndex] -> degrees[j_sub];
			capacities[k][j_sub] = IloNumVarArray(env, degree, 0, S);
		}
	}
	fprintf(stderr, "\rCapacities created. (%.2lf seconds)\n", double(clock() - timerStart) / CLOCKS_PER_SEC); timerStart = clock();
	// IloNumVarArray * capacities = new IloNumVarArray [G.V];
	// for (int j = 0; j < G.V; j++) {
	// 	capacities[j] = IloNumVarArray(env, G.NSize[j]);
	// 	for (int nIdx = 0; nIdx < G.NSize[j]; nIdx++) {
	// 		// int neighbor = G.N[j][nIdx];
	// 		char name[ 100 ];
	// 		sprintf(name, "C;%d;%d", j, nIdx);
	// 		capacities[j][nIdx] = IloNumVar(env, 0, G.V, name);
	// 	}
	// }

	// Objective: Maximize sum of patient-pair-wise scores per node
	for (int j = 0; j < G.V; j++) {
		// if ( j >= 1779 ) fprintf(stderr, "%d ", j);
		for (int i = 0; i < samples.indices.size(); i++) {
			bool iHasColor = geneAlterations[j].count(i);
			// if ( j >= 1779 ) fprintf(stderr, "%d ", iHasColor);
			for (int i1 = i+1; i1 < samples.indices.size(); i1++) {
				bool i1HasColor = geneAlterations[j].count(i1);
				// if ( j >= 1779 ) fprintf(stderr, "%d ", i1HasColor);
				if ( iHasColor && i1HasColor && (geneAlterations[j][i] & geneAlterations[j][i1]) ) objective += p[i][i1][j];
				else if ( iHasColor && i1HasColor && (geneAlterations[j][i] & geneAlterations[j][i1] == 0) ) objective -= p[i][i1][j];
				else ;	// the case of just one of the nodes being colorless is not penalized
			}
		}
		// cin.get();
	}
	// fprintf(stderr, "1337\n");
	model.add( IloMaximize(env, objective) );
	// fprintf(stderr, "1337\n");
	fprintf(stderr, "\tConstructed the objective function.\n");


	// Constraints (14), (15), (16) and (17)
	for (int i = 0; i < samples.indices.size(); i++) {
		for (int i1 = 0; i1 < samples.indices.size(); i1++) {
			for (int j = 0; j < G.V; j++) {
				// IloExpr e(env);
				// e += X[i1][j1][k];
				model.add( p[i][i1][j] - P[i] <= 0 );
				model.add( p[i][i1][j] - P[i1] <= 0 );
				model.add( p[i][i1][j] - X[j] <= 0 );
				model.add( p[i][i1][j] - P[i] - P[i1] - X[j] + 2 <= 0 );
			}
		}
	}
	fprintf(stderr, "\tConstructed constraint (14), (15), (16) and (17).\n");

	// Constraint (18) - Only one node can be the seed of the optimal subnetwork.
	{
		IloExpr e(env);
		for (int j = 0; j < G.V; j++) {
			e += S[j];
		}
		model.add( e == 1 );
	}
	fprintf(stderr, "\tConstructed constraint (18).\n");

	// constraint (19) - The seed is trivially chosen to be part of the subnetwork.
	for (int j = 0; j < G.V; j++) {
		model.add( X[j] >= S[j] );
	}
	fprintf(stderr, "\tConstructed constraint (19).\n");

	// constraint (20) - Sink edges all have capacity one. Handled in sinkEdges construction.
	fprintf(stderr, "\tConstructed constraint (20).\n");

	// constraint (21) - We let flow from the source go only through the node marked as the seed.
	for (int j = 0; j < G.V; j++) {
		model.add( sourceEdges[j] - G.V * S[j] <= 0 );
	}
	fprintf(stderr, "\tConstructed constraint (21).\n");

	// constraint (22) - If a node is chosen to be part of the subnetwork, there must be flow entering it.
	for (int j = 0; j < G.V; j++) {
		IloExpr e(env);
		for (int nIdx = 0; nIdx < G.NSize[j]; nIdx++) {
			// int neighbor = G.N[j][nIdx];
			e += capacities[j][nIdx];
		}
		model.add( e >= X[j] );
	}
	fprintf(stderr, "\tConstructed constraint (22).\n");

	// constraint (23) - And if there is flow entering a node, then it must be part of the subnetwork.
	for (int j = 0; j < G.V; j++) {
		IloExpr e(env);
		for (int nIdx = 0; nIdx < G.NSize[j]; nIdx++) {
			// int neighbor = G.N[j][nIdx];
			e += capacities[j][nIdx];
		}
		model.add( e <= G.V * X[j] );
	}
	fprintf(stderr, "\tConstructed constraint (23).\n");

	// constraint (24): Flow must be conserved.
	for (int j = 0; j < G.V; j++) {
		IloExpr e(env);
		e += sourceEdges[j];
		for ( auto edge : G.incomingEdges[j] ) {
			int sourceIdx = edge.first;
			int edgeOffset = edge.second;
			e += capacities[sourceIdx][edgeOffset];
		}
		for (int nIdx = 0; nIdx < G.NSize[j]; nIdx++) {
			// int neighbor = G.N[j][nIdx];
			e -= capacities[j][nIdx];
		}
		e -= sinkEdges[j];
		model.add(e == 0);
	}
	fprintf(stderr, "\tConstructed constraint (24).\n");

	// constrain (25): Maximum subnetwork size is M.
	{
		IloExpr e(env);
		for (int j = 0; j < G.V; j++) {
			e += X[j];
		}
		model.add( e <= maxSubnetworkSize );
	}
	fprintf(stderr, "\tConstructed constraint (25).\n");

	try {
		IloCplex cplex(model);
		cplex.exportModel("model.lp");
		fprintf(stderr, "ILP model file written to 'model.lp'.\n");
		cplex.solve();
		cplex.writeSolution("solution.txt");
		fprintf(stderr, "ILP solution file written to 'solution.txt'.\n");
		FILE * fout = fopen("output.txt", "w");
		for (int j = 0; j < G.V; j++) {
			if (cplex.getValue( S[j] ) == 1 ) {
				fprintf(fout, "Seed: %s\n", G.nodeNames[j].c_str());
				break;
			}
		}
		fprintf(fout, "Genes:");
		for (int j = 0; j < G.V; j++) {
			if (cplex.getValue( X[j] ) == 1 ) fprintf(fout, "\t%s", G.nodeNames[j].c_str());
		}
		fprintf(fout, "\nPatients:");
		for (int i = 0; i < samples.indices.size(); i++) {
			if (cplex.getValue( P[i] ) == 1 ) fprintf(fout, "\t%s", samples.names[i].c_str());
		}
		fprintf(fout, "\n");
		fclose(fout);
	}
	catch (IloException &ex) {}
}

int main( int argc, char * argv[] ) {
	printHeader( "MSC-NCI-motif beta" );
	// INPUT CHECK
	if (argc <= 1) {
		fprintf(stderr, "./motif -n [network] -l [alteration profiles] -s [maximum subnetwork size]\n\n");
		return 0;
	}
	char consoleFlags[] = {'n', 'l', 's', 0};
	unordered_map<char, string> consoleParameters;
	for (int i = 1; i < argc; i++) {
		if ( argv[i][0] == '-' && argv[i][1] && i + 1 < argc && argv[i + 1][0] != '-' ) {
			consoleParameters[ argv[i][1] ] = string( argv[i + 1] );
			i++;
		}
	}
	for (char * ptrFlag = consoleFlags; *ptrFlag; ptrFlag++) {
		if ( consoleParameters.find(*ptrFlag) == consoleParameters.end() ) {
			fprintf(stderr, "\n< Error > Missing value for parameter '%c'. Exiting program.\n", *ptrFlag);
			exit(0);
		}
	}
	int maxSubnetworkSize;
	sscanf(consoleParameters['s'].c_str(), "%d", &maxSubnetworkSize);
	printHeader("Reading Input");
	readUndirectedNetwork( consoleParameters['n'].c_str() );
	readAlterationProfiles( consoleParameters['l'].c_str() );
	printHeader("Solving the problem");
	runILPSolver(maxSubnetworkSize);
	return 0;
}
