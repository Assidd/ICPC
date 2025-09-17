#include "algotester.h"
#include <bits/stdc++.h>

using namespace std;

long long f(vector<long long>& cost, vector<vector<int>>& flowVerts, vector<vector<int>>& pipelines)
{
	long long result = 0;
    for (long long i = 0; i < int(pipelines.size()); ++i) 
	{
		vector<int>& flows = pipelines[i];
		for (long long j = 0; j < int(flows.size()); ++j)
		{
			vector<int>& verts = flowVerts[flows[j]];
			if (j == 0)
				result += cost[verts[0]];
			for (long long k = 1; k < int(verts.size()); ++k)
				result += cost[verts[k]];
		}
	}
	return result;
}


int main(int argc, char* argv[])
{
    auto [in, out, ans] = initScorer(argc, argv);
    
    // ---------- Read N, M, K ----------
    long long N = in.readInt();
    long long M = in.readInt();
    long long K = in.readInt();

    // ---------- Read vertices ----------
    vector<long long> cost(N);
    
    for (long long i = 0; i < N; ++i) 
	{
        cost[i] = in.readInt();
        in.readInt();
        int pz = in.readInt();
        for (long long j = 0; j < pz; ++j)
		    in.readInt();
	}
    // ---------- Read flows ----------
    vector<vector<int>> flowVerts(M);
    for (long long i = 0; i < M; ++i) 
	{
        in.readInt();
        long long L = in.readInt();

		flowVerts[i].resize(L);
        for (long long j = 0; j < L; ++j) 
			flowVerts[i][j] = in.readInt();
    }
    
    // ---------- Read pipelines ----------
    vector<vector<int>> pipelines(K);

    for (long long i = 0; i < K; ++i)
	{
        long long L = in.readInt();
        pipelines[i].resize(L);
        for (long long j = 0; j < L; ++j) 
			pipelines[i][j] = in.readInt();
    }
    
    // ---------- Calculate U_full ----------
    long long UFull = f(cost, flowVerts, pipelines);
    
    
    // ---------- Read Sub pipelines ----------
	int USubLen = out.readInt();
	
    vector<vector<int>> pipelinesSub(USubLen);
    for (long long i = 0; i < USubLen; ++i)
    {
		int idx = out.readInt();
		pipelinesSub[i] = pipelines[idx];
	}

    // ---------- Calculate U_sub ----------
    long long USub = f(cost, flowVerts, pipelinesSub);
	
    // ---------- Read New pipelines ----------
	int UNewLen = out.readInt();
	
    vector<vector<int>> pipelinesNew(UNewLen);
    
    for (long long i = 0; i < UNewLen; ++i)
	{
        long long L = out.readInt();
        pipelinesNew[i].resize(L);
        for (long long j = 0; j < L; ++j) 
			pipelinesNew[i][j] = out.readInt();
    }
	
    // ---------- Calculate U_new ----------
    long long UNew = f(cost, flowVerts, pipelinesNew);
    
    
    // ---------- Calculate scores ----------
    double score1 = max(0.0, double(UFull - USub) / UFull);
    double score2 = max(0.0, double(UFull - UNew) / UFull);
	
    long long result = 10'000'000 * (0.2 * score1 + 0.8 * score2);
    
    cout << result << '\n';
    
    return 0;
}


