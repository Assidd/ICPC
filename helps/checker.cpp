#include "algotester.h"
#include <bits/stdc++.h>

using namespace std;

void f(int N, int M, int K, const vector<long long>& w, const vector<long long>& loc, const vector<vector<int>>& flowVerts, const vector<vector<int>>& pipelines, vector<unordered_set<int>> preconditions)
{
	assert(M == w.size());
	assert(M == flowVerts.size());
	assert(K == pipelines.size());
	
    vector<int> flowStart(M), flowEnd(M);
    vector<long long> flowUseCnt(M, 0);
    for (int i = 0; i < M; ++i)
	{
        flowStart[i] = flowVerts[i].front();
        flowEnd[i]   = flowVerts[i].back();
	}
    vector<int> used(N, false);
    vector<int> last(N, -1);

	for (long long pi = 0; pi < K; ++pi)
	{
        const auto &ids = pipelines[pi];

        // 1) chaining
        for (size_t j = 0; j < ids.size(); ++j) 
        {
			flowUseCnt[ids[j]]++;
		}
        for (size_t j = 0; j + 1 < ids.size(); ++j) 
		{
            int a = ids[j], b = ids[j + 1];
            check(flowEnd[a] == flowStart[b], "Pipeline chaining error: end of flow != start of next flow");
        }

        // 2) expand vertex walk (avoid duplicating junction head)
        vector<int> walk;
        if (!ids.empty()) 
		{
            walk.insert(walk.end(), flowVerts[ids[0]].begin(), flowVerts[ids[0]].end());
            for (size_t j = 1; j < ids.size(); ++j) 
			{
                const auto &fv = flowVerts[ids[j]];
                walk.insert(walk.end(), fv.begin() + 1, fv.end());
            }
        }

        // 3) location rules per occurrence over the full walk
        for (size_t pos = 0; pos < walk.size(); ++pos) 
		{
            int v = walk[pos];
            if (loc[v] == 0) 
			{
                check(pos == 0, "Vertex with location=0 must appear only at the beginning of a pipeline");
            } else if (loc[v] == 1) 
			{
                check(!(pos == 0 || pos + 1 == walk.size()),
                      "Vertex with location=1 cannot be at the beginning or the end of a pipeline");
            } else if (loc[v] == 2) 
			{
                check(pos + 1 == walk.size(),
                      "Vertex with location=2 must appear only at the end of a pipeline");
            }
        }

		// 4) check for preconditions
		for (size_t pos = 0; pos < walk.size(); pos++)
		{
			last[walk[pos]] = pos;
		}
		
		for (int pos = 0; pos < int(walk.size()); pos++)
		{
			int u = walk[pos];
			if (last[u] != pos)
			{
				used[u] = true;
				continue;
			}
			
			vector<int> good;
			for (auto v : preconditions[u])
				if (used[v])
					good.push_back(v);
			for (auto v : good)
				preconditions[u].erase(v);
			used[u] = true;
		}
		
		for (size_t pos = 0; pos < walk.size(); pos++)
		{
			used[walk[pos]] = false;
			last[walk[pos]] = -1;
		}
    }
    
    // for (auto x: preconditions)
    // {
    //     cout << x.first << ' ' << x.second << endl;
    // }
    for (auto precondition : preconditions)
    {
		check(precondition.size() == 0, "Some of the preconditions are not met");
	}

    // ---------- Flow occurrence requirement ----------
    for (long long i = 0; i < M; ++i)
	{
        check(flowUseCnt[i] >= w[i], "Flow occurrence requirement not met");
    }
}


int main(int argc, char* argv[])
{
    auto [in, out, ans] = initChecker(argc, argv);
    
    // ---------- Read N, M, K ----------
    long long N = in.readInt();
    long long M = in.readInt();
    long long K = in.readInt();

    // ---------- Read vertices ----------
    vector<long long> loc(N);
	vector<unordered_set<int>> preconditions(N);
    
    for (long long i = 0; i < N; ++i) 
	{
        in.readInt();
        loc[i] = in.readInt();
        int pz = in.readInt();
        for (long long j = 0; j < pz; ++j)
        {
		    int v = in.readInt();
		    preconditions[i].insert(v);
		}
	}
    // ---------- Read flows ----------
    vector<long long> w(M);
    vector<vector<int>> flowVerts(M);
    for (long long i = 0; i < M; ++i) 
	{
        w[i] = in.readInt();
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
    
    
    
    // ---------- Read Sub pipelines ----------
	int USubLen = out.readInt(1, K, "USubLen");
	
	set<int> s;
    vector<vector<int>> pipelinesSub(USubLen);
    for (long long i = 0; i < USubLen; ++i)
    {
		int idx = out.readInt(0, K - 1);
		pipelinesSub[i] = pipelines[idx];
		s.insert(idx);
	}
	check(s.size() == USubLen, "not unique pipes");

    f(N, M, pipelinesSub.size(), w, loc, flowVerts, pipelinesSub, preconditions);
	
	int UNewLen = out.readInt(1, 1e5, "UNewLen");
	
    // ---------- Read New pipelines ----------
    vector<vector<int>> pipelinesNew(UNewLen);
    
    for (long long i = 0; i < UNewLen; ++i)
	{
        long long L = out.readInt(1, 1000, "Pipe size");
        pipelinesNew[i].resize(L);
        for (long long j = 0; j < L; ++j) 
			pipelinesNew[i][j] = out.readInt(0, M - 1, "flow index");
    }
	
    f(N, M, pipelinesNew.size(), w, loc, flowVerts, pipelinesNew, preconditions);
    
    out.readEof();
    
    return 0;
}




