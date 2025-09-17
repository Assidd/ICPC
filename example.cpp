#include <bits/stdc++.h>
using namespace std;

int main() 
{
	ios::sync_with_stdio(0);
	cin.tie(0);
	

    long long N, M, K;
    cin >> N >> M >> K;

    vector<long long> cost(N);
    vector<int> loc(N);
	map<pair<int, int>, int> preconditions;

    for (long long i = 0; i < N; ++i) 
	{
		cin >> cost[i] >> loc[i];
        long long L;
        cin >> L;
        
        vector<int> s(L);
        for (long long j = 0; j < L; ++j)
        {
			cin >> s[j];
			preconditions[make_pair(i, s[j])] = 0;
		}
    }

    vector<long long> w(M);
    vector<vector<int>> flowVerts(M);

    for (long long i = 0; i < M; ++i) 
	{
		cin >> w[i];
        long long L;
        cin >> L;
        
        flowVerts[i].resize(L);
        for (long long j = 0; j < L; ++j) 
			cin >> flowVerts[i][j];
    }

    vector<vector<int>> pipelines(K);
    vector<int> cntFlows(M);
    vector<int> first(N, -1), last(N, -1);

    for (long long i = 0; i < K; ++i)
	{
        long long L;
        cin >> L;
        pipelines[i].resize(L);
		vector<int> walk;
        for (long long j = 0; j < L; ++j) 
        {
			cin >> pipelines[i][j];
			cntFlows[pipelines[i][j]]++;
			for (auto v : flowVerts[pipelines[i][j]])
				walk.push_back(v);
		}
		
		int sz = walk.size();
        for (long long j = 0; j < sz; ++j) 
			last[walk[j]] = j;
        for (long long j = sz - 1; j >= 0; --j) 
			first[walk[j]] = j;
		for (auto& [p, cnt] : preconditions)
		{
			auto [u, v] = p;
			if (last[u] != -1 && last[v] != -1 && last[u] > first[v])
				cnt++;
		}
        for (long long j = 0; j < sz; ++j) 
			last[walk[j]] = first[walk[j]] = -1;
    }
    
    vector<int> idx;
    vector<int> cntW(M);
    for (long long i = 0; i < K; ++i)
    {
		bool ok = true;
		for (auto f : pipelines[i])
		{
			cntW[f]++;
			ok &= cntFlows[f] - cntW[f] >= w[f];
		}
			
        long long L = pipelines[i].size();
        
		vector<int> walk;
        for (long long j = 0; j < L; ++j) 
        {
			for (auto v : flowVerts[pipelines[i][j]])
				walk.push_back(v);
		}
		
		int sz = walk.size();
        for (long long j = 0; j < sz; ++j) 
			last[walk[j]] = j;
        for (long long j = sz - 1; j >= 0; --j) 
			first[walk[j]] = j;
		for (auto& [p, cnt] : preconditions)
		{
			auto [u, v] = p;
			if (last[u] != -1 && last[v] != -1 && last[u] > first[v])
				ok &= cnt > 1;
		}
        
		if (!ok)
		{
			idx.push_back(i);
		    for (auto f : pipelines[i])
		    	cntW[f] = 0;
		}
		else
		{
		    for (auto f : pipelines[i])
		    {
		    	cntFlows[f] -= cntW[f];
		    	cntW[f] = 0;
		    }
			for (auto& [p, cnt] : preconditions)
			{
				auto [u, v] = p;
				if (last[u] != -1 && last[v] != -1 && last[u] > first[v])
					cnt--;
			}
		}
        for (long long j = 0; j < sz; ++j) 
			last[walk[j]] = first[walk[j]] = -1;
	}

	cout << idx.size() << '\n';
	for (auto i : idx)
		cout << i << ' ';
	cout << '\n';
	cout << idx.size() << '\n';
	for (auto i : idx)
	{
		cout << pipelines[i].size();
		for (auto f : pipelines[i])
			cout << ' ' << f;
		cout << '\n';
	}
   
    return 0;
}





