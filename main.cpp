#include <bits/stdc++.h>
using namespace std;

struct FastScanner {
    static const int BUFSIZE = 1 << 22;
    int idx = 0, size = 0;
    char buf[BUFSIZE];
    inline char getch() {
        if (idx >= size) {
            size = (int)fread(buf, 1, BUFSIZE, stdin);
            idx = 0;
            if (size == 0) return 0;
        }
        return buf[idx++];
    }
    template<typename T>
    bool readInt(T &out) {
        char c = getch();
        if (!c) return false;
        T sign = 1, val = 0;
        while (c != '-' && (c < '0' || c > '9')) { c = getch(); if (!c) return false; }
        if (c == '-') { sign = -1; c = getch(); }
        for (; c >= '0' && c <= '9'; c = getch()) val = val * 10 + (c - '0');
        out = val * sign;
        return true;
    }
};

static inline uint64_t pack_pair(uint32_t v, uint32_t i) { return ((uint64_t)v << 32) | (uint64_t)i; }

struct OA64to32 {
    static constexpr uint64_t EMPTY = ~0ull;
    vector<uint64_t> key;
    vector<uint32_t> val;
    uint32_t fill = 0;
    static inline uint64_t h(uint64_t x) {
        x += 0x9e3779b97f4a7c15ull;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
        return x ^ (x >> 31);
    }
    inline void init(size_t want) {
        size_t cap = 1;
        while (cap < (want << 1)) cap <<= 1;
        key.assign(cap, EMPTY);
        val.assign(cap, 0);
        fill = 0;
    }
    inline void rehash() {
        vector<uint64_t> ok = move(key);
        vector<uint32_t> ov = move(val);
        key.assign(ok.size() << 1, EMPTY);
        val.assign(key.size(), 0);
        fill = 0;
        size_t mask = key.size() - 1;
        for (size_t i = 0; i < ok.size(); ++i) {
            uint64_t k = ok[i];
            if (k == EMPTY) continue;
            uint32_t v = ov[i];
            size_t p = h(k) & mask;
            while (key[p] != EMPTY) p = (p + 1) & mask;
            key[p] = k;
            val[p] = v;
            ++fill;
        }
    }
    inline uint32_t get_or_insert(uint64_t k, uint32_t nxt) {
        size_t mask = key.size() - 1;
        size_t p = h(k) & mask;
        while (true) {
            uint64_t cur = key[p];
            if (cur == EMPTY) {
                key[p] = k;
                val[p] = nxt;
                ++fill;
                if ((uint64_t)fill * 10 > key.size() * 7) rehash();
                return nxt;
            }
            if (cur == k) return val[p];
            p = (p + 1) & mask;
        }
    }
};

struct OASet64 {
    static constexpr uint64_t EMPTY = ~0ull;
    vector<uint64_t> tab;
    vector<uint32_t> used;
    uint32_t fill = 0;
    static inline uint64_t h(uint64_t x) {
        x += 0x9e3779b97f4a7c15ull;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ull;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebull;
        return x ^ (x >> 31);
    }
    inline void reserve_soft(size_t want) {
        size_t cap = 1;
        while (cap < (want << 1)) cap <<= 1;
        tab.assign(cap, EMPTY);
        used.clear();
        used.reserve(cap >> 2);
        fill = 0;
    }
    inline void ensure(size_t want) {
        if (tab.empty()) { reserve_soft(want); return; }
        size_t need = 1;
        while (need < (want << 1)) need <<= 1;
        if (need > tab.size()) reserve_soft(need);
    }
    inline void reset() {
        for (uint32_t pos : used) tab[pos] = EMPTY;
        used.clear();
        fill = 0;
    }
    inline void rehash() {
        vector<uint64_t> old = move(tab);
        tab.assign(old.size() << 1, EMPTY);
        used.clear();
        used.reserve(tab.size() >> 2);
        fill = 0;
        size_t mask = tab.size() - 1;
        for (uint64_t k : old) {
            if (k == EMPTY) continue;
            size_t p = h(k) & mask;
            while (tab[p] != EMPTY) p = (p + 1) & mask;
            tab[p] = k;
            used.push_back((uint32_t)p);
            ++fill;
        }
    }
    inline void insert(uint64_t k) {
        size_t mask = tab.size() - 1;
        size_t p = h(k) & mask;
        while (true) {
            uint64_t cur = tab[p];
            if (cur == EMPTY) {
                tab[p] = k;
                used.push_back((uint32_t)p);
                ++fill;
                if ((uint64_t)fill * 10 > tab.size() * 7) rehash();
                return;
            }
            if (cur == k) return;
            p = (p + 1) & mask;
        }
    }
};

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    FastScanner fs;
    int n, m, k;
    if (!fs.readInt(n)) return 0;
    fs.readInt(m);
    fs.readInt(k);

    vector<uint32_t> cost(n), p_off(n + 1, 0), p_data, p_sz(n, 0);
    p_data.reserve(1 << 20);
    for (int i = 0; i < n; i++) {
        uint32_t ci, li; int L;
        fs.readInt(ci); fs.readInt(li); fs.readInt(L);
        cost[i] = ci;
        p_off[i] = (uint32_t)p_data.size();
        for (int t = 0; t < L; t++) { uint32_t v; fs.readInt(v); p_data.push_back(v); }
        p_sz[i] = (uint32_t)L;
    }
    p_off[n] = (uint32_t)p_data.size();

    vector<uint32_t> w(m), flow_off(m + 1, 0), flow_data;
    flow_data.reserve((size_t)max(1, m * 6));
    vector<uint32_t> flow_first_cost(m, 0), flow_start(m, 0), flow_end(m, 0);
    vector<uint64_t> flow_rest_cost(m, 0);
    vector<uint32_t> flow_pair_w(m, 0);
    for (int j = 0; j < m; j++) {
        uint32_t wj; int L;
        fs.readInt(wj); fs.readInt(L);
        w[j] = wj;
        flow_off[j] = (uint32_t)flow_data.size();
        int fv = -1, lv = -1;
        uint64_t rest = 0;
        uint32_t pairs = 0;
        for (int t = 0; t < L; t++) {
            uint32_t v; fs.readInt(v);
            flow_data.push_back(v);
            if (t == 0) fv = (int)v; else { rest += cost[v]; pairs += p_sz[v]; }
            lv = (int)v;
        }
        flow_off[j + 1] = (uint32_t)flow_data.size();
        flow_first_cost[j] = cost[(size_t)fv];
        flow_rest_cost[j] = rest;
        flow_start[j] = (uint32_t)fv;
        flow_end[j] = (uint32_t)lv;
        flow_pair_w[j] = pairs;
    }

    vector<uint32_t> pipe_off(k + 1, 0), pipe_data;
    pipe_data.reserve((size_t)max(1, k * 10));
    for (int p = 0; p < k; p++) {
        int L; fs.readInt(L);
        pipe_off[p] = (uint32_t)pipe_data.size();
        for (int t = 0; t < L; t++) { uint32_t id; fs.readInt(id); pipe_data.push_back(id); }
    }
    pipe_off[k] = (uint32_t)pipe_data.size();

    vector<uint64_t> pipe_cost(k, 0);
    vector<uint32_t> flow_occ_full(m, 0);
    for (int p = 0; p < k; p++) {
        uint32_t st = pipe_off[p], en = pipe_off[p + 1];
        if (st == en) { pipe_cost[p] = 0; continue; }
        uint64_t c = 0;
        c += flow_first_cost[pipe_data[st]];
        for (uint32_t s = st; s < en; ++s) { uint32_t fj = pipe_data[s]; c += flow_rest_cost[fj]; flow_occ_full[fj]++; }
        pipe_cost[p] = c;
    }

    vector<uint32_t> seen(n, 0); uint32_t token = 1;
    OASet64 uniq; uniq.reserve_soft(1024);
    OA64to32 mp; mp.init(max<size_t>(1 << 20, p_data.size() << 1));
    vector<uint32_t> pair_counts; pair_counts.reserve(1 << 20);
    vector<uint32_t> pair_idx_off(k + 1, 0), pair_idx_data; pair_idx_data.reserve(1 << 20);

    vector<uint32_t> pipe_pair_cap(k, 0);
    for (int p = 0; p < k; p++) {
        uint32_t st = pipe_off[p], en = pipe_off[p + 1], cap = 0;
        for (uint32_t s = st; s < en; ++s) cap += flow_pair_w[pipe_data[s]];
        pipe_pair_cap[p] = max<uint32_t>(cap, 1u);
    }

    for (int p = 0; p < k; p++) {
        uint32_t st = pipe_off[p], en = pipe_off[p + 1];
        pair_idx_off[p] = (uint32_t)pair_idx_data.size();
        if (st == en) continue;
        ++token; if (token >= 0x7fffff00u) { fill(seen.begin(), seen.end(), 0u); token = 1; }
        uniq.ensure(pipe_pair_cap[p]);
        uniq.reset();
        for (uint32_t s = st; s < en; ++s) {
            uint32_t fj = pipe_data[s];
            uint32_t a = flow_off[fj], b = flow_off[fj + 1];
            uint32_t t0 = (s == st) ? 0u : 1u;
            for (uint32_t t = t0; t < b - a; ++t) {
                uint32_t i = flow_data[a + t];
                uint32_t ps = p_off[i], pe = p_off[i + 1];
                for (uint32_t it = ps; it < pe; ++it) {
                    uint32_t v = p_data[it];
                    if (seen[v] == token) uniq.insert(pack_pair(v, i));
                }
                seen[i] = token;
            }
        }
        for (uint32_t pos : uniq.used) {
            uint64_t key = uniq.tab[pos];
            uint32_t id = mp.get_or_insert(key, (uint32_t)pair_counts.size());
            if (id == pair_counts.size()) pair_counts.push_back(1u);
            else pair_counts[id] += 1u;
            pair_idx_data.push_back(id);
        }
    }
    pair_idx_off[k] = (uint32_t)pair_idx_data.size();

    vector<uint32_t> comp_off(k + 1, 0), comp_flow_id;
    vector<uint16_t> comp_flow_cnt;
    comp_flow_id.reserve(pipe_data.size()); comp_flow_cnt.reserve(pipe_data.size());
    {
        vector<uint32_t> tmp;
        for (int p = 0; p < k; p++) {
            comp_off[p] = (uint32_t)comp_flow_id.size();
            uint32_t st = pipe_off[p], en = pipe_off[p + 1];
            if (st == en) continue;
            tmp.assign(pipe_data.begin() + st, pipe_data.begin() + en);
            sort(tmp.begin(), tmp.end());
            for (size_t i = 0; i < tmp.size();) {
                size_t j = i + 1;
                while (j < tmp.size() && tmp[j] == tmp[i]) ++j;
                comp_flow_id.push_back(tmp[i]);
                comp_flow_cnt.push_back((uint16_t)(j - i));
                i = j;
            }
        }
        comp_off[k] = (uint32_t)comp_flow_id.size();
    }

    vector<uint32_t> flow_slack(m, 0);
    for (int j = 0; j < m; j++) {
        uint32_t occ = flow_occ_full[j], need = w[j];
        flow_slack[j] = occ > need ? (occ - need) : 0u;
    }

    vector<uint64_t> met_pair_sum(k, 0), met_flow_sum(k, 0);
    for (int p = 0; p < k; p++) {
        uint64_t s1 = 0;
        for (uint32_t it = pair_idx_off[p]; it < pair_idx_off[p + 1]; ++it) {
            uint32_t id = pair_idx_data[it];
            uint32_t c = pair_counts[id];
            if (c > 0) s1 += (uint64_t)(c - 1u);
        }
        met_pair_sum[p] = s1;
        uint64_t s2 = 0;
        for (uint32_t it = comp_off[p]; it < comp_off[p + 1]; ++it) s2 += (uint64_t)flow_slack[comp_flow_id[it]];
        met_flow_sum[p] = s2;
    }

    auto can_remove = [&](int p, const vector<uint32_t>& cur_flow_occ, const vector<uint32_t>& cur_pair_cnt) -> bool {
        if (pipe_off[p] == pipe_off[p + 1]) return true;
        for (uint32_t it = pair_idx_off[p]; it < pair_idx_off[p + 1]; ++it) {
            uint32_t id = pair_idx_data[it];
            if (cur_pair_cnt[id] <= 1u) return false;
        }
        for (uint32_t it = comp_off[p]; it < comp_off[p + 1]; ++it) {
            uint32_t f = comp_flow_id[it];
            uint32_t c = comp_flow_cnt[it];
            if ((uint64_t)cur_flow_occ[f] < (uint64_t)w[f] + (uint64_t)c) return false;
        }
        return true;
    };
    auto apply_remove = [&](int p, vector<uint32_t>& cur_flow_occ, vector<uint32_t>& cur_pair_cnt) {
        for (uint32_t it = comp_off[p]; it < comp_off[p + 1]; ++it) cur_flow_occ[comp_flow_id[it]] -= comp_flow_cnt[it];
        for (uint32_t it = pair_idx_off[p]; it < pair_idx_off[p + 1]; ++it) {
            uint32_t id = pair_idx_data[it];
            if (cur_pair_cnt[id]) --cur_pair_cnt[id];
        }
    };
    auto apply_add = [&](int p, vector<uint32_t>& cur_flow_occ, vector<uint32_t>& cur_pair_cnt) {
        for (uint32_t it = comp_off[p]; it < comp_off[p + 1]; ++it) cur_flow_occ[comp_flow_id[it]] += comp_flow_cnt[it];
        for (uint32_t it = pair_idx_off[p]; it < pair_idx_off[p + 1]; ++it) {
            uint32_t id = pair_idx_data[it];
            if (cur_pair_cnt[id] < 0xffffu) ++cur_pair_cnt[id];
        }
    };

    auto selection_cost = [&](const vector<int>& v) -> unsigned long long {
        unsigned long long s = 0;
        for (int p : v) s += pipe_cost[p];
        return s;
    };

    auto greedy_cover_param = [&](vector<int>& out, uint32_t PB_SCALE) {
        vector<uint8_t> need_pair(pair_counts.size(), 1u);
        vector<uint64_t> need_flow(m);
        for (int j = 0; j < m; j++) need_flow[j] = w[j];
        uint64_t pairs_left = need_pair.size();
        uint64_t flow_left_sum = 0;
        for (int j = 0; j < m; j++) flow_left_sum += need_flow[j];
        struct Node { uint64_t cost, gain; int id, ver; };
        auto worse = [&](const Node& a, const Node& b) {
            return (__int128)a.cost * b.gain > (__int128)b.cost * a.gain;
        };
        int curv = 1;
        auto gain_pair = [&](int p) -> uint64_t {
            uint64_t g = 0;
            for (uint32_t t = pair_idx_off[p]; t < pair_idx_off[p + 1]; ++t) {
                uint32_t id = pair_idx_data[t];
                if (need_pair[id]) ++g;
            }
            return g;
        };
        auto gain_flow = [&](int p) -> uint64_t {
            uint64_t g = 0;
            for (uint32_t t = comp_off[p]; t < comp_off[p + 1]; ++t) {
                uint32_t f = comp_flow_id[t];
                uint32_t c = comp_flow_cnt[t];
                uint64_t add = min<uint64_t>(need_flow[f], (uint64_t)c);
                g += add;
            }
            return g;
        };
        auto calc_gain = [&](int p) -> uint64_t {
            uint64_t gp = gain_pair(p), gf = gain_flow(p);
            uint64_t W = max<uint64_t>(1, (pairs_left ? (flow_left_sum / max<uint64_t>(1, pairs_left)) : 1));
            __int128 g = (__int128)gf + (((__int128)W * gp * PB_SCALE) >> 10);
            return g > 0 ? (uint64_t)min<__int128>(g, (__int128)ULLONG_MAX) : 0ull;
        };
        vector<Node> heap; heap.reserve(k);
        for (int p = 0; p < k; p++) {
            uint64_t g = calc_gain(p);
            if (g) heap.push_back({pipe_cost[p], g, p, curv});
        }
        make_heap(heap.begin(), heap.end(), worse);
        vector<char> picked(k, 0);
        while ((pairs_left > 0 || flow_left_sum > 0) && !heap.empty()) {
            pop_heap(heap.begin(), heap.end(), worse);
            Node it = heap.back(); heap.pop_back();
            if (picked[it.id]) continue;
            if (it.ver != curv) {
                uint64_t ng = calc_gain(it.id);
                if (ng) {
                    heap.push_back({pipe_cost[it.id], ng, it.id, curv});
                    push_heap(heap.begin(), heap.end(), worse);
                }
                continue;
            }
            uint64_t ng = calc_gain(it.id);
            if (!ng) continue;
            if (ng != it.gain) {
                heap.push_back({it.cost, ng, it.id, curv});
                push_heap(heap.begin(), heap.end(), worse);
                continue;
            }
            int p = it.id;
            picked[p] = 1;
            out.push_back(p);
            for (uint32_t t = pair_idx_off[p]; t < pair_idx_off[p + 1]; ++t) {
                uint32_t id = pair_idx_data[t];
                if (need_pair[id]) {
                    need_pair[id] = 0;
                    --pairs_left;
                }
            }
            for (uint32_t t = comp_off[p]; t < comp_off[p + 1]; ++t) {
                uint32_t f = comp_flow_id[t];
                uint32_t c = comp_flow_cnt[t];
                uint64_t take = min<uint64_t>(need_flow[f], (uint64_t)c);
                if (take) {
                    need_flow[f] -= take;
                    flow_left_sum -= take;
                }
            }
            ++curv;
        }
        if (out.empty()) {
            int best = 0;
            uint64_t bc = pipe_cost[0];
            for (int p = 1; p < k; p++) if (pipe_cost[p] < bc) { bc = pipe_cost[p]; best = p; }
            out.push_back(best);
        }
        vector<uint32_t> occ(m, 0);
        vector<uint32_t> pc(pair_counts.size(), 0);
        for (int p : out) {
            for (uint32_t t = comp_off[p]; t < comp_off[p + 1]; ++t) occ[comp_flow_id[t]] += comp_flow_cnt[t];
            for (uint32_t t = pair_idx_off[p]; t < pair_idx_off[p + 1]; ++t) ++pc[pair_idx_data[t]];
        }
        vector<int> ord = out;
        sort(ord.begin(), ord.end(), [&](int a, int b) {
            if (pipe_cost[a] != pipe_cost[b]) return pipe_cost[a] > pipe_cost[b];
            return a > b;
        });
        vector<char> keep(k, 0);
        for (int p : out) keep[p] = 1;
        for (int p : ord) {
            if (!keep[p]) continue;
            if (can_remove(p, occ, pc)) { apply_remove(p, occ, pc); keep[p] = 0; }
        }
        out.clear();
        for (int i = 0; i < k; i++) if (keep[i]) out.push_back(i);
        if (out.empty()) {
            int best = 0;
            uint64_t bc = pipe_cost[0];
            for (int p = 1; p < k; p++) if (pipe_cost[p] < bc) { bc = pipe_cost[p]; best = p; }
            out.push_back(best);
        }
    };

    auto greedy_variants = [&](vector<int>& best_out) {
        static const uint32_t BIAS[3] = {1024, 1152, 1280};
        vector<int> cand, tmp;
        unsigned long long bestC = ULLONG_MAX;
        for (uint32_t b : BIAS) {
            tmp.clear();
            greedy_cover_param(tmp, b);
            unsigned long long c = selection_cost(tmp);
            if (c < bestC) { bestC = c; cand = tmp; }
        }
        best_out.swap(cand);
    };

    auto run_ord_remove = [&](const vector<int>& ord, vector<int>& keep_out) {
        vector<uint32_t> cur_flow_occ = flow_occ_full, cur_pair_cnt(pair_counts.size(), 0);
        for (uint32_t id = 0; id < pair_counts.size(); ++id) cur_pair_cnt[id] = pair_counts[id];
        vector<char> keep(k, 1);
        for (int p : ord) {
            if (keep[p] && can_remove(p, cur_flow_occ, cur_pair_cnt)) { apply_remove(p, cur_flow_occ, cur_pair_cnt); keep[p] = 0; }
        }
        bool changed = true;
        while (changed) {
            changed = false;
            for (int p = 0; p < k; p++) if (keep[p] && can_remove(p, cur_flow_occ, cur_pair_cnt)) { apply_remove(p, cur_flow_occ, cur_pair_cnt); keep[p] = 0; changed = true; }
        }
        int cnt = 0; for (int i = 0; i < k; i++) if (keep[i]) cnt++;
        if (cnt == 0) {
            int best = 0; uint64_t bc = pipe_cost[0];
            for (int p = 1; p < k; p++) if (pipe_cost[p] < bc) { bc = pipe_cost[p]; best = p; }
            keep[best] = 1;
        }
        keep_out.clear();
        for (int i = 0; i < k; i++) if (keep[i]) keep_out.push_back(i);
    };

    vector<int> sub_idx_a; greedy_variants(sub_idx_a);

    vector<int> sub_idx_b1;
    {
        vector<int> base(k); iota(base.begin(), base.end(), 0);
        auto by_combo = [&](int a, int b) {
            __int128 ka = (__int128)pipe_cost[a] * 2 + (__int128)met_pair_sum[a] + (__int128)met_flow_sum[a];
            __int128 kb = (__int128)pipe_cost[b] * 2 + (__int128)met_pair_sum[b] + (__int128)met_flow_sum[b];
            if (ka != kb) return ka > kb;
            return a < b;
        };
        sort(base.begin(), base.end(), by_combo);
        run_ord_remove(base, sub_idx_b1);
        if (sub_idx_b1.empty()) {
            int best = 0; uint64_t bc = pipe_cost[0];
            for (int p = 1; p < k; p++) if (pipe_cost[p] < bc) { bc = pipe_cost[p]; best = p; }
            sub_idx_b1.push_back(best);
        }
    }

    vector<int> sub_idx_b2;
    {
        vector<int> base(k); iota(base.begin(), base.end(), 0);
        auto by_cost = [&](int a, int b) { if (pipe_cost[a] != pipe_cost[b]) return pipe_cost[a] > pipe_cost[b]; return a < b; };
        auto by_pair = [&](int a, int b) {
            if (met_pair_sum[a] != met_pair_sum[b]) return met_pair_sum[a] > met_pair_sum[b];
            if (pipe_cost[a] != pipe_cost[b]) return pipe_cost[a] > pipe_cost[b];
            return a < b;
        };
        auto by_flow = [&](int a, int b) {
            if (met_flow_sum[a] != met_flow_sum[b]) return met_flow_sum[a] > met_flow_sum[b];
            if (pipe_cost[a] != pipe_cost[b]) return pipe_cost[a] > pipe_cost[b];
            return a < b;
        };
        vector<int> oc = base, op = base, of = base;
        sort(oc.begin(), oc.end(), by_cost);
        sort(op.begin(), op.end(), by_pair);
        sort(of.begin(), of.end(), by_flow);
        vector<int> k1, k2, k3; run_ord_remove(oc, k1); run_ord_remove(op, k2); run_ord_remove(of, k3);
        unsigned long long c1 = selection_cost(k1), c2 = selection_cost(k2), c3 = selection_cost(k3);
        sub_idx_b2 = k1;
        if (c2 < selection_cost(sub_idx_b2)) sub_idx_b2 = k2;
        if (c3 < selection_cost(sub_idx_b2)) sub_idx_b2 = k3;
        if (sub_idx_b2.empty()) {
            int best = 0; uint64_t bc = pipe_cost[0];
            for (int p = 1; p < k; p++) if (pipe_cost[p] < bc) { bc = pipe_cost[p]; best = p; }
            sub_idx_b2.push_back(best);
        }
    }

    auto cost_sum = [&](const vector<int>& v) -> unsigned long long { return selection_cost(v); };
    vector<int> sub_idx = sub_idx_a;
    if (cost_sum(sub_idx_b1) < cost_sum(sub_idx)) sub_idx = sub_idx_b1;
    if (cost_sum(sub_idx_b2) < cost_sum(sub_idx)) sub_idx = sub_idx_b2;

    vector<int> new_idx; greedy_variants(new_idx);
    if (new_idx.empty()) {
        int best = 0; uint64_t bc = pipe_cost[0];
        for (int p = 1; p < k; p++) if (pipe_cost[p] < bc) { bc = pipe_cost[p]; best = p; }
        new_idx.push_back(best);
    }
    if (sub_idx.empty()) {
        int best = 0; uint64_t bc = pipe_cost[0];
        for (int p = 1; p < k; p++) if (pipe_cost[p] < bc) { bc = pipe_cost[p]; best = p; }
        sub_idx.push_back(best);
    }

    vector<char> keep(k, 0);
    for (int p : sub_idx) keep[p] = 1;
    vector<uint32_t> occ(m, 0), pc(pair_counts.size(), 0);
    for (int p = 0; p < k; p++) if (keep[p]) {
        for (uint32_t t = comp_off[p]; t < comp_off[p + 1]; ++t) occ[comp_flow_id[t]] += comp_flow_cnt[t];
        for (uint32_t t = pair_idx_off[p]; t < pair_idx_off[p + 1]; ++t) ++pc[pair_idx_data[t]];
    }

    vector<uint32_t> flow_use_off(m + 1, 0);
    for (int p = 0; p < k; p++) if (pipe_off[p] != pipe_off[p + 1]) for (uint32_t t = comp_off[p]; t < comp_off[p + 1]; ++t) ++flow_use_off[comp_flow_id[t]];
    for (int j = 1; j <= m; j++) flow_use_off[j] += flow_use_off[j - 1];
    vector<uint32_t> flow_use_data(flow_use_off[m], 0);
    {
        vector<uint32_t> cur = flow_use_off;
        for (int p = 0; p < k; p++) if (pipe_off[p] != pipe_off[p + 1]) {
            for (uint32_t t = comp_off[p]; t < comp_off[p + 1]; ++t) {
                uint32_t f = comp_flow_id[t];
                flow_use_data[--cur[f]] = p;
            }
        }
    }

    vector<uint32_t> def_flow(m, 0);
    vector<uint8_t> need_pair(pair_counts.size(), 0);
    vector<uint32_t> def_flow_idx; def_flow_idx.reserve(1024);
    vector<uint32_t> def_pair_idx; def_pair_idx.reserve(1024);
    vector<uint8_t> cand_mark(k, 0);
    const int MAX_TRIES = 32;
    const int CAND_LIMIT = 512;
    const int ADD_LIMIT = 2;
    vector<int> selected;
    selected.reserve(k);
    for (int p = 0; p < k; p++) if (keep[p]) selected.push_back(p);
    sort(selected.begin(), selected.end(), [&](int a, int b) {
        if (pipe_cost[a] != pipe_cost[b]) return pipe_cost[a] > pipe_cost[b];
        return a > b;
    });
    int tries = 0;
    for (int idx = 0; idx < (int)selected.size() && tries < MAX_TRIES; ++idx) {
        int p = selected[idx];
        if (!keep[p]) continue;
        bool removable = true;
        for (uint32_t t = pair_idx_off[p]; t < pair_idx_off[p + 1]; ++t) {
            uint32_t id = pair_idx_data[t];
            if (pc[id] <= 1u) { removable = false; break; }
        }
        if (removable) {
            for (uint32_t t = comp_off[p]; t < comp_off[p + 1] && removable; ++t) {
                uint32_t f = comp_flow_id[t];
                uint32_t c = comp_flow_cnt[t];
                if ((uint64_t)occ[f] < (uint64_t)w[f] + (uint64_t)c) removable = false;
            }
        }
        if (removable) {
            apply_remove(p, occ, pc);
            keep[p] = 0;
            ++tries;
            continue;
        }
        def_flow_idx.clear();
        def_pair_idx.clear();
        for (uint32_t t = comp_off[p]; t < comp_off[p + 1]; ++t) {
            uint32_t f = comp_flow_id[t];
            uint32_t c = comp_flow_cnt[t];
            uint32_t need = (occ[f] > c) ? (w[f] > (occ[f] - c) ? w[f] - (occ[f] - c) : 0u) : w[f];
            if (need) {
                if (def_flow[f] == 0) def_flow_idx.push_back(f);
                def_flow[f] = need;
            }
        }
        for (uint32_t t = pair_idx_off[p]; t < pair_idx_off[p + 1]; ++t) {
            uint32_t id = pair_idx_data[t];
            if (pc[id] <= 1u) {
                if (!need_pair[id]) def_pair_idx.push_back(id);
                need_pair[id] = 1;
            }
        }
        if (def_flow_idx.empty() && def_pair_idx.empty()) {
            apply_remove(p, occ, pc);
            keep[p] = 0;
            ++tries;
            continue;
        }
        vector<int> cand; cand.reserve(512);
        for (uint32_t f : def_flow_idx) {
            uint32_t st = flow_use_off[f], en = flow_use_off[f + 1];
            int added = 0;
            for (uint32_t pos = st; pos < en; ++pos) {
                int q = flow_use_data[pos];
                if (q == p || keep[q]) continue;
                if (!cand_mark[q]) { cand_mark[q] = 1; cand.push_back(q); if ((int)cand.size() >= CAND_LIMIT) break; }
                if (++added >= 256) break;
            }
            if ((int)cand.size() >= CAND_LIMIT) break;
        }
        for (int q : cand) cand_mark[q] = 0;
        if (cand.empty()) {
            for (uint32_t f : def_flow_idx) def_flow[f] = 0;
            for (uint32_t id : def_pair_idx) need_pair[id] = 0;
            continue;
        }
        auto eval_gain = [&](int q) -> pair<uint64_t, uint64_t> {
            uint64_t gflow = 0, gpair = 0;
            for (uint32_t t = comp_off[q]; t < comp_off[q + 1]; ++t) {
                uint32_t f = comp_flow_id[t];
                uint32_t c = comp_flow_cnt[t];
                uint32_t need = def_flow[f];
                if (need) gflow += min<uint64_t>(need, (uint64_t)c);
            }
            for (uint32_t t = pair_idx_off[q]; t < pair_idx_off[q + 1]; ++t) {
                uint32_t id = pair_idx_data[t];
                if (need_pair[id]) ++gpair;
            }
            return {gflow, gpair};
        };
        vector<int> added;
        uint64_t needPairs = def_pair_idx.size();
        uint64_t sumNeedFlow = 0;
        for (uint32_t f : def_flow_idx) sumNeedFlow += def_flow[f];
        int steps = 0;
        while ((needPairs > 0 || sumNeedFlow > 0) && steps < ADD_LIMIT) {
            int best = -1; __int128 bestK = -1;
            for (int q : cand) {
                auto gg = eval_gain(q);
                uint64_t gf = gg.first, gp = gg.second;
                if (gf == 0 && gp == 0) continue;
                uint64_t W = max<uint64_t>(1, (needPairs ? (sumNeedFlow / max<uint64_t>(1, needPairs)) : 1));
                __int128 gain = (__int128)gf + ((__int128)W * gp);
                __int128 key = (gain > 0) ? (((__int128)pipe_cost[q]) << 24) / gain : (__int128)((uint64_t)-1);
                if (best == -1 || key < bestK) { bestK = key; best = q; }
            }
            if (best == -1) break;
            added.push_back(best);
            for (uint32_t t = comp_off[best]; t < comp_off[best + 1]; ++t) {
                uint32_t f = comp_flow_id[t];
                uint32_t c = comp_flow_cnt[t];
                uint32_t take = min<uint32_t>(def_flow[f], c);
                if (take) { def_flow[f] -= take; sumNeedFlow -= take; }
            }
            for (uint32_t t = pair_idx_off[best]; t < pair_idx_off[best + 1]; ++t) {
                uint32_t id = pair_idx_data[t];
                if (need_pair[id]) { need_pair[id] = 0; if (needPairs) --needPairs; }
            }
            ++steps;
        }
        bool ok = (needPairs == 0 && sumNeedFlow == 0);
        for (uint32_t f : def_flow_idx) def_flow[f] = 0;
        for (uint32_t id : def_pair_idx) need_pair[id] = 0;
        if (!ok) continue;
        unsigned long long addCost = 0;
        for (int q : added) addCost += pipe_cost[q];
        if (addCost >= pipe_cost[p]) continue;
        apply_remove(p, occ, pc);
        keep[p] = 0;
        for (int q : added) if (!keep[q]) { apply_add(q, occ, pc); keep[q] = 1; }
        int checks = 0; bool changed = true;
        while (changed && checks < 32) {
            changed = false; ++checks;
            for (int r : selected) {
                if (!keep[r]) continue;
                if (can_remove(r, occ, pc)) { apply_remove(r, occ, pc); keep[r] = 0; changed = true; }
            }
        }
        ++tries;
    }

    vector<int> final_keep; final_keep.reserve(k);
    for (int p = 0; p < k; p++) if (keep[p]) final_keep.push_back(p);
    if (final_keep.empty()) {
        int best = 0; uint64_t bc = pipe_cost[0];
        for (int p = 1; p < k; p++) if (pipe_cost[p] < bc) { bc = pipe_cost[p]; best = p; }
        final_keep.push_back(best);
    }

    cout << (int)final_keep.size() << "\n";
    for (size_t i = 0; i < final_keep.size(); ++i) { if (i) cout << ' '; cout << final_keep[i]; }
    cout << "\n";

    cout << (int)new_idx.size() << "\n";
    for (int p : new_idx) {
        uint32_t st = pipe_off[p], en = pipe_off[p + 1];
        cout << (int)(en - st);
        for (uint32_t s = st; s < en; ++s) cout << ' ' << pipe_data[s];
        cout << "\n";
    }
    return 0;
}

