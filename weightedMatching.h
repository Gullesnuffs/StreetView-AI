/**
 * See https://github.com/kth-competitive-programming/kactl
 * Author: Stanford
 * Date: Unknown
 * Source: Stanford Notebook
 * Description: Min cost bipartite matching. Negate costs for max cost.
 * Time: O(N^3)
 * Status: tested during ICPC 2015
 */

#include <bits/stdc++.h>
#define rep(i, a, b) for(int i = a; i < (b); ++i)
#define trav(a, x) for(auto& a : x)
#define all(x) x.begin(), x.end()
#define sz(x) (int)(x).size()

using namespace std;

typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<double> vd;

double MinCostMatching(const vector<vd>& cost, vi& L, vi& R);
