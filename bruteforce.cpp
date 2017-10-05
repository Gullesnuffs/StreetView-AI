#include <bits/stdc++.h>

using namespace std;

struct Junction{
	double lat, lng;
};

struct Street{
	int from, to;
	bool directed;
	double cost, length;
};

struct TestCase{
	vector<Junction> junctions;
	vector<Street> streets;
};

TestCase parseTestCase(){
	int N, M, T, C, S;
	TestCase data;
	data.junctions.resize(N);
	data.streets.resize(M);
	for(int i = 0; i < N; i++) {
		cin >> data.junction[i].lat >> data.junction[i].lng;
	}
	for(int i = 0; i < M; i++) {

	}
}

int main(){
}
