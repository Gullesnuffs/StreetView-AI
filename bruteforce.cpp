#include <bits/stdc++.h>

using namespace std;

struct Junction{
	double lat, lng;
};

struct Street{
	int from, to;
	bool directed;
	double duration, length;
	int index;

	int other(int current) {
		return from == current ? to : from;
	}
};

struct TestCase{
	vector<Junction> junctions;
	vector<vector<Street>> outEdges;
	vector<Street> streets;
	int cars;
	int startIndex;
	int timeLimit;
};

TestCase parseTestCase() {
	TestCase data;
	int N, M;
	cin >> N >> M >> data.timeLimit >> data.cars >> data.startIndex;
	data.junctions.resize(N);
	data.streets.resize(M);

	for (auto& junction : data.junctions) {
		cin >> junction.lat >> junction.lng;
	}

	for (int i = 0; i < M; i++) {
		auto& street = data.streets[i];
		street.index = i;
		int direction;
		cin >> street.from >> street.to >> direction >> street.duration >> street.length;
		street.directed = direction == 1;
		data.outEdges[street.from].push_back(street);
		if (!street.directed) data.outEdges[street.to].push_back(street);
	}
	return data;
}

int main(){
}
