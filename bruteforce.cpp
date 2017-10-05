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

struct Car{
	vector<int> junctions;
};

struct Solution{
	vector<Car> cars;

	void print() {
		cout << cars.size() << endl;
		for(Car car : cars) {
			cout << car.junctions.size() << endl;
			for(int j : car.junctions)
				cout << j << endl;
		}
	}
};

struct State{
	int currentCar;
	int currentCarLocation;
	vector<bool> covered;
	Solution solution;

	bool operator<(const State &other) const {
		if(currentCar != other.currentCar)
			return currentCar < other.currentCar;
		if(currentCarLocation != other.currentCarLocation) {
			return currentCarLocation < other.currentCarLocation;
		}
		for(int i = 0; i < (int)covered.size(); i++) {
			if(covered[i] != other.streetCovers[i])
				return covered[i];
		}
		return 0;
	}
};

map<State, double> dp;
priority_queue<pair<double, State> > pq;

void addState(double timeRemaining, State s) {
	auto it = dp.find(s);
	if(it == dp.end() || timeRemaining > it->first) {
		dp[s] = timeRemaining;
		pq.push(s);
	}
}

void expandState(double timeRemaining, State s, TestCase data) {
	if(s.currentCar+1 < data.cars) {
		State newState = s;
		newState.currentCar++;
		newState.currentCarLocation = data.startIndex;
		newState.solution.cars[newState.car].junctions.push_back(data.startIndex);
		double newTimeRemaining = (data.cars - newState.car) * data.timeLimit;
		addState(newTimeRemaining, newState);
	}
	double carTimeRemaining = timeRemaining - (data.cars - newState.car + 1) * data.timeLimit;
	for(auto edge : data.outEdges[s.currentCarLocation]){
		if(edge.duration < carTimeRemaining) {
			int to = edge.other(s.currentCarLocation);
			State newState = s;
			newState.currentCarLocation = to;
			newState.covered[edge.index] = true;
			newState.solution.cars[newState.car].junctions.push_back(to);
			addState(timeRemaining - edge.duration, newState);
		}
	}
}

Solution bruteforce(TestCase data) {
	State init;
	init.currentCar = 0;
	init.covered.resize(data.streets.size());
	pq.push(make_pair(data.timeLimit * data.cars, init));
	while(!pq.empty()) {
		State s = pq.top();
		pq.pop();
		expandState(dp[s], s);
	}
}

int main(){
}
