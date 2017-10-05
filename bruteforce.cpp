#include <bits/stdc++.h>

using namespace std;

struct Junction{
	double lat, lng;
};

struct Street{
	int from, to;
	bool directed;
	int duration, length;
	int index;

	int other(int current) const {
		return to + from - current;
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
	data.outEdges.resize(N);

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
			if(covered[i] != other.covered[i])
				return covered[i];
		}
		return 0;
	}

	int score(const TestCase& data) {
		int ret = 0;
		for(int i = 0; i < (int)data.streets.size(); i++) {
			ret += covered[i] * data.streets[i].length;
		}
		return ret;
	}
};

map<State, int> dp;
priority_queue<pair<int, State> > pq;

void addState(int timeRemaining, const State& s) {
	auto it = dp.find(s);
	if(it == dp.end() || timeRemaining > it->second) {
		dp[s] = timeRemaining;
		pq.push(make_pair(timeRemaining, s));
	}
}

void expandState(int timeRemaining, const State& s, const TestCase& data) {
	if(s.currentCar+1 < data.cars) {
		State newState = s;
		newState.currentCar++;
		newState.currentCarLocation = data.startIndex;
		newState.solution.cars[newState.currentCar].junctions.push_back(data.startIndex);
		int newTimeRemaining = (data.cars - newState.currentCar) * data.timeLimit;
		addState(newTimeRemaining, newState);
	}
	int carTimeRemaining = timeRemaining - (data.cars - s.currentCar - 1) * data.timeLimit;
	//cerr << s.currentCarLocation << endl;
	for(auto& edge : data.outEdges[s.currentCarLocation]){
		if(edge.duration < carTimeRemaining) {
			int to = edge.other(s.currentCarLocation);
			State newState = s;
			newState.currentCarLocation = to;
			newState.covered[edge.index] = true;
			newState.solution.cars[newState.currentCar].junctions.push_back(to);
			addState(timeRemaining - edge.duration, newState);
		}
	}
}

Solution bruteforce(TestCase data) {
	State init;
	init.currentCar = 0;
	init.currentCarLocation = data.startIndex;
	init.covered.resize(data.streets.size());
	init.solution.cars.resize(data.cars);
	init.solution.cars[0].junctions.push_back(data.startIndex);
	addState(data.timeLimit * data.cars, init);
	State bestState = init;
	int bestSolutionScore = 0;

	int iteration = 0;
	while(!pq.empty()) {
		auto t = pq.top();
		State s = t.second;
		pq.pop();
		if(dp[s] != t.first)
			continue;
		if (iteration % 1000 == 0) {
			cerr << t.first << endl;
		}
		iteration++;
		int score = s.score(data);
		if(score >= bestSolutionScore) {
			bestSolutionScore = score;
			bestState = s;
		}
		expandState(t.first, s, data);
	}
	return bestState.solution;
}

int main(){
	auto testCase = parseTestCase();
	auto solution = bruteforce(testCase);
	solution.print();
}
