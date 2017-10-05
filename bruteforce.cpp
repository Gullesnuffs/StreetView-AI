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

	int score(const TestCase& data) const {
		int ret = 0;
		for(int i = 0; i < (int)data.streets.size(); i++) {
			ret += covered[i] * data.streets[i].length;
		}
		return ret;
	}

	int scoreUpperBound(const TestCase& data, int remainingTime) const {
		// More exact but slower upper bound
		/*vector<int> bestScore(remainingTime+1);
		bestScore[0] = score(data);
		for(int i = 0; i < (int)data.streets.size(); i++) {
			if(covered[i])
				continue;
			const Street& s = data.streets[i];
			for(int j = remainingTime; j >= s.duration; --j) {
				bestScore[j] = max(bestScore[j], bestScore[j-s.duration] + s.length);
			}
		}
		int ret = bestScore[0];
		for(int i = 1; i < remainingTime+1; i++)
			ret = max(ret, bestScore[i]);
		return ret;*/

		vector<pair<double, int> > remainingStreets;
		for(int i = 0; i < (int)data.streets.size(); i++) {
			if(covered[i])
				continue;
			double averageValue = ((double)data.streets[i].length) / data.streets[i].duration;
			remainingStreets.emplace_back(averageValue, i);
		}
		sort(remainingStreets.begin(), remainingStreets.end());
		int ub = score(data);
		int t = remainingTime;
		for(int i = (int)remainingStreets.size()-1; i >= 0; --i) {
			const Street& s = data.streets[remainingStreets[i].second];
			if(s.duration <= t){
				t -= s.duration;
				ub += s.length;
			}
			else {
				return ub + floor((s.length * t + 0.0001) / s.duration);
			}
		}
		return ub;
	}
};

struct QueueEntry {
	int remainingTime;
	int upperBound;
	State s;

	QueueEntry(int _remainingTime, int _upperBound, const State _s) : 
		remainingTime(_remainingTime),
		upperBound(_upperBound),
		s(_s) {
	}

	bool operator<(const QueueEntry& other) const {
		if(upperBound != other.upperBound)
			return upperBound < other.upperBound;
		if(remainingTime != other.remainingTime)
			return remainingTime < other.remainingTime;
		return s < other.s;
	}
};

map<State, int> dp;
priority_queue<QueueEntry> pq;

void addState(int remainingTime, const State& s, const TestCase& data) {
	auto it = dp.find(s);
	if(it == dp.end() || remainingTime > it->second) {
		dp[s] = remainingTime;
		pq.push(QueueEntry(remainingTime, s.scoreUpperBound(data, remainingTime), s));
	}
}

void expandState(int remainingTime, const State& s, const TestCase& data) {
	if(s.currentCar+1 < data.cars) {
		State newState = s;
		newState.currentCar++;
		newState.currentCarLocation = data.startIndex;
		newState.solution.cars[newState.currentCar].junctions.push_back(data.startIndex);
		int newTimeRemaining = (data.cars - newState.currentCar) * data.timeLimit;
		addState(newTimeRemaining, newState, data);
		if(newState.scoreUpperBound(data, newTimeRemaining) > s.scoreUpperBound(data, remainingTime)){
			assert(0);
		}
	}
	int carTimeRemaining = remainingTime - (data.cars - s.currentCar - 1) * data.timeLimit;
	//cerr << s.currentCarLocation << endl;
	for(auto& edge : data.outEdges[s.currentCarLocation]){
		if(edge.duration < carTimeRemaining) {
			int to = edge.other(s.currentCarLocation);
			State newState = s;
			newState.currentCarLocation = to;
			newState.covered[edge.index] = true;
			newState.solution.cars[newState.currentCar].junctions.push_back(to);
			if(newState.scoreUpperBound(data, remainingTime - edge.duration) > s.scoreUpperBound(data, remainingTime)){
				assert(0);
			}
			addState(remainingTime - edge.duration, newState, data);
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
	addState(data.timeLimit * data.cars, init, data);
	State bestState = init;
	int bestSolutionScore = 0;

	int iteration = 0;
	while(!pq.empty()) {
		auto cur = pq.top();
		pq.pop();
		if(dp[cur.s] != cur.remainingTime)
			continue;
		if (iteration % 1000 == 0) {
			cerr << bestSolutionScore << " - " << cur.upperBound << endl;
			//cerr << cur.remainingTime << endl;
		}
		iteration++;
		int score = cur.s.score(data);
		if(score >= bestSolutionScore) {
			bestSolutionScore = score;
			bestState = cur.s;
		}
		if(cur.upperBound <= bestSolutionScore) {
			break;
		}
		expandState(cur.remainingTime, cur.s, data);
	}
	cerr << "Best score: " << bestState.score(data) << endl;
	return bestState.solution;
}

int main(){
	auto testCase = parseTestCase();
	auto solution = bruteforce(testCase);
	solution.print();
}
