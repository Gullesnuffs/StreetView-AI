#include <bits/stdc++.h>
#include "weightedMatching.h"

using namespace std;

struct Junction{
	double lat, lng;
};

struct Street{
	int from, to;
	bool directed;
	int duration, length;
	int index;
	vector<int> innerJunctions;

	int other(int current) const {
		return to + from - current;
	}

	bool operator< (const Street &other) const {
		if(index != other.index) {
			return index < other.index;
		}
		return from < other.from;
	}
};

struct TestCase{
	vector<Junction> junctions;
	vector<vector<Street>> outEdges;
	vector<vector<Street>> inEdges;
	vector<Street> streets;
	vector<Street> orderedStreets;
	int cars;
	int startIndex;
	int timeLimit;

	Street getEdge(int from, int to) const {
		for(auto e : outEdges[from]) {
			if(e.other(from) == to)
				return e;
		}
		assert(0);
	}
};

void calculateOrderedStreets(TestCase& data) {
	vector<pair<double,int>> partial;
	for(int i = 0; i < (int)data.streets.size(); i++) {
		double averageValue = ((double)data.streets[i].length) / data.streets[i].duration;
		partial.emplace_back(averageValue, i);
	}
	sort(partial.begin(), partial.end());
	reverse(partial.begin(), partial.end());
	for (auto p : partial) {
		data.orderedStreets.push_back(data.streets[p.second]);
	}
}

void fillInAndOutEdges(TestCase& data) {
	for (auto& street : data.streets) {
		data.outEdges[street.from].push_back(street);
		data.inEdges[street.to].push_back(street);
		if (!street.directed) {
			data.outEdges[street.to].push_back(street);
			data.inEdges[street.from].push_back(street);
		}
	}
}

TestCase parseTestCase() {
	TestCase data;
	int N, M;
	cin >> N >> M >> data.timeLimit >> data.cars >> data.startIndex;
	data.junctions.resize(N);
	data.streets.resize(M);
	data.outEdges.resize(N);
	data.inEdges.resize(N);

	for (auto& junction : data.junctions) {
		cin >> junction.lat >> junction.lng;
	}

	set<pair<int, int> > edges;
	for (int i = 0; i < M; i++) {
		auto& street = data.streets[i];
		street.index = i;
		int direction;
		cin >> street.from >> street.to >> direction >> street.duration >> street.length;
		if(edges.count(make_pair(street.from, street.to))) {
			cerr << "Duplicate edges detected!" << endl;
			cerr << street.from << " " << street.to << " " << direction << " " << street.duration << " " << street.length << endl;
			assert(0);
		}
		edges.insert(make_pair(street.from, street.to));
		edges.insert(make_pair(street.to, street.from));
		street.directed = direction == 1;
	}
	fillInAndOutEdges(data);
	calculateOrderedStreets(data);
	for (int i = 0; i < N; i++) {
		vector<bool> vis(N);
		queue<int> q;
		q.push(i);
		vis[i] = 1;
		int visN = 1;
		while(!q.empty()) {
			int cur = q.front();
			q.pop();
			for(auto e : data.outEdges[cur]) {
				int to = e.other(cur);
				if(!vis[to]) {
					vis[to] = true;
					++visN;
					q.push(to);
				}
			}
		}
		if(visN < N) {
			cerr << "Not every node is reachable from node " << i << endl;
			assert(0);
		}
	}
	return data;
}

bool hasStreetBetween(const TestCase& data, int junction1, int junction2) {
	for (Street s : data.outEdges[junction1]) {
		if (s.other(junction1) == junction2) return true;
	}
	return false;
}

Street streetBetween(const TestCase& data, int junction1, int junction2) {
	for (Street s : data.outEdges[junction1]) {
		if (s.other(junction1) == junction2) return s;
	}
	cerr << junction1 << " " << junction2 << endl;
	assert(false);
}

struct Car{
	vector<int> junctions;
};

struct PartialSolution {
	int car;
	int junction;
	PartialSolution* previous;

	PartialSolution(PartialSolution* previous, int car, int junction) : previous(previous), car(car), junction(junction) {
	 }
};

struct Solution{
	vector<Car> cars;

	Solution() {}

	Solution(TestCase& data, PartialSolution* partial) {
		cars.resize(data.cars);
		while(partial != nullptr) {
			cars[partial->car].junctions.push_back(partial->junction);
			partial = partial->previous;
		}
		for(auto& car : cars) {
			reverse(car.junctions.begin(), car.junctions.end());
		}
	}

	void print() {
		cout << cars.size() << endl;
		for(Car car : cars) {
			cout << car.junctions.size() << endl;
			for(int j : car.junctions)
				cout << j << endl;
		}
	}

	int score(const TestCase& data) const {
		vector<bool> covered(data.streets.size());
		for(const Car& car : cars) {
			for (int i = 0; i < car.junctions.size()-1; i++) {
				int a = car.junctions[i];
				int b = car.junctions[i+1];
				for (auto& edge : data.streets) {
					if (edge.other(a) == b && edge.other(b) == a && !covered[edge.index]) {
						covered[edge.index] = true;
						break;
					}
				}
				//covered[streetBetween(data, car.junctions[i], car.junctions[i+1]).index] = true;
			}
		}

		int ret = 0;
		for(int i = 0; i < (int)data.streets.size(); i++) {
			ret += covered[i] * data.streets[i].length;
		}
		return ret;
	}
};

struct BruteForceState{
	int currentCar;
	int currentCarLocation;
	vector<bool> covered;
	PartialSolution* solution;

	bool operator<(const BruteForceState &other) const {
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
		int ub = score(data);

		int t = remainingTime;
		for (auto& s : data.orderedStreets) {
			if(covered[s.index]) continue;

			if(s.duration <= t){
				t -= s.duration;
				ub += s.length;
			}
			else {
				return ub + (int) floor((s.length * t + 0.0001) / s.duration);
			}
		}
		return ub;
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

		int ub = score(data);

		int t = remainingTime;
		for (auto& s : data.orderedStreets) {
			if(covered[s.index]) continue;

			if(s.duration <= t){
				t -= s.duration;
				ub += s.length;
			}
			else {
				return ub + (int) floor((s.length * t + 0.0001) / s.duration);
			}
		}
		return ub;
	}
};

struct QueueEntry {
	int remainingTime;
	int upperBound;
	BruteForceState s;

	QueueEntry(int _remainingTime, int _upperBound, const BruteForceState _s) :
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

TestCase compressTestCase(TestCase& data, Solution& solution) {
	vector<int> split(data.junctions.size());


	vector<bool> covered(data.streets.size());
	int maxAdditional = 1000;
	for(const Car& car : solution.cars) {
		for (int i = 0; i < car.junctions.size()-1; i++) {
			covered[streetBetween(data, car.junctions[i], car.junctions[i+1]).index] = true;
		}
	}

	for (int i = 0; i < covered.size(); i++) {
		if (!covered[i]) {
			maxAdditional--;
			if (maxAdditional <= 0) {
				covered[i] = true;
			}
		}
	}

	for(auto& edge : data.streets) {
		if (!covered[edge.index]) {
			split[edge.from] = 10;
			split[edge.to] = 10;
		}
	}

	int splitEvery = 10000;
	split[data.startIndex] = 10;
	for(auto& car : solution.cars) {
		int dist = 0;
		for (int i = 0; i < (int)car.junctions.size(); i++) {
			int index = car.junctions[i];
			if ((split[index] > 0 && dist >= splitEvery/2) || dist >= splitEvery) {
				split[index]++;
				dist = 0;
			}
			dist++;
		}
		split[car.junctions.back()] = 10;
	}

	// Remove splits only visited once (optional, but usually good)
	for (int& x : split) {
		if (x <= 1) x = 0;
	}

	int numSplits = 0;
	for (int x : split) if(x > 0) numSplits++;

	TestCase newData;
	newData.junctions = vector<Junction>(numSplits);
	newData.outEdges = vector<vector<Street>>(numSplits);
	newData.inEdges = vector<vector<Street>>(numSplits);
	newData.cars = data.cars;
	newData.timeLimit = data.timeLimit;
	vector<bool> taken(data.streets.size());

	vector<int> split2nodeIndex = vector<int>(split.size(), -1);
	int nodeIndex = 0;
	for (int i = 0, j = 0; i < (int)split2nodeIndex.size(); i++) {
		if (split[i] > 0) {
			split2nodeIndex[i] = j;
			newData.junctions[j] = data.junctions[i];
			j++;
		}
	}

	newData.startIndex = split2nodeIndex[data.startIndex];

	for(auto street : data.streets) {
		if (!covered[street.index]) {
			cout << "DEBUG:LINE " << 0 << " " << street.from << " " << street.to << endl;
			street.innerJunctions.push_back(street.from);
			street.innerJunctions.push_back(street.to);
			street.from = split2nodeIndex[street.from];
			street.to = split2nodeIndex[street.to];
			assert(street.from != -1);
			assert(street.to != -1);
			street.index = newData.streets.size();
			newData.streets.push_back(street);
		}
	}

	for(auto& car : solution.cars) {
		int dist = 0;
		int i = 0;
		while (i < (int)car.junctions.size() - 1) {
			Street combinedStreet;
			int origFrom = car.junctions[i];
			combinedStreet.from = split2nodeIndex[car.junctions[i]];

			vector<Street> span;
			while(!split[car.junctions[i+1]]) {
				combinedStreet.innerJunctions.push_back(car.junctions[i]);
				span.push_back(streetBetween(data, car.junctions[i], car.junctions[i+1]));
				i++;
			}
			combinedStreet.innerJunctions.push_back(car.junctions[i]);
			span.push_back(streetBetween(data, car.junctions[i], car.junctions[i+1]));
			i++;
			combinedStreet.innerJunctions.push_back(car.junctions[i]);

			assert(combinedStreet.innerJunctions.front() == origFrom);
			assert(combinedStreet.innerJunctions.back() == car.junctions[i]);

			combinedStreet.duration = 0;
			combinedStreet.length = 0;
			combinedStreet.directed = false;
			combinedStreet.to = split2nodeIndex[car.junctions[i]];
			combinedStreet.index = newData.streets.size();
			for (int k = 0; k < combinedStreet.innerJunctions.size() - 1; k++) {
				streetBetween(data, combinedStreet.innerJunctions[k], combinedStreet.innerJunctions[k+1]);
			}

			for (Street& s : span) {
				if (!taken[s.index]) {
					combinedStreet.length += s.length;
				}
				taken[s.index] = true;
				combinedStreet.duration += s.duration;
				combinedStreet.directed |= s.directed;
			}
			cout << "DEBUG:LINE " << 0 << " " << origFrom << " " << car.junctions[i] << endl;

			newData.streets.push_back(combinedStreet);
		}
	}

	fillInAndOutEdges(newData);
	calculateOrderedStreets(newData);
	return newData;
}

Solution expandSolution (TestCase& originalTestCase, TestCase& expandedTestCase, Solution& solution) {
	vector<bool> covered(expandedTestCase.streets.size());
	Solution newSolution;
	newSolution.cars.resize(solution.cars.size());
	for(int c = 0; c < (int)solution.cars.size(); c++) {
		const Car& car = solution.cars[c];
		auto& newJunctions = newSolution.cars[c].junctions;

		for (int i = 0; i < car.junctions.size()-1; i++) {
			int a = car.junctions[i];
			int b = car.junctions[i+1];
			Street bestEdge;
			bool found = false;
			for (auto& edge : expandedTestCase.orderedStreets) {
				if (edge.other(a) == b && edge.other(b) == a) {
					bestEdge = edge;
					found = true;
					if (!covered[edge.index]) break;
				}
			}

			assert(found);
			covered[bestEdge.index] = true;
			if (a == bestEdge.from) {
				for (auto it = bestEdge.innerJunctions.begin(); it != bestEdge.innerJunctions.end(); it++) {
					if (newJunctions.size() > 0) {
						assert(hasStreetBetween(originalTestCase, newJunctions.back(), *it));
					}
					newJunctions.push_back(*it);
				}
			} else {
				for (auto it = bestEdge.innerJunctions.rbegin(); it != bestEdge.innerJunctions.rend(); it++) {
					if (newJunctions.size() > 0) {
						assert(hasStreetBetween(originalTestCase, newJunctions.back(), *it));
					}
					newJunctions.push_back(*it);
				}
			}

			// Add the last junction in the path
			if (i != (int)car.junctions.size()-2) {
				newJunctions.pop_back();
			}
		}
	}

	for (auto& c : newSolution.cars) {
		for (auto j : c.junctions) {
			assert(j >= 0 && j < originalTestCase.junctions.size());
		}
		for (int i = 0; i < c.junctions.size()-1; i++) {
			assert(hasStreetBetween(originalTestCase, c.junctions[i], c.junctions[i+1]));
		}
	}
	return newSolution;
}

map<BruteForceState, int> dp;
priority_queue<QueueEntry> pq;

void addState(int remainingTime, const BruteForceState& s, const TestCase& data) {
	auto it = dp.find(s);
	if(it == dp.end() || remainingTime > it->second) {
		dp[s] = remainingTime;
		pq.push(QueueEntry(remainingTime, s.scoreUpperBound(data, remainingTime), s));
	}
}

void expandState(int remainingTime, const BruteForceState& s, const TestCase& data) {
	if(s.currentCar+1 < data.cars) {
		BruteForceState newState = s;
		newState.currentCar++;
		newState.currentCarLocation = data.startIndex;
		newState.solution = new PartialSolution(s.solution, newState.currentCar, data.startIndex);
		//.cars[newState.currentCar].junctions.push_back(data.startIndex);
		int newTimeRemaining = (data.cars - newState.currentCar) * data.timeLimit;
		addState(newTimeRemaining, newState, data);
		#if DEBUG
		assert(newState.scoreUpperBound(data, newTimeRemaining) <= s.scoreUpperBound(data, remainingTime));
		#endif
	}
	int carTimeRemaining = remainingTime - (data.cars - s.currentCar - 1) * data.timeLimit;
	//cerr << s.currentCarLocation << endl;
	for(auto& edge : data.outEdges[s.currentCarLocation]){
		if(edge.duration < carTimeRemaining) {
			int to = edge.other(s.currentCarLocation);
			BruteForceState newState = s;
			newState.currentCarLocation = to;
			newState.covered[edge.index] = true;
			//newState.solution.cars[newState.currentCar].junctions.push_back(to);
			newState.solution = new PartialSolution(s.solution, newState.currentCar, to);
			#if DEBUG
			assert(newState.scoreUpperBound(data, remainingTime - edge.duration) <= s.scoreUpperBound(data, remainingTime));
			#endif
			addState(remainingTime - edge.duration, newState, data);
		}
	}
}

clock_t w0;
clock_t w1;
clock_t w2;
clock_t w3;
clock_t w4;

// Exact exponential-time solution based on A*
Solution bruteforce(TestCase data) {
	BruteForceState init;
	init.currentCar = 0;
	init.currentCarLocation = data.startIndex;
	init.covered.resize(data.streets.size());
	init.solution = new PartialSolution(nullptr, 0, data.startIndex);
	//init.solution.cars.resize(data.cars);
	//init.solution.cars[0].junctions.push_back(data.startIndex);
	addState(data.timeLimit * data.cars, init, data);
	BruteForceState bestState = init;
	int bestSolutionScore = 0;

	int iteration = 0;
	while(!pq.empty()) {
		auto t0 = clock();
		auto cur = pq.top();
		pq.pop();
		auto t1 = clock();
		w0 += t1 - t0;
		if(dp[cur.s] != cur.remainingTime)
			continue;
		if (iteration % 1000 == 0) {
			cerr << bestSolutionScore << " - " << cur.upperBound << endl;
			cerr << w0 << " " << w1 << " " << w2 << " " << w3 << " " << w4 << endl;
			//cerr << cur.remainingTime << endl;
		}
		iteration++;
		int score = cur.s.score(data);
		auto t2 = clock();
		w1 += t2 - t1;
		if(score >= bestSolutionScore) {
			bestSolutionScore = score;
			bestState = cur.s;
		}
		auto t3 = clock();
		w2 += t3 - t2;
		if(cur.upperBound <= bestSolutionScore) {
			break;
		}
		expandState(cur.remainingTime, cur.s, data);
		auto t4 = clock();
		w3 += t4 - t3;
	}
	cerr << "Best score: " << bestState.score(data) << endl;
	return Solution(data, bestState.solution);
}

int countOut(const TestCase& data, const State& s, map<Street, int>& source, map<Street, int>& destination, int node) {
	int ret = 0;
	for(auto e : data.outEdges[node]) {
		if(s.covered[e.index])
			continue;
		if(source.count(e) && source[e] != node)
			continue;
		ret++;
	}
	return ret;
}

int countIn(const TestCase& data, const State& s, map<Street, int>& source, map<Street, int>& destination, int node) {
	int ret = 0;
	for(auto e : data.inEdges[node]) {
		if(s.covered[e.index])
			continue;
		if(destination.count(e) && destination[e] != node)
			continue;
		ret++;
	}
	return ret;
}

struct Path {
	int start;
	vector<Street> edges;

	Path() {
	}

	Path(int _start, vector<Street> _edges) {
		start = _start;
		edges = _edges;
	}

	int duration() {
		int ret = 0;
		for(auto e : edges) {
			ret += e.duration;
		}
		return ret;
	}
};

struct EulerGraph {
	map<Street, int> source;
	map<Street, int> destination;
	map<Street, int> extraEdges;
	map<int, map<Street, int> > extraOutEdges;
	map<int, map<Street, int> > extraInEdges;
	TestCase data;
	State state;

	EulerGraph(const TestCase& _data, const State& _state) {
		data = _data;
		state = _state;
	}

	int outDegree(int node) {
		int ret = 0;
		for(auto e : data.outEdges[node]) {
			if(state.covered[e.index])
				continue;
			if(e.directed || (source.count(e) && source[e] == node)) {
				++ret;
			}
		}
		for(auto it : extraOutEdges[node]) {
			ret += it.second;
		}
		return ret;
	}

	int inDegree(int node) {
		int ret = 0;
		for(auto e : data.inEdges[node]) {
			if(state.covered[e.index])
				continue;
			if(e.directed || (destination.count(e) && destination[e] == node)) {
				++ret;
			}
		}
		for(auto it : extraInEdges[node]) {
			ret += it.second;
		}
		return ret;
	}

	int relativeDegree(int node) {
		return outDegree(node) - inDegree(node);
	}

	int totalCost() {
		int ret = 0;
		for(auto e : extraEdges) {
			ret += e.first.duration * e.second;
		}
		return ret;
	}

	bool isOk() {
		for(int node = 0; node < (int)data.outEdges.size(); ++node) {
			if(relativeDegree(node))
				return false;
		}
		return true;
	}

	void addExtraEdge(Street s, int src, int dst) {
		if(src != s.from) {
			swap(s.from, s.to);
		}
		assert(s.from == src && s.to == dst);
		extraEdges[s]++;
		extraOutEdges[src][s]++;
		extraInEdges[dst][s]++;
		//cerr << "Added edge " << src << " " << dst << endl;
	}

	vector<Street> getCycle() {
		int curNode = data.startIndex;
		vector<Street> ret;
		vector<vector<Street> > edges(data.outEdges.size());
		set<int> unvisited;
		unvisited.insert(data.startIndex);
		for(auto e : data.streets) {
			if(state.covered[e.index]) {
				continue;
			}
			unvisited.insert(e.from);
			unvisited.insert(e.to);
			if(e.directed) {
				edges[e.from].push_back(e);
			}
			else {
				edges[source[e]].push_back(e);
			}
		}
		for(auto it : extraEdges) {
			for(int i = 0; i < it.second; i++) {
				edges[it.first.from].push_back(it.first);
			}
		}
		vector<vector<Street>::iterator> its;
		for(int i = 0; i < (int)data.outEdges.size(); i++) {
			its.push_back(edges[i].begin());
		}
		int startNode = data.startIndex;
		while(true) {
			priority_queue<pair<int, int> > q;
			map<int, int> minDis;
			map<int, int> parent;
			map<int, Street> parentEdge;
			q.push(make_pair(0, startNode));
			minDis[startNode] = 0;
			parent[startNode] = -1;
			bool hasUnvisited = false;
			while(!q.empty()){
				auto cur = q.top();
				q.pop();
				int d = -cur.first;
				int node = cur.second;
				if(minDis[node] < d)
					continue;
				if(unvisited.count(node)) {
					vector<Street> v;
					while(node != startNode){
						Street e = parentEdge[node];
						v.push_back(e);
						node = parent[node];
					}
					reverse(v.begin(), v.end());
					for(auto e : v)
						ret.push_back(e);
					startNode = cur.second;
					hasUnvisited = true;
					break;
				}
				for(auto e : data.outEdges[node]) {
					int to = e.other(node);
					int newDis = d + e.duration;
					auto it = minDis.find(to);
					if(it == minDis.end() || newDis < it->second) {
						q.push(make_pair(-newDis, to));
						minDis[to] = newDis;
						parent[to] = node;
						parentEdge[to] = e;
					}
				}
			}
			if(!hasUnvisited)
				break;
			vector<int> s = {startNode};
			vector<Street> E;
			vector<Street> v;
			while(!s.empty()) {
				int x = s.back();
				unvisited.erase(x);
				auto& it = its[x], end = edges[x].end();
				if(it == end) { 
					s.pop_back(); 
					if(!E.empty()) {
						v.push_back(E.back()); 
						E.pop_back();
					}
				}
				else { 
					s.push_back(it->other(x)); 
					E.push_back(*it);
					++its[x];
				}
			}
			reverse(all(v));
			for(auto e : v)
				ret.push_back(e);
		}
		return ret;
	}

	void greedy() {
		cerr << "Using greedy " << endl;
		for(auto e : data.streets) {
			if(e.directed)
				continue;
			if(relativeDegree(e.from) > relativeDegree(e.to)) {
				source[e] = e.to;
				destination[e] = e.from;
			}
			else {
				source[e] = e.from;
				destination[e] = e.to;
			}
		}
		for(int curNode = -1; curNode < (int)data.outEdges.size(); curNode++) {
			bool needOutEdge = false;
			bool initialCheck = false;
			if(curNode == -1) {
				curNode = data.startIndex;
				initialCheck = true;
				if(outDegree(curNode) == 0) {
					needOutEdge = true;
				}
			}
			while(relativeDegree(curNode) < 0 || needOutEdge) {
				needOutEdge = false;
				priority_queue<pair<int, int> > q;
				map<int, int> minDis;
				map<int, int> parent;
				map<int, Street> parentEdge;
				q.push(make_pair(0, curNode));
				minDis[curNode] = 0;
				parent[curNode] = -1;
				while(!q.empty()){
					auto cur = q.top();
					q.pop();
					int d = -cur.first;
					int node = cur.second;
					if(minDis[node] < d)
						continue;
					if(relativeDegree(node) > 0) {
						while(node != curNode){
							Street e = parentEdge[node];
							addExtraEdge(e, parent[node], node);
							node = parent[node];
						}
						break;
					}
					for(auto e : data.outEdges[node]) {
						int to = e.other(node);
						int newDis = d + e.duration;
						auto it = minDis.find(to);
						if(it == minDis.end() || newDis < it->second) {
							q.push(make_pair(-newDis, to));
							minDis[to] = newDis;
							parent[to] = node;
							parentEdge[to] = e;
						}
					}
				}
			}
			if(initialCheck)
				curNode = -1;
		}
	}

	map<int, Path> computeDistances(int from, set<int> targets) {
		priority_queue<pair<int, int> > q;
		map<int, int> minDis;
		map<int, int> parent;
		map<int, Street> parentEdge;
		map<int, Path> paths;
		q.push(make_pair(0, from));
		minDis[from] = 0;
		parent[from] = -1;
		set<int> remainingTargets = targets;
		while(!q.empty()){
			auto cur = q.top();
			q.pop();
			int d = -cur.first;
			int node = cur.second;
			if(minDis[node] < d)
				continue;
			if(remainingTargets.count(node)) {
				vector<Street> path;
				int tmpNode = node;
				while(tmpNode != from){
					path.push_back(parentEdge[tmpNode]);
					tmpNode = parent[tmpNode];
				}
				reverse(path.begin(), path.end());
				paths[node] = Path(from, path);
			}
			for(auto e : data.outEdges[node]) {
				int to = e.other(node);
				int newDis = d + e.duration;
				auto it = minDis.find(to);
				if(it == minDis.end() || newDis < it->second) {
					q.push(make_pair(-newDis, to));
					minDis[to] = newDis;
					parent[to] = node;
					if(to == e.from) {
						swap(e.from, e.to);
					}
					parentEdge[to] = e;
				}
			}
			remainingTargets.erase(node);
			if(!remainingTargets.size()) {
				break;
			}
		}
		return paths;
	}

	int bestCost;
	map<Street, int> bestSource;
	map<Street, int> bestDestination;
	vector<Street> bestExtra;

	void optimize(set<Street> undirected, set<int> nodes) {
		if(undirected.size()) {
			Street s = *undirected.begin();
			assert(!s.directed);
			undirected.erase(s);
			source[s] = s.from;
			destination[s] = s.to;
			optimize(undirected, nodes);
			source[s] = s.to;
			destination[s] = s.from;
			optimize(undirected, nodes);
		}
		else {
			map<int, map<int, Path> > paths;
			vector<int> indices;
			vector<int> sources;
			vector<int> sinks;
			for(int node : nodes) {
				paths[node] = computeDistances(node, nodes);
				int r = relativeDegree(node);
				if(node == data.startIndex && outDegree(node) == 0 && inDegree(node) == 0) {
					sinks.push_back(node);
					indices.push_back(node);
					sources.push_back(node);
					indices.push_back(node);
				}
				if(r > 0) {
					while(r--) {
						sinks.push_back(node);
						indices.push_back(node);
					}
				}
				else {
					while(r++) {
						sources.push_back(node);
						indices.push_back(node);
					}
				}
			}
			assert(sources.size() == sinks.size());
			vector<vector<double> > costs(sources.size());
			for(int i = 0; i < (int)sources.size(); i++) {
				for(int j = 0; j < (int)sinks.size(); j++) {
					costs[i].push_back(paths[sources[i]][sinks[j]].duration());
				}
			}
			vector<int> L;
			vector<int> R;
			int C = (int)MinCostMatching(costs, L, R);
			vector<Street> extraEdges;
			for(int i = 0; i < (int)L.size(); i++) {
				for(auto e : paths[sources[i]][sinks[L[i]]].edges) {
					extraEdges.push_back(e);
				}
			}
			if(C < bestCost) {
				bestCost = C;
				bestSource = source;
				bestDestination = destination;
				bestExtra = extraEdges;
			}
		}
	}

	void removeAllExtra(Street e) {
		for(int i = 0; i < 2; i++) {
			swap(e.from, e.to);
			if(extraEdges.count(e) && extraEdges[e]) {
				//cerr << "Removed edge " << e.from << " " << e.to << " " << extraEdges[e] << " times" << endl;
				extraEdges.erase(e);
				extraOutEdges[e.from].erase(e);
				extraOutEdges[e.to].erase(e);
				extraInEdges[e.from].erase(e);
				extraInEdges[e.to].erase(e);
			}
		}
	}

	void optimize() {
		for(int centerNode = 0; centerNode < (int)data.outEdges.size(); ++centerNode) {
			priority_queue<pair<int, int> > q;
			map<int, int> minDis;
			set<int> nodes;
			set<Street> undirected;
			minDis[centerNode] = 0;
			q.push(make_pair(0, centerNode));
			while(!q.empty() && nodes.size() < 8) {
				auto cur = q.top();
				q.pop();
				int d = -cur.first;
				int node = cur.second;
				if(minDis[node] < d)
					continue;
				nodes.insert(node);
				for(auto e : data.outEdges[node]) {
					int to = e.other(node);
					int newDis = d + e.duration;
					auto it = minDis.find(to);
					if(it == minDis.end() || newDis < it->second) {
						q.push(make_pair(-newDis, to));
						minDis[to] = newDis;
					}
					if(!e.directed && undirected.size() < 5 && nodes.count(node) && nodes.count(to)) {
						undirected.insert(e);
					}
				}
				for(auto e : data.inEdges[node]) {
					int to = e.other(node);
					int newDis = d + e.duration;
					auto it = minDis.find(to);
					if(it == minDis.end() || newDis < it->second) {
						q.push(make_pair(-newDis, to));
						minDis[to] = newDis;
					}
				}
			}
			int previousTotal = totalCost();
			assert(isOk());
			/*cerr << "Optimizing over nodes";
			for(auto node : nodes)
				cerr << " " << node;
			cerr << endl;*/
			/*for(auto e : data.streets) {
				if(nodes.count(e.from) && nodes.count(e.to)) {
					if(!e.directed) {
						undirected.insert(e);
					}
					removeAllExtra(e);
				}
			}*/
			for(auto e : undirected) {
				removeAllExtra(e);
			}
			int intermediateCost = totalCost();
			bestCost = 1000000000;
			optimize(undirected, nodes);
			source = bestSource;
			destination = bestDestination;
			//cerr << "Adding extra" << endl;
			for(auto e : bestExtra) {
				addExtraEdge(e, e.from, e.to);
			}
			//cerr << "Finished adding extra" << endl;
			int newTotal = totalCost();
			assert(isOk());
			if(newTotal > previousTotal) {
				cerr << "ERROR! New total is " << newTotal << endl;
				cerr << "Old was " << previousTotal << endl;
				cerr << "Intermediate was " << intermediateCost << endl;
				assert(0);
			}
			if(newTotal < previousTotal) {
				cerr << newTotal << endl;
			}
		}
	}
};

int totalScore;

State extendSolution(const TestCase& data, State s) {
	map<Street, int> source;
	map<Street, int> destination;
	vector<int> outDegree(data.junctions.size());
	vector<int> inDegree(data.junctions.size());
	vector<int> totDegree(data.junctions.size());
	for(auto e : data.streets) {
		if(s.covered[e.index])
			continue;
		totDegree[e.from]++;
		totDegree[e.to]++;
	}
	for(auto e : data.streets) {
		if(s.covered[e.index])
			continue;
		if(e.directed)
			continue;
		int rel1 = countOut(data, s, source, destination, e.from) - 
		countIn(data, s, source, destination, e.from);
		int rel2 = countOut(data, s, source, destination, e.to) - 
		countIn(data, s, source, destination, e.to);
		if(rel1 > rel2) {
			source[e] = e.to;
			destination[e] = e.from;
		}
		else{
			source[e] = e.from;
			destination[e] = e.to;
		}
	}
	for(int i = 0; i < (int)data.junctions.size(); i++) {
		int totOut = 0;
		int totIn = 0;
		for(auto e : data.outEdges[i]){
			if(s.covered[e.index])
				continue;
			++totOut;
		}
		for(auto e : data.inEdges[i]){
			if(s.covered[e.index])
				continue;
			++totIn;
		}
		outDegree[i] = totOut;
		inDegree[i] = totIn;
		assert(totDegree[i] >= outDegree[i]);
	}
	int remainingTime = data.timeLimit;
	Car c = s.solution.cars[s.currentCar];
	for(int i = 0; i < (int)c.junctions.size()-1; ++i) {
		remainingTime -= data.getEdge(c.junctions[i], c.junctions[i+1]).duration;
	}
	assert(remainingTime >= 0);
	while(true) {
		double bestEdgeValue = -1;
		Street bestEdge;
		for(auto e : data.outEdges[s.currentCarLocation]) {
			if(s.covered[e.index] || e.duration > remainingTime)
				continue;
			double value = (rand()%1000)+10000;
			//value *= ((double)e.length) / e.duration + 1;
			//value /= (1 + totDegree[e.other(s.currentCarLocation)]);
			if(e.directed)
				value *= 5;
			else {
				if(source[e] == s.currentCarLocation) {
					value *= 1;
				}
				else {
					value /= 1;
				}
			}
			/*if(!e.directed && outDegree[e.other(s.currentCarLocation)] < inDegree[e.other(s.currentCarLocation)]) {
				value /= 1.3;
			}*/
			if(value > bestEdgeValue) {
				bestEdgeValue = value;
				bestEdge = e;
			}
		}
		if(bestEdgeValue >= 0) {
			int to = bestEdge.other(s.currentCarLocation);
			assert(totDegree[s.currentCarLocation] >= outDegree[s.currentCarLocation]);
			assert(totDegree[to] >= outDegree[to]);
			--totDegree[s.currentCarLocation];
			--totDegree[to];
			--outDegree[s.currentCarLocation];
			--inDegree[to];
			if(!bestEdge.directed) {
				--inDegree[s.currentCarLocation];
				--outDegree[to];
			}
			assert(totDegree[s.currentCarLocation] >= outDegree[s.currentCarLocation]);
			assert(totDegree[to] >= outDegree[to]);
			s.currentCarLocation = to;
			s.covered[bestEdge.index] = true;
			s.solution.cars[s.currentCar].junctions.push_back(to);
			remainingTime -= bestEdge.duration;
			//cerr << "Good " << e.duration << endl;
			continue;
		}
		priority_queue<pair<int, int> > q;
		map<int, int> minDis;
		map<int, int> parent;
		map<int, Street> parentEdge;
		q.push(make_pair(0, s.currentCarLocation));
		minDis[s.currentCarLocation] = 0;
		parent[s.currentCarLocation] = -1;
		bool failed = true;
		while(!q.empty()){
			auto cur = q.top();
			q.pop();
			int d = -cur.first;
			int node = cur.second;
			if(minDis[node] < d)
				continue;
			assert(totDegree[node] >= outDegree[node]);
			if(outDegree[node] && /* && outDegree[node] >= inDegree[node] && */totDegree[node]%2 >= 1 && node != s.currentCarLocation) {
				vector<int> path;
				while(node != s.currentCarLocation){
					path.push_back(node);
					node = parent[node];
				}
				reverse(path.begin(), path.end());
				//cerr << "Path length " << path.size() << endl;
				for(int i = 0; i < (int)path.size(); i++) {
					int to = path[i];
					const Street& edge = parentEdge[to];
					if(edge.duration > remainingTime)
						break;
					failed = false;
					if(!s.covered[edge.index]) {
						assert(s.currentCarLocation == edge.other(to));
						assert(totDegree[s.currentCarLocation] >= outDegree[s.currentCarLocation]);
						assert(totDegree[to] >= outDegree[to]);
						--totDegree[s.currentCarLocation];
						--totDegree[to];
						--outDegree[s.currentCarLocation];
						--inDegree[to];
						if(!edge.directed) {
							--inDegree[s.currentCarLocation];
							--outDegree[to];
						}
						assert(totDegree[s.currentCarLocation] >= outDegree[s.currentCarLocation]);
						assert(totDegree[to] >= outDegree[to]);
					}
					s.currentCarLocation = to;
					s.covered[edge.index] = true;
					s.solution.cars[s.currentCar].junctions.push_back(to);
					remainingTime -= edge.duration;
					//cerr << "Bad  " << edge.duration << endl;
				}
				break;
			}
			for(auto e : data.outEdges[node]) {
				int to = e.other(node);
				int newDis = d;
				//if(s.covered[e.index])
					newDis += e.duration;
				auto it = minDis.find(to);
				if(it == minDis.end() || newDis < it->second) {
					q.push(make_pair(-newDis, to));
					minDis[to] = newDis;
					parent[to] = node;
					parentEdge[to] = e;
				}
			}
		}
		if(failed)
			break;
	}
	return s;
}


// Solution based on constructing Eulerian paths
Solution eulerianSolver(TestCase data) {
	State s;
	s.currentCar = 0;
	s.covered.resize(data.streets.size());
	s.solution.cars.resize(data.cars);
	/*map<int, int> numIn;
	map<int, int> numOut;
	map<int, int> numUndirected;
	for(auto e : data.streets) {
		if(e.directed){
			++numOut[e.from];
			++numIn[e.to];
		}
		else{
			++numUndirected[e.from];
			++numUndirected[e.to];
		}
	}
	for(int i = 0; i < (int)data.inEdges.size(); i++) {
		cerr << numIn[i] << " " << numOut[i] << " " << numUndirected[i] << endl;
	}*/
	for(int c = 0; c < data.cars; c++) {
		//cerr << "Car " << c << endl;
		if(false){
			s.currentCar = c;
			s.solution.cars[s.currentCar].junctions.push_back(data.startIndex);
			EulerGraph eulerGraph(data, s);
			eulerGraph.greedy();
			eulerGraph.optimize();
			vector<Street> sequence = eulerGraph.getCycle();
			int curNode = data.startIndex;
			double remainingTime = data.timeLimit;
			bool enoughTime = true;
			for(auto e : sequence) {
				if(e.duration > remainingTime) {
					enoughTime = false;
					break;
				}
				remainingTime -= e.duration;
				curNode = e.other(curNode);
				s.covered[e.index] = true;
				s.solution.cars[s.currentCar].junctions.push_back(curNode);
			}
			if(enoughTime) {
				cerr << "Warning! Car " << c << " finished its cycle" << endl;
			}
			continue;
		}
		s.currentCar = c;
		s.currentCarLocation = data.startIndex;
		s.solution.cars[s.currentCar].junctions.push_back(data.startIndex);
		s = extendSolution(data, s);
		//cerr << "Remaining " << remainingTime << endl;
		/*for(int i = 0; i < (int)data.junctions.size(); i++) {
			int totOut = additionalOutEdges.size();
			int totIn = additionalInEdges.size();
			for(auto e : data.outEdges[i]){
				if(covered[e.index])
					continue;
				++totOut;
			}
			for(auto e : data.inEdges[i]){
				if(covered[e.index])
					continue;
				++totIn;
			}
			relDegree[i] = totOut - totIn;
		}
		for(int i = 0; i < (int)data.junctions.size(); i++) {
			while(relDegree[i]) {
				priority_queue<pair<int, int> > q;
				map<int, int> minDis;
				q.push(make_pair(0, i));
				minDis[i] = 0;
				while(!q.empty()){
					auto cur = q.top();
					q.pop();
					int d = -cur.first;
					int node = cur.second;
					if(relDegree[node]){

						break;
					}
					if(minDis[node] < d)
						continue;
					for(auto e : data.outEdges[node]) {
						int to = e.other(node);
						int newDis = d + e.duration;
						auto it = minDis.find(to);
						if(it == minDis.end() || newDis < it->second) {
							q.push(make_pair(-newDis, to));
							M[to] = newDis;
						}
					}
				}
			}
		}*/
	}
	cerr << "Score: " << s.score(data) << endl;
	totalScore = s.score(data);
	return s.solution;
}

int checkSolution(const TestCase& data, const Solution& solution) {
	set<Street> doneStreets;
	int totDuration = 0;
	for(auto c : solution.cars) {
		int duration = 0;
		for(int i = 1; i < (int)c.junctions.size(); ++i) {
			auto e = data.getEdge(c.junctions[i-1], c.junctions[i]);
			duration += e.duration;
			doneStreets.insert(e);
		}
		if(duration > data.timeLimit)
			assert(0);
		totDuration += duration;
	}
	int ret = 0;
	for(auto e : doneStreets)
		ret += e.length;
	cerr << totDuration << " ";
	return ret;
}

Solution optimizeSolution(const TestCase& data, Solution solution) {
	while(true) {
		EulerGraph graph(data, State());
		map<Street, int> timesTraversed;
		for(auto c : solution.cars) {
			for(int i = 1; i < (int)c.junctions.size(); ++i) {
				timesTraversed[data.getEdge(c.junctions[i-1], c.junctions[i])]++;
			}
		}
		bool improved = false;
		for(int cind = 0; cind < (int)solution.cars.size(); ++cind) {
			auto& c = solution.cars[cind];
			for(int i = 0; i < (int)c.junctions.size(); ++i) {
				int duration = 0;
				set<Street> newEdges;
				for(int j = i; j < (int)c.junctions.size()-1; ++j) {
					auto e = data.getEdge(c.junctions[j], c.junctions[j+1]);
					if(newEdges.count(e))
						break;
					newEdges.insert(e);
					duration += e.duration;
					assert(timesTraversed[e]);
					if(timesTraversed[e] == 1)
						break;
					set<int> target;
					target.insert(c.junctions[j+1]);
					Path p = graph.computeDistances(c.junctions[i], target)[c.junctions[j+1]];
					if(p.duration() < duration) {
						/*cerr << checkSolution(data, solution) << endl;
						cerr << "Decreased time by " << duration - p.duration() << endl;*/
						improved = true;
						vector<int> newJunctions;
						for(int k = 0; k <= i; ++k) {
							newJunctions.push_back(c.junctions[k]);
						}
						for(int k = 0; k < (int)c.junctions.size()-1; ++k) {
							--timesTraversed[data.getEdge(c.junctions[k], c.junctions[k+1])];
						}
						int cur = newJunctions.back();
						for(int k = 0; k < (int)p.edges.size(); ++k) {
							cur = p.edges[k].other(cur);
							newJunctions.push_back(cur);
						}
						for(int k = j+2; k < (int)c.junctions.size(); ++k) {
							newJunctions.push_back(c.junctions[k]);
						}
						c.junctions = newJunctions;
						for(int k = 0; k < (int)c.junctions.size()-1; ++k) {
							++timesTraversed[data.getEdge(c.junctions[k], c.junctions[k+1])];
						}
						break;
					}
				}
			}
			cerr << checkSolution(data, solution) << endl;
			State s;
			s.solution = solution;
			s.currentCar = cind;
			s.currentCarLocation = c.junctions.back();
			s.covered.resize(data.streets.size());
			for(auto c2 : solution.cars) {
				for(int i = 1; i < (int)c2.junctions.size(); ++i) {
					s.covered[data.getEdge(c2.junctions[i-1], c2.junctions[i]).index] = true;
				}
			}
			s = extendSolution(data, s);
			solution = s.solution;
			cerr << checkSolution(data, solution) << endl;
		}
		if(!improved)
			break;
	}
	return solution;
}

void convertToPDDL(const TestCase& data) {
	ofstream domain;
	domain.open("domain.pddl");
	domain << "(define (domain metrictest)" << endl;
    domain << "\t(:requirements :strips :action-costs)" << endl;
    domain << "\t(:predicates (done ?x) (at ?y))" << endl;
  	domain << endl;
   	domain << "\t(:functions " << endl;
    domain << "\t\t((total-cost) - number)" << endl;
    domain << "\t)" << endl;

	int edgeCount = 0;
	int actionCount = 0;
	for(auto e : data.streets) {
		domain << "\t(:action do" << actionCount++ << endl;
		domain << "\t\t:precondition (at node" << e.from << ")" << endl;
		domain << "\t\t:effect (and (done edge" << edgeCount << ") (increase (total-cost) " << e.duration << ") (not (at node" << e.from << ")) (at node" << e.to << "))" << endl;
		domain << "\t)" << endl;
		if(!e.directed) {
			domain << "\t(:action do" << actionCount++ << endl;
			domain << "\t\t:precondition (at node" << e.to << ")" << endl;
			domain << "\t\t:effect (and (done edge" << edgeCount << ") (increase (total-cost) " << e.duration << ") (not (at node" << e.to << ")) (at node" << e.from << "))" << endl;
			domain << "\t)" << endl;
		}
		++edgeCount;
	}
	domain << ")" << endl << endl;
	domain.close();

	ofstream problem;
	problem.open("problem.pddl");

	problem << "(define (problem problem_1)" << endl;
    problem << "\t(:domain metrictest)" << endl;
    problem << "\t(:objects foo)" << endl;
    problem << "\t(:init (= (total-cost) 0) (at node0))" << endl;
    problem << "\t(:goal (and";
	for(int i = 0; i < edgeCount; ++i) {
		problem << " (done edge" << i << ")";
	}
	problem << "))" << endl;
    problem << "\t(:metric minimize (total-cost))" << endl;
    problem << ")" << endl;
	problem.close();

}

int main(){
	auto testCase = parseTestCase();
	convertToPDDL(testCase);
	auto solution = Solution();
	//solution = optimizeSolution(testCase, solution);
	auto bestSolution = solution;
	ll sumScores = 0;
	ll bestScore = 0;
	ll numScores = 0;
	for(int i = 0; i < 5; i++) {
		solution = eulerianSolver(testCase);
		if(totalScore > 1915000) {
			solution = optimizeSolution(testCase, solution);
		}
		totalScore = checkSolution(testCase, solution);
		sumScores += totalScore;
		if(totalScore > bestScore) {
			bestScore = totalScore;
			bestSolution = solution;
		}
		++numScores;
		cerr << "Best: " << bestScore << endl;
		cerr << "Average: " << (sumScores)/numScores << endl << endl;
	}
	/*auto solution = bruteforce(testCase);
	auto bestSolution = solution;*/
	bestSolution.print();
}
