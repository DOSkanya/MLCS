#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>
#include <map>

#define DEBUG
//#undef DEBUG

#define SEQUENCE_LENGTH 120
#define BASES_COUNT 4

enum class NitrogenousBases
{
	A, G, C, T
};

class MatchPoint {
public:
	MatchPoint(int matchPointLength, int initValue) {
		positions.resize(matchPointLength, initValue);
	}

	MatchPoint(MatchPoint& mp, int index, int*** successorTable) {
		positions.resize(mp.positions.size(), 0);
		for (int i = 0; i < positions.size(); i++) {
			positions[i] = successorTable[i][index][mp.positions[i]];
		}
	}

	bool isValid() {
		for (int i = 0; i < positions.size(); i++) {
			if (positions[i] == 0) return false;
		}
		return true;
	}

	bool operator < (const MatchPoint& mp) const {
		if (this->positions.size() == mp.positions.size()) {
			for (int i = 0; i < this->positions.size(); i++) {
				if (mp.positions[i] > this->positions[i]) return true;
				else if (mp.positions[i] < this->positions[i]) return false;
				else if (i == this->positions.size() - 1 && mp.positions[i] == this->positions[i]) return false;
			}
		}
		else if (this->positions.size() < mp.positions.size()) return true;
		else return false;
	}

	std::vector<int> positions;
};

class LeveledDAGNode {
public:
	LeveledDAGNode(int matchPointLength, int initValue) : matchPoint(matchPointLength, initValue), incidentEdgeCount(0) {}
	void GenerateSuccessor(int*** successorTable, std::map<MatchPoint, LeveledDAGNode*>& DAG, 
							std::vector<LeveledDAGNode*>& nextLevel, LeveledDAGNode* endNode) {
		for (int i = 0; i < BASES_COUNT; i++) {
			MatchPoint newMatchPoint(matchPoint, i, successorTable);
			if (newMatchPoint.isValid()) {
				auto it = DAG.find(newMatchPoint);
				if (it == DAG.end()) {
					LeveledDAGNode* newDAGNode = new LeveledDAGNode(matchPoint.positions.size(), 0);
					newDAGNode->matchPoint = newMatchPoint;
					std::string LCS;
					switch (i) {
					case 0: LCS = "a"; newDAGNode->symbol = 'a'; break;
					case 1: LCS = "g"; newDAGNode->symbol = 'g'; break;
					case 2: LCS = "c"; newDAGNode->symbol = 'c'; break;
					case 3: LCS = "t"; newDAGNode->symbol = 't'; break;
					}
					newDAGNode->partialLCS.push_back(LCS);
					newDAGNode->incidentEdgeCount++;
					DAG.insert(std::make_pair(newDAGNode->matchPoint, newDAGNode));
					successors.push_back(newDAGNode);
					nextLevel.push_back(newDAGNode);
				}
				else {
					successors.push_back(it->second);
					it->second->incidentEdgeCount++;
				}
			}
		}

		if (successors.size() == 0) {
			successors.push_back(endNode);
			endNode->incidentEdgeCount++;
		}
	}

	void RemovedFromDAG(std::map<MatchPoint, LeveledDAGNode*>& DAG) {
		for (int i = 0; i < successors.size(); i++) {
			int s_length = successors[i]->partialLCS.size() == 0 ? 0 : successors[i]->partialLCS[0].length();
			int p_length = partialLCS.size() == 0 ? 0 : partialLCS[0].length();
			
			if (p_length >= s_length) {
				successors[i]->partialLCS.clear();
				successors[i]->partialLCS.shrink_to_fit();
				for (int j = 0; j < partialLCS.size(); j++) {
					std::string str = partialLCS[j] + successors[i]->symbol;
					successors[i]->partialLCS.push_back(str);
				}
			}
			else if (p_length + 1 == s_length && p_length != 0) {
				for (int j = 0; j < partialLCS.size(); j++) {
					std::string str = partialLCS[j] + successors[i]->symbol;
					successors[i]->partialLCS.push_back(str);
				}
			}

			successors[i]->incidentEdgeCount--;
		}
	}

	void PrintMLCS() {
		for (int i = 0; i < partialLCS.size(); i++) {
			std::cout << partialLCS[i] << std::endl;
		}
	}

	int MLCSLength() {
		return partialLCS[0].length();
	}

	MatchPoint matchPoint;
	int incidentEdgeCount;
	char symbol;
private:
	std::vector<LeveledDAGNode*> successors;
	std::vector<std::string> partialLCS;
};

char** ReadSequenceFromFile(int number) {
	std::ifstream file = std::ifstream("DNA Sequence.txt", std::ios::in);
	if (!file.is_open()) {
		std::cout << "ERROR::Fail to load Sequence." << std::endl;
	}

	char** sequence = new char* [number];
	for (int i = 0; i < number; i++) {
		sequence[i] = new char[SEQUENCE_LENGTH];
	}

	for (int i = 0; i < number; i++) {
		char* temp = new char[SEQUENCE_LENGTH + 50];
		file.getline(temp, SEQUENCE_LENGTH + 50);
		memcpy(sequence[i], temp, SEQUENCE_LENGTH);
	}
	
	file.close();
	return sequence;
}

int*** ConstructSuccessorTable(int sequenceCount, char** sequence) {
	// Allocate memory required for successor table
	int*** successorTable = new int** [sequenceCount];
	for (int i = 0; i < sequenceCount; i++) {
		successorTable[i] = new int* [BASES_COUNT];
		for (int j = 0; j < BASES_COUNT; j++) {
			successorTable[i][j] = new int[SEQUENCE_LENGTH + 1];
			// Initialize table with all zero
			memset(successorTable[i][j], 0, sizeof(int) * (SEQUENCE_LENGTH + 1));
		}
	}

	// Fill the successor table
	for (int i = 0; i < sequenceCount; i++) {
		// Only scan the sequence once
		int pointer[BASES_COUNT] = { 0 };
		for (int j = 0; j < SEQUENCE_LENGTH; j++) {
			NitrogenousBases base = NitrogenousBases::A;
			switch (sequence[i][j]) {
			case 'a': base = NitrogenousBases::A; break;
			case 'g': base = NitrogenousBases::G; break;
			case 'c': base = NitrogenousBases::C; break;
			case 't': base = NitrogenousBases::T; break;
			}
			
			while (pointer[(int)base] <= j) {
				successorTable[i][(int)base][pointer[(int)base]] = j + 1;
				pointer[(int)base]++;
			}
		}
	}

	return successorTable;
}

void RemoveOutdatedNode(std::map<MatchPoint, LeveledDAGNode*>& DAG, LeveledDAGNode* endNode) {
	std::vector<LeveledDAGNode*> noIncidentEdgeNode;
	for (auto it = DAG.begin(); it != DAG.end(); it++) {
		if (it->second->incidentEdgeCount == 0 && it->second != endNode) noIncidentEdgeNode.push_back(it->second);
	}

	for (int i = 0; i < noIncidentEdgeNode.size(); i++) {
		noIncidentEdgeNode[i]->RemovedFromDAG(DAG);
		DAG.erase(noIncidentEdgeNode[i]->matchPoint);
		delete noIncidentEdgeNode[i];
	}
	
}

void SwapBetweenLevels(std::vector<LeveledDAGNode*>& currentLevel, std::vector<LeveledDAGNode*>& nextLevel) {
	currentLevel.swap(nextLevel);
	nextLevel.clear();
}

int main(int argc, char** argv) {
#ifdef DEBUG
	int sequenceCount = 4;
#else
	int sequenceCount = std::atoi(argv[1]);
#endif
	
	if (sequenceCount < 3 || sequenceCount > 10) {
		std::cout << "ERROR::Should be bigger than 3 and less or equal to 10." << std::endl;
	}

	char** sequence = ReadSequenceFromFile(sequenceCount);

	int*** successorTable = ConstructSuccessorTable(sequenceCount, sequence);

	std::map<MatchPoint, LeveledDAGNode*> leveledDAG;
	std::vector<LeveledDAGNode*> currentLevel;
	std::vector<LeveledDAGNode*> nextLevel;

	// Initialized before constructing leveled DAG
	LeveledDAGNode* sourceNode = new LeveledDAGNode(sequenceCount, 0);
	LeveledDAGNode* endNode = new LeveledDAGNode(sequenceCount, std::numeric_limits<int>::max());
	leveledDAG.insert(std::make_pair(sourceNode->matchPoint, sourceNode));
	leveledDAG.insert(std::make_pair(endNode->matchPoint, endNode));
	currentLevel.push_back(sourceNode);

	while (currentLevel.size() != 0) {
		for (int i = 0; i < currentLevel.size(); i++) {
			currentLevel[i]->GenerateSuccessor(successorTable, leveledDAG, nextLevel, endNode);
		}
		RemoveOutdatedNode(leveledDAG, endNode);
		SwapBetweenLevels(currentLevel, nextLevel);
	}

	while (leveledDAG.size() > 1) {
		RemoveOutdatedNode(leveledDAG, endNode);
	}

	endNode->PrintMLCS();
	std::cout << "Length: " << endNode->MLCSLength() << std::endl;

	return 0;
}