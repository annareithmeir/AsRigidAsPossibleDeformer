#define NOMINMAX
#include <Windows.h>
#include <iostream>
#include <fstream>
#include <chrono>

#include "Eigen.h"
#include "SimpleMesh.h"
#include <Eigen/Sparse>
#include "GUI.h"
#include "ArapDeformer.h"


using namespace std;

vector<int> pointsT300() {
	vector<int> fixedPoints = { 58,31,32,104,106,167,177,181,226,258,261,262,271,283 };
	return fixedPoints;
}

vector<int> pointsT750() {
	vector<int> fixedPoints = { 111,46,199,216,240,261,386,397,414,481,492,518,579,589,595,635,646,682,695,738,740 };
	return fixedPoints;
}

vector<int> pointsT1000() {
	vector<int> fixedPoints = { 68,2,65,67,136,181,216,254,280,296,336,364,392,414,456,457,515,570,738,773,886,934,944,991 };
	return fixedPoints;
}

vector<int> pointsT1500() {
	vector<int> fixedPoints = { 146,1,96,102,108,159,192,270,275,398,417,435,495,499,503,537,605,673,754,832,1029,1099,1113,1224,1233,1236,1237,1331,1333,1406,1407,1488 };
	return fixedPoints;
}

vector<int> pointsT3000() {
	vector<int> fixedPoints = { 178,4,13,205,215,355,430,647,679,757,779,826,830,832,841,908,989,1045,1048,1202,1290,1303,1426,1433,1724,1788,1822,2026,2037,2081,2082,2246,2264,2491,2635,2684,2818,2841,2972 };
	return fixedPoints;
}

vector<int> squarePoints() {
	vector<int> fixedPoints;
	fixedPoints.push_back(220);
	for (int i = 0; i <= 20; i++) {
		fixedPoints.push_back(i);
		fixedPoints.push_back(420+i);
	}
	return fixedPoints;
}

vector<int> testSmall() {
	vector<int> fixedPoints;
	fixedPoints.push_back(9);
	fixedPoints.push_back(0);
	fixedPoints.push_back(5);
	return fixedPoints;
}

void runTest(vector<int> fixedPoints, int handleID, Vector3d handleMoved, vector<Vector3d> target_positions, int numIterations, int weight_type, int estimation_type, string filenameMesh, string filename) {
	SimpleMesh source_mesh;
	if (!source_mesh.loadMesh(filenameMesh)) {
		std::cout << "Mesh file wasn't read successfully at location: " << filenameMesh << std::endl;
		return;
	}

	ArapDeformer deformer(&source_mesh, weight_type, estimation_type);

	auto t1 = std::chrono::high_resolution_clock::now();
	deformer.applyDeformation(fixedPoints, handleID, handleMoved, numIterations);
	auto t2 = std::chrono::high_resolution_clock::now();

	std::chrono::duration<double> eps = t2 - t1;
	std::cout << "Deformation completed in " << eps.count() << " seconds." << std::endl;

	string weight = "";
	switch (weight_type)
	{
	case 0: weight = "uniform"; break;
	case 1: weight = "constant"; break;
	case 2: weight = "cotangent"; break;
	default:
		weight = "cotangent"; break;
	}

	string estimation = "";
	switch (estimation_type)
	{
	case 0: estimation = "QR"; break;
	case 1: estimation = "LU"; break;
	case 2: estimation = "DENSE"; break;
	default:
		estimation = "LU"; break;
	}

	string jump_type = "none";
	if (handleMoved == target_positions[0]) jump_type = "small";
	else if (handleMoved == target_positions[1]) jump_type = "large";


	string testname = filename + "_" + jump_type + "_" + weight + "_" + estimation + "_" + to_string(numIterations) + "_iterations";
	deformer.m_mesh.writeMesh("../results/" + testname + ".off");
	std::cout << "Deformation for " << testname << " completed in " << eps.count() << " seconds." << std::endl;
	std::cout << "##############################################" << std::endl;
}

int main() {
	/***
	 *  SimpleMesh.m_vertices -> original positions
	 *  SimpleMesh.m_verticesPrime -> computed new position in deformed mesh
	 * 	ARAP.deform -> estimateVertices() löst das LGS Lx=b, wobei L=SimpleMesh.m_systemMatrix, x=pprime, b=b im code
	 *
	 * 	test.off
	 * 	0--1--2--3--4  --> y
	 *  | /| /| /| /|
	 * 	5--6--7--8--9
	 *
	 * 	|
	 *  x
	 *
	 ***/

	int numIterations = 3;
	int testcase = 0;
	bool testing = false;
	bool weight_test = false;

	// Load the source and target mesh.
	bool isInvalid = true;
	std::string filenameMesh, filename;
	do {
		isInvalid = true;
		filenameMesh = std::string("../data/");
		std::cout << "Insert filename: " << std::endl;
		std::cin >> filename;
		filenameMesh.append(filename);
		filenameMesh.append(".off");

		std::ifstream file(filenameMesh);
		if (file.is_open()) {
			isInvalid = false;
		}
		if (isInvalid) {
			std::cout << "Mesh file does not exist." << std::endl;
			std::cin.clear();
			std::cin.ignore();
		}
	} while (isInvalid);

	int weight_type, estimation_type;

	do {
		isInvalid = true;
		std::cout << "Choose type of weighting function (0 = uniform, 1 = constant, 2 = cotangent)" << std::endl;
		std::cin >> weight_type;

		if (!std::cin.fail()) {
			if (weight_type >= 0 && weight_type <= 2) {
				isInvalid = false;
			}
		}
		if (isInvalid) {
			std::cout << "Invalid input!" << std::endl;
			std::cin.clear();
			std::cin.ignore();
		}
	} while (isInvalid);

	do {
		isInvalid = true;
		std::cout << "Choose type of matrix decomposition (0 = sparse QR, 1 = sparse LU, 2 = non sparse matrices)" << std::endl;
		std::cin >> estimation_type;

		if (!std::cin.fail()) {
			if (estimation_type >= 0 && estimation_type <= 2) {
				isInvalid = false;
			}
		}
		if (isInvalid) {
			std::cout << "Invalid input!" << std::endl;
			std::cin.clear();
			std::cin.ignore();
		}
	} while (isInvalid);

	if (!testing) {
		GUI* gui = new GUI(filenameMesh, numIterations, weight_type, estimation_type);
	}
	else {

		vector<int> fixedPoints;

		Vector4i iter(3, 10, 25, 50);

		if (filename.compare("armadillo_300") == 0) testcase = 0;
		else if (filename.compare("armadillo_750") == 0) testcase = 1;
		else if (filename.compare("armadillo_1k") == 0) testcase = 2;
		else if (filename.compare("armadillo_1500") == 0) testcase = 3;
		else if (filename.compare("armadillo_6k") == 0) testcase = 4;
		else if (filename.compare("square_21") == 0) testcase = 5;
		else if (filename.compare("test") == 0) testcase = 6;


		//´set fixed points
		switch (testcase)
		{
		case 0: fixedPoints = pointsT300(); break;
		case 1: fixedPoints = pointsT750(); break;
		case 2: fixedPoints = pointsT1000(); break;
		case 3: fixedPoints = pointsT1500(); break;
		case 4: fixedPoints = pointsT3000(); break;
		case 5: fixedPoints = squarePoints(); break;
		case 6: fixedPoints = testSmall(); break;
		default:
			fixedPoints = pointsT300(); break;
		}

		//we always wrote the handle first into the vector
		int handleID = fixedPoints[0];

		vector<Vector3d> targetPositions = { Vector3d(-35.093834, 91.842247,-34.709782), Vector3d(93.464325, 124.944626, -40.314186) };
		if (!weight_test) {
			for (Vector3d handleMoved : targetPositions) {
				for (int i = 0; i < iter.size(); i++) {
					//for (int w = 0; w < 3; w++) {
					for (int e = 0; e < 3; e++) {
						runTest(fixedPoints, handleID, handleMoved, targetPositions, iter[i], weight_type, e, filenameMesh, filename);
					}
					//}
				}
			}
		}
		else {
			targetPositions = {Vector3d(0.5, 0, 0.1), Vector3d(0.5, 0, 1)};//{ Vector3d(2, 5, 0), Vector3d(2, 3, 1) }; //
			for (Vector3d handleMoved : targetPositions) {
				for (int i = 0; i < iter.size(); i++) {
					for (int w = 0; w < 3; w++) {
						runTest(fixedPoints, handleID, handleMoved, targetPositions, iter[i], w, estimation_type, filenameMesh, filename);
					}
				}
			}
		}
	}
	return 0;
}