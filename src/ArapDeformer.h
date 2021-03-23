#pragma once

#include <iostream>
#include <fstream>
#include <Eigen/Sparse>

#include "Eigen.h"
#include "SimpleMesh.h"
// #include <omp.h>
using namespace std;

#define THRESHOLD 1.0e-3f // Threshold for ending optimization

//The fixed points and handles are represented as constraints
struct Constraint {
	int vertexID;
	Vector3d position;
};


class ArapDeformer {

public:
    ~ArapDeformer();
    ArapDeformer();
	ArapDeformer(SimpleMesh *mesh, int weight_type, int estimation_type);
    // void initDeformation(vector<int> fixed_points);
    void setHandleConstraint(int handleID, Vector3d newHandlePosition);
    void applyDeformation(vector<int> fixed_points, int handleID, Vector3d handleNewPosition, int iterations);
    SimpleMesh m_mesh;

private:
    void estimateRotation();
    void updateB();
    void estimateVertices();
    double calculateEnergy();
    void buildWeightMatrix(); 
    void calculateSystemMatrix();
    Vector3d getConstraintI(int id);
    bool isInConstraints(int i);
    void updateSystemMatrix();
    void setWeightType(int weight_type);
    void setDecompositionType(int estimation_type);


    vector<MatrixXd> m_cell_rotations;
    MatrixXd m_system_matrix;
    MatrixXd m_system_matrix_original;
    MatrixXd m_weight_matrix;
    SparseMatrix<double> m_system_matrix_sparse;
    int m_num_v;
    int m_num_p;
    MatrixXd m_b;
    int m_handle_id;
    Vector3d m_new_handle_position;
    vector<Constraint> m_constraints;
    bool use_uniform_weights;
    bool use_constant_weights;
    bool use_cotangent_weights;
    bool use_sparse_qr;
    bool use_sparse_lu;
    bool use_sparse_matrices;

};

