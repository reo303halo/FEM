#pragma once
#define _USE_MATH_DEFINES

#include <cmath>
#include <vector>
#include <list>

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <ttl/halfedge/HeTriang.h>
#include <ttl/halfedge/HeDart.h>
#include <ttl/halfedge/HeTraits.h>


enum ProblemType {
    LAPLACE,
    POISSON,
    HELMHOLTZ,
    EIGENVALUE
};

class FEMobject
{
public:
    FEMobject();

    // -------------------------------
    // Problem setup
    // -------------------------------
    void setProblemType(ProblemType pt);
    void setBoundaryFunctions(double(*gD)(double,double), double(*gN)(double,double), double kappa);
    void setSourceFunction(double(*f)(double,double));
    void setExactSolution(double(*exact)(double,double));
    void setHelmholtzParameter(double lambda);
    
    // -------------------------------
    // Mesh generation
    // -------------------------------
    void circleMesh(int n, int m, double radius);
    void squareMesh(int n, double h, Eigen::Vector2d origin);

    // -------------------------------
    // PDE & FEM matrices
    // -------------------------------
    Eigen::SparseMatrix<double> massMat(std::list<hed::Edge*> trilist, int np);
    Eigen::SparseMatrix<double> stiffMat(std::list<hed::Edge*> trilist, int np);
    Eigen::VectorXd loadVect(std::list<hed::Edge*> trilist, int np);
    Eigen::SparseMatrix<double> RobinMat(std::list<hed::Dart> boundary, int np);
    Eigen::VectorXd RobinVect(std::list<hed::Dart> boundary, int np);

    // -------------------------------
    // Solve
    // -------------------------------
    void solve();
    void assignContiguousIDs();
    int getDoFs() const; // the number of degrees of freedom (#DoF) for the current mesh
    double getError(); // L2 error

    // -------------------------------
    // Functions for boundary/forcing
    // -------------------------------
    double kappa(double x, double y);
    double gN(double x, double y); // Neumann BC
    double gD(double x, double y); // Dirichlet BC
    double f(double x, double y);
    double uexact(double x, double y);

    // -------------------------------
    // Utilities
    // -------------------------------
    double triarea(hed::Node* N1, hed::Node* N2, hed::Node* N3);
    Eigen::Vector2<Eigen::Vector3d> gradients(Eigen::Vector3d x, Eigen::Vector3d y, double area);

    // -------------------------------
    // Visualization
    // -------------------------------
    void visualization(const std::string &filename);

public:
    hed::Triangulation triang;
    ProblemType Ptype;
    Eigen::VectorXd zeta;
    std::unordered_map<hed::Node*, int> ID;

    // Mesh data as simple arrays (useful for assembly and exporting)
    std::vector<Eigen::Vector2d> V;        // vertex coordinates (x,y)
    std::vector<Eigen::Vector3i> T;        // triangles (indices into V, 0-based)
    std::vector<int> boundaryFlag;         // 1 if boundary node, 0 if interior


private:
    double gD_value;
    double gN_value;
    double f_value;
    double robin_kappa;
    double lambda_param;
};
