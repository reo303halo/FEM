#include "FEMobject.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <algorithm>
#include <Eigen/Eigenvalues>


// -------------------------------
// Constructor
// -------------------------------
FEMobject::FEMobject()
    : gD_value(0.0), gN_value(0.0), f_value(0.0),
      robin_kappa(0.0), lambda_param(0.0), Ptype(LAPLACE)
{
    // Initialize boundaryFlag and other members if needed
    boundaryFlag.clear();
    zeta = Eigen::VectorXd();
}


// -------------------------------
// Problem setup
// -------------------------------
void FEMobject::setProblemType(ProblemType pt){
    Ptype = pt;
}

void FEMobject::setHelmholtzParameter(double lambda) {
    lambda_param = lambda;
}


// -------------------------------
// Mesh generation
// -------------------------------
void FEMobject::circleMesh(int n, int m, double r){
    std::vector<hed::Node*>* nodes = new std::vector<hed::Node*>;;
    hed::Node* p = new hed::Node(0,0);
    nodes->push_back(p);

    for (int j = 0; j < n; j++) {
        double dist = r*(j+1)/n;
        for (int i = 0; i < m*(j+1); i++) {

            double a = (i*2*M_PI)/(m*(j+1));
            double x = cos(a) * dist;
            double y = sin(a) * dist;
            nodes->push_back(new hed::Node(x, y));
        }
    }
    triang.createDelaunay(nodes->begin(), nodes->end());
}

void FEMobject::squareMesh(int n, double h, Eigen::Vector2d origin)
{
    double ox = origin(0);
    double oy = origin(1);
    std::vector<hed::Node*>* nodes = new std::vector<hed::Node*>;
    
    for (int j = 0; j <= n; j++) {
        for (int i = 0; i <= n; i++) {
            float x = ox + i * h / n;
            float y = oy + j * h / n;
            nodes->push_back(new hed::Node(x, y));
        }
    }
    triang.createDelaunay(nodes->begin(), nodes->end());
}


// -------------------------------
// PDE & FEM matrices
// -------------------------------
Eigen::SparseMatrix<double>
FEMobject::massMat(std::list<hed::Edge*> trilist, int np)
{
    Eigen::SparseMatrix<double> M(np, np);

    for (auto K = trilist.begin(); K != trilist.end(); ++K) {
        hed::Edge* edge = *K;

        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();

        Eigen::Vector3<Eigen::Index> loc2glb = {
            static_cast<Eigen::Index>(ID[N1]),
            static_cast<Eigen::Index>(ID[N2]),
            static_cast<Eigen::Index>(ID[N3])
        };

        Eigen::Vector3d x(N1->x(), N2->x(), N3->x());
        Eigen::Vector3d y(N1->y(), N2->y(), N3->y());

        double area = triarea(N1, N2, N3);

        Eigen::Matrix3d mat;
        mat << 2, 1, 1,
               1, 2, 1,
               1, 1, 2;
        Eigen::Matrix3d MK = 1./12*mat*area;

        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                M.coeffRef(loc2glb[i], loc2glb[j]) += MK(i, j);
            }
        }
    }
    return M;
}


Eigen::SparseMatrix<double>
FEMobject::stiffMat(std::list<hed::Edge*> trilist, int np) {
    Eigen::SparseMatrix<double> A(np, np); // initialization of the stiffness matrix A

    for (auto K = trilist.begin(); K != trilist.end(); K++) { // Loop over triangles
        hed::Edge* edge = *K;

        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();

        Eigen::Vector3<Eigen::Index> loc2glb = {
            static_cast<Eigen::Index>(ID[N1]),
            static_cast<Eigen::Index>(ID[N2]),
            static_cast<Eigen::Index>(ID[N3])
        };

        Eigen::Vector3d x(N1->x(), N2->x(), N3->x());
        Eigen::Vector3d y(N1->y(), N2->y(), N3->y());

        double area = triarea(N1, N2, N3);
        Eigen::Vector2<Eigen::Vector3d> gradphi = gradients(x, y, area);

        Eigen::Vector3d b = gradphi(0);
        Eigen::Vector3d c = gradphi(1);
        Eigen::Matrix3d AK = (b * b.transpose() + c * c.transpose()) * area;

        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                A.coeffRef(loc2glb[i], loc2glb[j]) += AK(i, j);
            }
        }
    }
    return A;
}

Eigen::VectorXd FEMobject::loadVect(std::list<hed::Edge*> trilist, int np) {
    Eigen::VectorXd b = Eigen::VectorXd::Zero(np);

    for (auto K = trilist.begin(); K != trilist.end(); ++K) {
        hed::Edge* edge = *K;

        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();
        
        Eigen::Vector3<Eigen::Index> loc2glb = {
            static_cast<Eigen::Index>(ID[N1]),
            static_cast<Eigen::Index>(ID[N2]),
            static_cast<Eigen::Index>(ID[N3])
        };

        Eigen::Vector3d x(N1->x(), N2->x(), N3->x());
        Eigen::Vector3d y(N1->y(), N2->y(), N3->y());

        double area = triarea(N1, N2, N3);
        double xc = (x(0) + x(1) + x(2)) / 3;
        double yc = (y(0) + y(1) + y(2)) / 3;
        double F = f(xc, yc);

        Eigen::Vector3d vec = {1, 1, 1};
        Eigen::Vector3d bK = F/3 * vec * area;
        
        for (int j = 0; j < 3; ++j) {
            b.coeffRef(loc2glb(j)) += bK(j);
        }
    }
    return b;
}


Eigen::SparseMatrix<double>
FEMobject::RobinMat(std::list<hed::Dart> boundary, int np) {
    Eigen::SparseMatrix<double> R(np, np);
    hed::Edge* edge = triang.getBoundaryEdge();
    hed::Dart b_dart(edge);
    ttl::getBoundary(b_dart, boundary);
    std::list<hed::Dart>::iterator E;
    
    for (E = boundary.begin(); E != boundary.end(); E++){
        hed::Node* N1 = E->getNode();
        hed::Node* N2 = E->getOppositeNode();
        
        Eigen::Vector2<Eigen::Index> loc2glb = {
            static_cast<Eigen::Index>(ID[N1]),
            static_cast<Eigen::Index>(ID[N2])
        };
        
        Eigen::Vector2d x = {N1->x(), N2->x()};
        Eigen::Vector2d y = {N1->y(), N2->y()};
        
        double len = sqrt((x(0) - x(1))*(x(0) - x(1)) + (y(0) - y(1))*(y(0) - y(1)));
        
        double xc = (x(0) + x(1))/2;
        double yc = (y(0) + y(1))/2;
        double k = FEMobject::kappa(xc, yc);
        
        Eigen::Matrix2d mat;
        mat << 2, 1,
               1, 2;
        Eigen::Matrix2d RE = k/6*mat*len;
        
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                R.coeffRef(loc2glb(i), loc2glb(j)) += RE(i, j);
            }
        }
    }
    return R;
}

Eigen::VectorXd FEMobject::RobinVect(std::list<hed::Dart> boundary, int np) {
    Eigen::VectorXd r = Eigen::VectorXd::Zero(np);
    hed::Edge* edge = triang.getBoundaryEdge();
    hed::Dart b_dart(edge);
    ttl::getBoundary(b_dart, boundary);
    std::list<hed::Dart>::iterator E;

    for (E = boundary.begin(); E != boundary.end(); E++){
        hed::Node* N1 = E->getNode();
        hed::Node* N2 = E->getOppositeNode();

        Eigen::Vector2<Eigen::Index> loc2glb = {
            static_cast<Eigen::Index>(ID[N1]),
            static_cast<Eigen::Index>(ID[N2])
        };

        Eigen::Vector2d x = {N1->x(), N2->x()};
        Eigen::Vector2d y = {N1->y(), N2->y()};

        double len = sqrt((x(0) - x(1))*(x(0) - x(1)) + (y(0) - y(1))*(y(0) - y(1)));

        double xc = (x(0) + x(1))/2;
        double yc = (y(0) + y(1))/2;
        double tmp = kappa(xc, yc)*gD(xc, yc)+gN(xc, yc);

        Eigen::Vector2d vec = {1, 1};
        Eigen::Vector2d rE = tmp*vec*len/2;

        for (int i = 0; i < 2; i++) {
            r.coeffRef(loc2glb(i)) += rE(i);
        }
    }
    return r;
}


// -------------------------------
// Solve
// -------------------------------
void FEMobject::solve() {
    // Get node list and define np (number of degrees of freedom)
    auto nodes_ptr = triang.getNodes();
    int np = static_cast<int>(nodes_ptr->size());

    // Build contiguous IDs and check they match np
    assignContiguousIDs();
    if ((int)ID.size() != np) {
        std::cerr << "Warning: ID.size() (" << ID.size()
                  << ") != np (" << np << ")\n";
    }

    // Initialize solution vector with correct size BEFORE any indexing
    zeta = Eigen::VectorXd::Zero(np);

    // -------------------------------
    // Get triangle and boundary lists
    // -------------------------------
    std::list<hed::Edge*> trilist = triang.getLeadingEdges();
    std::list<hed::Dart> boundary;
    hed::Edge* edge = triang.getBoundaryEdge();
    hed::Dart b_dart(edge);
    ttl::getBoundary(b_dart, boundary);

    // -------------------------------
    // Assemble FEM matrices (these functions must use ID[node] internally)
    // -------------------------------
    Eigen::SparseMatrix<double> A = stiffMat(trilist, np);  // stiffness matrix
    Eigen::SparseMatrix<double> M = massMat(trilist, np);   // mass matrix
    Eigen::SparseMatrix<double> R = RobinMat(boundary, np); // Robin boundary
    Eigen::VectorXd r = RobinVect(boundary, np);            // Robin vector
    Eigen::VectorXd b = loadVect(trilist, np);              // domain load vector

    // -------------------------------
    // Solve based on problem type
    // -------------------------------
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    switch (Ptype) {
        case LAPLACE: {
            solver.compute(A + R);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Factorization failed in LAPLACE\n";
                return;
            }
            zeta = solver.solve(r);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Solve failed in LAPLACE\n";
                return;
            }

            // assign solution back to nodes using ID map
            for (auto node : *nodes_ptr) {
                int idx = ID[node];
                if (idx < 0 || idx >= zeta.size()) {
                    std::cerr << "Index out-of-range when writing solution to node (LAPLACE): " << idx << "\n";
                    continue;
                }
                node->init(node->x(), node->y(), zeta(idx));
            }
            break;
        }

        case POISSON: {
            Eigen::SparseMatrix<double> S = A + R;
            solver.compute(S);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Factorization failed in POISSON solve\n";
                return;
            }
            zeta = solver.solve(b + r);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Linear solve failed in POISSON\n";
                return;
            }

            for (auto node : *nodes_ptr) {
                int idx = ID[node];
                if (idx < 0 || idx >= zeta.size()) {
                    std::cerr << "Index out-of-range when writing solution to node (POISSON): " << idx << "\n";
                    continue;
                }
                node->init(node->x(), node->y(), zeta(idx));
            }
            break;
        }

        case HELMHOLTZ: {
            Eigen::SparseMatrix<double> H = A + R - lambda_param * M;
            solver.compute(H);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Helmholtz: factorization failed\n";
                return;
            }
            zeta = solver.solve(b + r);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Helmholtz: solve failed\n";
                return;
            }
            for (auto node : *nodes_ptr) {
                int idx = ID[node];
                if (idx < 0 || idx >= zeta.size()) {
                    std::cerr << "Index out-of-range when writing solution to node (HELMHOLTZ): " << idx << "\n";
                    continue;
                }
                node->init(node->x(), node->y(), zeta(idx));
            }
            break;
        }

        case EIGENVALUE: {
            Eigen::MatrixXd A_dense = Eigen::MatrixXd(A + R);
            Eigen::MatrixXd M_dense = Eigen::MatrixXd(M);

            Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> ges;
            ges.compute(A_dense, M_dense);
            if (ges.info() != Eigen::Success) {
                std::cerr << "Eigenvalue solver failed\n";
                return;
            }

            Eigen::MatrixXd eigvecs = ges.eigenvectors();
            // pick principal eigenvector (first column)
            if (eigvecs.cols() > 0) {
                Eigen::VectorXd v0 = eigvecs.col(0);
                // ensure sizes match
                if (v0.size() == zeta.size()) {
                    zeta = v0;
                } else {
                    std::cerr << "Eigenvector size mismatch: " << v0.size() << " vs zeta " << zeta.size() << "\n";
                    zeta = Eigen::VectorXd::Zero(zeta.size());
                }
            }

            for (auto node : *nodes_ptr) {
                int idx = ID[node];
                if (idx < 0 || idx >= zeta.size()) {
                    std::cerr << "Index out-of-range when writing eigenvector to node: " << idx << "\n";
                    continue;
                }
                node->init(node->x(), node->y(), zeta(idx));
            }
            break;
        }

        default:
            std::cerr << "Unknown problem type in solve()\n";
            break;
    }
}


void FEMobject::assignContiguousIDs() {
    ID.clear();
    auto nodes = triang.getNodes();
    int k = 0;
    for (auto* node : *nodes) {
        ID[node] = k++;
    }
}

int FEMobject::getDoFs() const {
    return static_cast<int>(triang.getNodes()->size());
}

double FEMobject::getError() {
    double error = 0.0;
    std::list<hed::Edge*> trilist = triang.getLeadingEdges();

    // quick safety: make sure zeta has a reasonable size
    if (zeta.size() == 0) {
        std::cerr << "getError(): zeta is empty — run solve() first.\n";
        return std::numeric_limits<double>::infinity();
    }

    for (hed::Edge* edg : trilist) {
        hed::Node* N1 = edg->getSourceNode();
        hed::Node* N2 = edg->getTargetNode();
        hed::Node* N3 = edg->getNextEdgeInFace()->getTargetNode();

        double x1 = N1->x(); double x2 = N2->x(); double x3 = N3->x();
        double y1 = N1->y(); double y2 = N2->y(); double y3 = N3->y();

        double x_c = (x1 + x2 + x3) / 3.0;
        double y_c = (y1 + y2 + y3) / 3.0;

        // Convert nodes to our contiguous indices
        auto it1 = ID.find(N1);
        auto it2 = ID.find(N2);
        auto it3 = ID.find(N3);

        if (it1 == ID.end() || it2 == ID.end() || it3 == ID.end()) {
            std::cerr << "getError(): node not in ID map (skipping triangle)\n";
            continue;
        }

        int id1 = it1->second;
        int id2 = it2->second;
        int id3 = it3->second;

        // Debug print if any index out of bounds
        if (id1 < 0 || id1 >= zeta.size() || id2 < 0 || id2 >= zeta.size() || id3 < 0 || id3 >= zeta.size()) {
            std::cerr << "getError(): index out of range: id1=" << id1
                      << " id2=" << id2 << " id3=" << id3
                      << " zeta.size()=" << zeta.size() << "\n";
            // skip this triangle instead of crashing (temporary — you can abort after debugging)
            continue;
        }

        double u_h = (zeta[id1] + zeta[id2] + zeta[id3]) / 3.0;
        double u_exact = uexact(x_c, y_c);

        double area = triarea(N1, N2, N3);
        double diff = u_exact - u_h;
        error += (diff * diff) * area;
    }
    return std::sqrt(error);
}



// -------------------------------
// Functions for boundary/forcing
// -------------------------------
double FEMobject::kappa(double x, double y) {
    switch (Ptype) {
        case LAPLACE:
            return std::pow(10.0, 6.0);
            break;
            
        case POISSON:
            return std::pow(10.0, 6.0);
            break;
            
        case HELMHOLTZ:
            if (x > 0.0) {
                return 0.0;
            } else { return std::pow(10.0, 6.0);}
            break;
            
        case EIGENVALUE:
            if (x > 0.0) {
                return 0.0;
            } else { return std::pow(10.0, 6.0);}
            break;
            
        default:
            break;
    }
    return 0;
}

double FEMobject::gN(double x, double y) {
    switch (Ptype) {
        case LAPLACE:
            return 0.0;
            break;
            
        case POISSON:
            return 0.0;
            break;
            
        case HELMHOLTZ:
            return 0.0;
            break;
            
        case EIGENVALUE:
            return 0.0;
            break;
            
        default:
            break;
    }
    return 0;
}

double FEMobject::gD(double x, double y) {
    switch (Ptype) {
        case LAPLACE:
            return cos(4*atan2(y,x));
            break;
            
        case POISSON:
            return std::pow(y, 2.0)/2;
            break;
            
        case HELMHOLTZ:
            return 0.25;
            break;
            
        case EIGENVALUE:
            return 0.0;
            break;
            
        default:
            break;
    }
    return 0;
}

double FEMobject::f(double x, double y) {
    switch (Ptype) {
        case LAPLACE:
            return 0.0;
            break;
            
        case POISSON:
            return 1.0;
            break;
            
        case HELMHOLTZ:
            return 0.0;
            break;
            
        case EIGENVALUE:
            return 0.0;
            break;
            
        default:
            break;
    }
    return 0;
}

double FEMobject::uexact(double x, double y) {
    switch (Ptype) {
        case LAPLACE:
            return std::cos(4.0 * std::atan2(y, x));

        case POISSON:
            return 0.5 * (1.0 - x*x);

        case HELMHOLTZ: {
            double k = std::sqrt(lambda_param);
            return 0.25 * (std::cos(k * x) + std::tan(k) * std::sin(k * x));
        }

        case EIGENVALUE:
            return 0.0;

        default:
            return 0.0;
    }
}


// -------------------------------
// Utilities
// -------------------------------
double FEMobject::triarea(hed::Node* N1, hed::Node* N2, hed::Node* N3) {
    Eigen::Vector3d a(N2->x()-N1->x(), N2->y()-N1->y(), 0);
    Eigen::Vector3d b(N3->x()-N1->x(), N3->y()-N1->y(), 0);
    double area = 0.5 * (b.cross(a)).norm();
    return area;
}

Eigen::Vector2<Eigen::Vector3d>
FEMobject::gradients(Eigen::Vector3d x, Eigen::Vector3d y, double area) {
    Eigen::Vector3d b((y(1)-y(2))/(2*area),(y(2)-y(0))/(2*area),(y(0)-y(1))/(2*area));
    Eigen::Vector3d c((x(2)-x(1))/(2*area),(x(0)-x(2))/(2*area),(x(1)-x(0))/(2*area));
    Eigen::Vector2<Eigen::Vector3d> gradphi = {b, c};
    return gradphi;
}


// -------------------------------
// Visualization
// -------------------------------
void FEMobject::visualization(const std::string &filename) {
    std::ofstream objfile;
    objfile.open(filename);
    std::map<int, int > indexmap;
    int current_ind = 1;
    
    std::unordered_map<hed::Node*, int> nodeIndex;
    std::unordered_map<hed::Edge*, int> edgeNormalIndex;
    
    std::list<hed::Node*>* nodelist = triang.getNodes();
    std::list<hed::Node*>::iterator L;
    int idx = 1;
    for (L = nodelist->begin(); L != nodelist->end(); L++) {
        hed::Node* node = *L;
        objfile << "v " << node->x() << " " << node->z() << " " << node->y() << "\n";
        nodeIndex[node] = idx;
        idx++;
    }
    
    std::list<hed::Edge*> trilist = triang.getLeadingEdges();
    std::list<hed::Edge*>::iterator K;
    
    int normind = 1;
    for (K = trilist.begin(); K != trilist.end(); K++) {
        hed::Edge* edge = *K;
        
        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();
        
        Eigen::Vector3d a(N2->x()-N1->x(), N2->y() - N1->y(), N2->z() - N1->z());
        Eigen::Vector3d b(N3->x()-N1->x(), N3->y() - N1->y(), N3->z() - N1->z());
        Eigen::Vector3d normal = a.cross(b);
        objfile << "vn " << normal.x() << " " << normal.z() << " " << normal.y() << "\n";
        edgeNormalIndex[edge] = normind;
        normind++;
    }
    
    for (K = trilist.begin(); K != trilist.end(); K++) {
        hed::Edge* edge = *K;
        
        hed::Node* N1 = edge->getSourceNode();
        hed::Node* N2 = edge->getTargetNode();
        hed::Node* N3 = edge->getNextEdgeInFace()->getTargetNode();
        int ni = edgeNormalIndex[edge];
        
        objfile << "f " <<
        nodeIndex[N1] << "//" << ni << " " <<
        nodeIndex[N2] << "//" << ni << " " <<
        nodeIndex[N3] << "//" << ni << "\n";
    }
    objfile.close();
}
