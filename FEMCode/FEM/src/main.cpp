#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include "FEMobject.h"

int main() {
    int solution = 1;   // 1 = Laplace, 2 = Poisson, 3 = Helmholtz, 4 = Eigenvalue
    FEMobject model;
    
    const int refinements = 5;
    int n = 8;           // initial mesh parameter
    int m = 4;           // circleMesh argument
    double r = 1.0;

    std::vector<int> dofs_list;
    std::vector<double> error_list;

    switch(solution)
    {
        case 1: {
            // -------------------------------
            // Solution 1: Circle mesh, Laplace
            for (int k = 0; k < refinements; ++k) {
                model.circleMesh(n, m, r);
                model.setProblemType(ProblemType::LAPLACE);
                
                model.solve();
                
                //std::string filename = "Laplace_" + std::to_string(n_square) + ".obj";
                //model.visualization(filename);
                
                int dofs = model.getDoFs();
                double err = model.getError();
                
                dofs_list.push_back(dofs);
                error_list.push_back(err);
                
                std::cout << "Refinement " << k
                << ": DoFs = " << dofs
                << ", L2 error = " << err << std::endl;
                
                n *= 2;   // refine by doubling resolution
            }
            break;
        }
        case 2: {
            // -------------------------------
            // Solution 2: Circle mesh, Poisson
            for (int k = 0; k < refinements; ++k) {
                model.circleMesh(n, m, r);
                model.setProblemType(ProblemType::POISSON);
                
                model.solve();
                
                //std::string filename = "Poisson_" + std::to_string(n_square) + ".obj";
                //model.visualization(filename);
                
                int dofs = model.getDoFs();
                double err = model.getError();
                
                dofs_list.push_back(dofs);
                error_list.push_back(err);
                
                std::cout << "Refinement " << k
                << ": DoFs = " << dofs
                << ", L2 error = " << err << std::endl;
                
                n *= 2;   // refine by doubling resolution
            }
            break;
        }
        case 3: {
            // -------------------------------
            // Solution 3: Square mesh, Helmholtz
            for (int k = 0; k < refinements; ++k) {
                int n_square = 8 * (1 << k);
                model.squareMesh(n_square, 1.0, {0.0, -0.5});
                model.setHelmholtzParameter(81);
                model.setProblemType(ProblemType::HELMHOLTZ);
                
                model.solve();
                
                //std::string filename = "Helmholtz_" + std::to_string(n_square) + ".obj";
                //model.visualization(filename);
                
                int dofs = model.getDoFs();
                double err = model.getError();
                
                dofs_list.push_back(dofs);
                error_list.push_back(err);
                
                std::cout << "Refinement " << k
                << ": DoFs = " << dofs
                << ", L2 error = " << err << std::endl;
            }
            break;
        }
        case 4: {
            // -------------------------------
            // Solution 4: Square mesh, Eigenvalue
            /*
            model.squareMesh(6, 1.0, {0.0, -0.5});
            model.setProblemType(ProblemType::EIGENVALUE);
            model.solve();
            model.visualization("Eigenvalue_6.obj");
            
            model.squareMesh(13, 1.0, {0.0, -0.5});
            model.setProblemType(ProblemType::EIGENVALUE);
            model.solve();
            model.visualization("Eigenvalue_13.obj");
             */

            model.squareMesh(16, 1.0, {0.0, -0.5});
            model.setProblemType(ProblemType::EIGENVALUE);
            model.solve();
            model.visualization("Eigenvalue_16.obj");
           
            break;
        }
        default:
            std::cerr << "Invalid mode selected!\n";
            break;
    }
    
    // Print to console and writes.
    if (solution != 4) {
        std::cout << "\nData for plotting (DoFs, Error):\n";
        
        std::string fileName = "convergence_solution_" + std::to_string(solution) + ".dat";
        std::ofstream fout(fileName);
        if (!fout.is_open()) {
            std::cerr << "Could not open convergence.dat for writing!\n";
        }
        
        for (int i = 0; i < dofs_list.size(); ++i) {
            
            std::cout << dofs_list[i] << " " << error_list[i] << "\n";
            
            if (fout.is_open()) {
                fout << dofs_list[i] << " " << error_list[i] << "\n";
            }
        }
        if (fout.is_open()) fout.close();
    }

    return 0;
}
