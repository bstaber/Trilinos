/*
Brian Staber (brian.staber@gmail.com)
*/

#include "plasticitySmallStrains.hpp"

/*plasticitySmallStrains::plasticitySmallStrains(){

}*/

plasticitySmallStrains::~plasticitySmallStrains(){

}

plasticitySmallStrains::plasticitySmallStrains(Epetra_Comm & comm, Teuchos::ParameterList & parameterlist){

  Delta         = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "delta");
  iter_min      = Teuchos::getParameter<int>(parameterlist.sublist("Newton"),    "iterMin");
  iter_max      = Teuchos::getParameter<int>(parameterlist.sublist("Newton"),    "iterMax");
  nb_bis_max    = Teuchos::getParameter<int>(parameterlist.sublist("Newton"),    "nbBisMax");
  norm_inf_tol  = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "NormFTol");
  norm_inf_max  = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "NormFMax");
  eps           = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "eps");
  success       = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "success_parameter");
  failure       = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "failure_parameter");
  bc_disp       = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "bc_disp");
  pressure_load = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "pressure_load");
  tol           = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "tol");

  Krylov = &parameterlist.sublist("Krylov");

  std::string mesh_file = Teuchos::getParameter<std::string>(parameterlist.sublist("Mesh"), "mesh_file");
  Mesh                  = new mesh(comm, mesh_file, 1.0);
  Comm                  = Mesh->Comm;

  StandardMap           = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
  OverlapMap            = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
  ImportToOverlapMap    = new Epetra_Import(*OverlapMap,*StandardMap);

  constructGaussMap(*GaussMap);

}

int plasticitySmallStrains::incremental_bvp(bool print){

    Epetra_Time Time(*Comm);

    Epetra_LinearProblem problem;
    AztecOO solver;

    Epetra_FECrsMatrix stiffness(Copy,*FEGraph);
    Epetra_FEVector    rhs(*StandardMap);
    Epetra_Vector      lhs(*StandardMap);
    Epetra_Vector      x(*StandardMap);
    Epetra_Vector      y(*StandardMap);
    Epetra_Vector      eto(*GaussMap);
    Epetra_Vector      ep_old(*GaussMap);

    double delta = Delta;
    double Assemble_time;
    double Aztec_time;
    double displacement;
    double norm_inf_rhs;
    double time_max = 1.0;
    double Krylov_res = 0;
    int FLAG1, FLAG2, FLAG3, nb_bis, iter;
    int Krylov_its = 0;

    time = 0.0;

    std::string solver_its = "GMRES_its";
    std::string solver_res = "GMRES_res";

    // need to add a first elastic prediction

    FLAG1=1;
    while (FLAG1==1){

        FLAG1=0;
        FLAG2=1;
        FLAG3=1;
        nb_bis = 0;
        time += delta;
        y = x;

        while (FLAG2==1){
            FLAG2=0;
            if(time-time_max>eps){
                if(time_max+delta-time>eps){
                    delta=time_max+delta-time;
                    time=time_max;
                }
                else{
                    FLAG1=0;
                    break;
                }
            }
            //interface->pressure_load = time*pressure_load;

            iter = 0;
            while (FLAG3==1){

                FLAG3=0;
                iter++;

                if(iter>iter_max){
                    if (MyPID==0) std::cout << "Iteration limit exceeds.\n";
                    return 1;
                }

                Time.ResetStartTime();
                get_matrix_and_rhs(eto, ep_old, stiffness, rhs);
                Assemble_time = Time.ElapsedTime();

                displacement = (iter==1) ? (delta*bc_disp) : (0.0);

                apply_dirichlet_conditions(stiffness, rhs, displacement);

                if(iter>1){

                    rhs.NormInf(&norm_inf_rhs);

                    if (MyPID==0 && print){
                        if(iter>2){
                            std::cout << "\t\t\t" << iter << "\t" << norm_inf_rhs << "\t" << Krylov_its << "\t\t" << Krylov_res << "\t\t" << Assemble_time << "\t\t\t" << Aztec_time << "\n"; /* << "\t\t\t" << time*pressure_load << "\n"; */
                        }
                        else{
                            std::cout << "\n Time" << "\t" << "Timestep" << "\t" << "Iter" << "\t" << "NormInf" << "\t\t" << solver_its << "\t" << solver_res << "\t\t" << "assemble_time" << "\t\t" << "aztec_time" << "\t\t" << "\n"; /*"pressure_load [kPa]" << "\n";*/
                            std::cout << " " << time << "\t" << delta << "\t\t" << iter << "\t" << norm_inf_rhs << "\t" << Krylov_its << "\t\t" << Krylov_res << "\t\t" << Assemble_time << "\t\t\t" << Aztec_time << "\t\t\t" << "\n"; /*time*pressure_load << "\n";*/
                        }
                    }

                    if (norm_inf_rhs<norm_inf_tol){
                        if (iter<=iter_min) delta = success*delta;
                        if (iter>=iter_max) delta = delta/failure;
                        FLAG1=1;
                        break;
                    }

                    if (norm_inf_rhs>norm_inf_max||iter==iter_max){
                        nb_bis++;
                        if (nb_bis<nb_bis_max){
                            delta /= failure;
                            time -= delta;
                            x = y;
                            if (MyPID==0) std::cout << "Bisecting increment: " << nb_bis << "\n";
                        }
                        else{
                            if (MyPID==0) std::cout << "Bisection number exceeds.\n";
                            return 2;
                        }
                        FLAG2=1;
                        FLAG3=1;
                        break;
                    }

                }

                lhs.PutScalar(0.0);
                problem.SetOperator(&stiffness);
                problem.SetLHS(&lhs);
                problem.SetRHS(&rhs);
                solver.SetProblem(problem);
                solver.SetParameters(*Krylov);

                Time.ResetStartTime();
                solver.Iterate(2000,tol);
                Aztec_time = Time.ElapsedTime();
                Krylov_its = solver.NumIters();
                Krylov_res = solver.TrueResidual();
                x.Update(1.0,lhs,1.0);

                // compute eto here ?

                FLAG3=1;
            }
        }
    // compute ep here ?
    }
    return 0;
}

void plasticitySmallStrains::get_matrix_and_rhs(Epetra_Vector & eto, Epetra_Vector & ep_old, Epetra_FECrsMatrix & K, Epetra_FEVector & F){
  
}

void plasticitySmallStrains::constructGaussMap(Epetra_Map & GaussMap){
  int e_gid;
  int n_local_cells = Mesh->n_local_cells;
  int n_gauss_cells = Mesh->n_gauss_cells;
  std::vector<int> local_gauss_points(n_local_cells*n_gauss_cells);
  for (unsigned int e_lid=0; e_lid<n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];
      for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
          local_gauss_points[e_lid*n_gauss_cells+gp] = e_gid*n_gauss_cells+gp;
      }
  }
  GaussMap = Epetra_Map(-1, n_local_cells*n_gauss_cells, &local_gauss_points[0], 0, *Comm);
}

void plasticitySmallStrains::create_FECrsGraph(){
    FEGraph = new Epetra_FECrsGraph(Copy,*StandardMap,100);
    int eglob, node;
    int *index;
    index = new int [3*Mesh->el_type];

    for (int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
        eglob = Mesh->local_cells[e_lid];
        for (int inode=0; inode<Mesh->el_type; ++inode){
            node = Mesh->cells_nodes[Mesh->el_type*eglob+inode];
            for (int ddl=0; ddl<3; ++ddl){
                index[3*inode+ddl] = 3*node+ddl;
            }
        }
        for (int i=0; i<3*Mesh->el_type; ++i){
            for (int j=0; j<3*Mesh->el_type; ++j){
                FEGraph->InsertGlobalIndices(1, &index[i], 1, &index[j]);
            }
        }
    }
    Comm->Barrier();
    FEGraph->GlobalAssemble();
    delete[] index;
}
