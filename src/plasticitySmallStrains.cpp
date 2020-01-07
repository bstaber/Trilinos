/*
Brian Staber (brian.staber@gmail.com)
*/

#include "plasticitySmallStrains.hpp"

template <typename T>
std::vector<T> linspace(T a, T b, size_t N);

template<typename T>
std::ostream &operator <<(std::ostream &os, const std::vector<T> &v);

plasticitySmallStrains::plasticitySmallStrains(){

}

plasticitySmallStrains::~plasticitySmallStrains(){

}

void plasticitySmallStrains::initialize(Epetra_Comm & comm, Teuchos::ParameterList & parameterlist){

  Delta         = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "delta");
  iter_min      = Teuchos::getParameter<int>(parameterlist.sublist("Newton"),    "iterMin");
  iter_max      = Teuchos::getParameter<int>(parameterlist.sublist("Newton"),    "iterMax");
  nb_bis_max    = Teuchos::getParameter<int>(parameterlist.sublist("Newton"),    "nbBisMax");
  norm_inf_tol  = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "NormFTol");
  norm_inf_max  = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "NormFMax");
  eps           = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "eps");
  success       = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "success_parameter");
  failure       = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "failure_parameter");
  tol           = Teuchos::getParameter<double>(parameterlist.sublist("Newton"), "tol");

  bc_disp       = Teuchos::getParameter<double>(parameterlist.sublist("BoundaryConditions"), "bc_disp");

  Direct = &parameterlist.sublist("Direct");
  Krylov = &parameterlist.sublist("Krylov");
  Newton = &parameterlist.sublist("Newton");
  BCs    = &parameterlist.sublist("BoundaryConditions");

  //Direct->set("MaxProcs",Comm->NumProc());

  CpuTime = new Epetra_Time(*Comm);

  //std::string mesh_file = Teuchos::getParameter<std::string>(parameterlist.sublist("Mesh"), "mesh_file");
  //double scaling = Teuchos::getParameter<double>(parameterlist.sublist("Mesh"), "scaling");
  Mesh  = new mesh(comm, parameterlist);
  Comm  = Mesh->Comm;
  MyPID = Comm->MyPID();

  StandardMap           = new Epetra_Map(-1,3*Mesh->n_local_nodes_without_ghosts,&Mesh->local_dof_without_ghosts[0],0,*Comm);
  OverlapMap            = new Epetra_Map(-1,3*Mesh->n_local_nodes,&Mesh->local_dof[0],0,*Comm);
  ImportToOverlapMap    = new Epetra_Import(*OverlapMap,*StandardMap);

  create_FECrsGraph();
  constructScalarGaussMap();
  setup_dirichlet_conditions();

  epcum = new Epetra_Vector(*GaussMap,true);
  sig = new Epetra_MultiVector(*GaussMap,6,true);
  eto = new Epetra_MultiVector(*GaussMap,6,true);

  epcum_converged = new Epetra_Vector(*GaussMap,true);
  sig_converged   = new Epetra_MultiVector(*GaussMap,6,true);
  eto_converged   = new Epetra_MultiVector(*GaussMap,6,true);

  tgm = new Epetra_MultiVector(*GaussMap,21,true);

  ELASTICITY.Reshape(6,6);
}

int plasticitySmallStrains::sequence_bvp(bool print){
  double delta = Delta;
  double displacement;
  double norm_inf_rhs;
  double time_max = 1.0;
  double Krylov_res = 0;
  int iter, Krylov_its;
  std::string solver_its = "GMRES_its";
  std::string solver_res = "GMRES_res";

  double nRes, nRes0;

  int Nincr = Teuchos::getParameter<int>(*Newton,"Nincr");
  double norm_l2_tol = Teuchos::getParameter<double>(*Newton,"norm_l2_tol");

  Epetra_LinearProblem problem;

  //this part is not useful yet, always using AztecOO
  std::string linear_solver;
  if (Newton->isParameter("linear_solver")) linear_solver = Teuchos::getParameter<std::string>(*Newton,"linear_solver");
  else linear_solver = "sparse_iterative";

  AztecOO solver;
  Amesos  solver_amesos;

  Epetra_FECrsMatrix A(Copy,*FEGraph);
  Epetra_FEVector Res(*StandardMap,true);

  Epetra_Vector u(*StandardMap,true);
  Epetra_Vector du(*StandardMap,true);
  Epetra_Vector Du(*StandardMap,true);
  Epetra_Vector DuOverlaped(*OverlapMap,true);

  std::vector<double> uload = linspace<double>(0.0,bc_disp,Nincr);
  uload.erase(uload.begin());
  integrate_constitutive_problem(DuOverlaped);

  for (unsigned int i=0; i<uload.size(); ++i){
    assemble_system(DuOverlaped,A,Res);
    if (i==0) apply_dirichlet_conditions(A,Res,uload[i]);
    else apply_dirichlet_conditions(A,Res,uload[i]-uload[i-1]);
    Res.Norm2(&nRes0);
    nRes = nRes0;
    Du.PutScalar(0.0);
    if (MyPID==0) std::cout << "Increment:" << i+1 << std::setw(16) << "Residual" << std::endl;
    if (MyPID==0) std::cout << std::setw(24) << nRes0 << std::endl;
    iter = 0;
    while ((nRes/nRes0>norm_l2_tol) & (iter<iter_max)) {
      //solve(problem, solver, A, Res, du);
      //AztecOO_solver(A, Res, du);
      //Amesos_solver(A, Res, du);
      if (linear_solver=="sparse_iterative") {
        AztecOO_solver(A, Res, du);
      }
      else if (linear_solver=="sparse_direct") {
        Amesos_solver(A, Res, du);
      }
      else {
        Amesos_solver(A, Res, du);
      }
      Du.Update(1.0,du,1.0);
      DuOverlaped.Import(Du,*ImportToOverlapMap,Insert);
      integrate_constitutive_problem(DuOverlaped);
      assemble_system(DuOverlaped,A,Res);
      apply_dirichlet_conditions(A,Res,0.0);
      Res.Norm2(&nRes);
      if (MyPID==0) std::cout << std::setw(24) << nRes << std::endl;
      iter++;
    }
    u.Update(1.0,Du,1.0);
    (*sig_converged)   = (*sig);
    (*epcum_converged) = (*epcum);
  }
  print_solution(u,"/Users/brian/Documents/GitHub/TrilinosUQComp/results/plasticity/plate/plate_u.mtx");
  return 0;
}

int plasticitySmallStrains::incremental_bvp(bool print){

    Epetra_LinearProblem problem;
    AztecOO solver;

    Epetra_FECrsMatrix A(Copy,*FEGraph);
    Epetra_FEVector Res(*StandardMap,true);

    Epetra_Vector u(*StandardMap,true);
    Epetra_Vector u_converged(*StandardMap,true);
    Epetra_Vector du(*StandardMap,true);
    Epetra_Vector Du(*StandardMap,true);
    Epetra_Vector DuOverlaped(*OverlapMap,true);

    double delta = Delta;
    double norm_inf_rhs;
    double time_max = 1.0;
    int FLAG1, FLAG2, FLAG3, nb_bis, iter;
    double Krylov_res = 0;
    int Krylov_its = 0;
    std::string solver_its = "GMRES_its";
    std::string solver_res = "GMRES_res";

    std::string linear_solver;
    if (Newton->isParameter("linear_solver")) linear_solver = Teuchos::getParameter<std::string>(*Newton,"linear_solver");
    else linear_solver = "sparse_iterative";

    time = 0.0;
    integrate_constitutive_problem(DuOverlaped);

    FLAG1=1;
    while (FLAG1==1){

        FLAG1=0;
        FLAG2=1;
        FLAG3=1;
        nb_bis = 0;
        time += delta;

        u.Update(1.0,Du,1.0);
        u_converged = u;
        (*sig_converged)   = (*sig);
        (*epcum_converged) = (*epcum);

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

            assemble_system(DuOverlaped,A,Res);
            apply_dirichlet_conditions(A,Res,bc_disp*delta);
            Du.PutScalar(0.0);
            iter = 0;
            norm_inf_rhs = 1.0;
            if (MyPID==0) std::cout << "\n Time" << "\t" << "Timestep" << "\t" << "Iter" << "\t" << "NormInf" << std::endl;
            if (MyPID==0) std::cout << " " << time << "\t" << delta << "\t\t" << iter << std::endl;

            while (FLAG3==1){
                //solve(problem, solver, A, Res, du);
                //AztecOO_solver(A, Res, du);
                //Amesos_solver(A, Res, du);
                if (linear_solver=="sparse_iterative") {
                  AztecOO_solver(A, Res, du);
                }
                else if (linear_solver=="sparse_direct") {
                  Amesos_solver(A, Res, du);
                }
                else {
                  Amesos_solver(A, Res, du);
                }
                Du.Update(1.0,du,1.0);
                DuOverlaped.Import(Du,*ImportToOverlapMap,Insert);
                integrate_constitutive_problem(DuOverlaped);
                assemble_system(DuOverlaped,A,Res);
                apply_dirichlet_conditions(A,Res,0.0);

                FLAG3=0;
                iter++;
                if(iter>iter_max){
                    if (MyPID==0) std::cout << "Iteration limit exceeds.\n";
                    return 1;
                }

                Res.NormInf(&norm_inf_rhs);
                if (MyPID==0) std::cout << " " << time << "\t" << delta << "\t\t" << iter << "\t" << norm_inf_rhs << std::endl;

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
                        u  = u_converged;
                        (*sig)   = (*sig_converged);
                        (*epcum) = (*epcum_converged);
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
                FLAG3=1;
            }
        }
    }
    print_solution(u,"/Users/brian/Documents/GitHub/TrilinosUQComp/results/plasticity/plate/plate_u.mtx");
    return 0;
}

void plasticitySmallStrains::solve(Epetra_LinearProblem & problem_, AztecOO & solver_,
                                   Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & lhs_){

  lhs_.PutScalar(0.0);

  problem_.SetOperator(&A);
  problem_.SetLHS(&lhs_);
  problem_.SetRHS(&b);
  solver_.SetProblem(problem_);
  solver_.SetParameters(*Krylov);

  CpuTime->ResetStartTime();
  solver_.Iterate(2000,tol);
  Aztec_time = CpuTime->ElapsedTime();
}

void plasticitySmallStrains::AztecOO_solver(Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & xx) {

  Epetra_LinearProblem problem;
  AztecOO solver;

  xx.PutScalar(0.0);

  problem.SetOperator(&A);
  problem.SetLHS(&xx);
  problem.SetRHS(&b);
  solver.SetProblem(problem);
  solver.SetParameters(*Krylov);

  CpuTime->ResetStartTime();
  solver.Iterate(2000,tol);
  Aztec_time = CpuTime->ElapsedTime();
}

void plasticitySmallStrains::Amesos_solver(Epetra_FECrsMatrix & A, Epetra_FEVector & b, Epetra_Vector & xx) {

  std::string amesos_type = Teuchos::getParameter<std::string>(*Direct, "solver_type");
  bool isAvailable;

  Amesos_BaseSolver * Solver;
  Amesos Factory;
  /*isAvailable = Factory.Query("Lapack"); if(isAvailable){ std::cout << "Lapack: yes.\n"; } else{ std::cout << "Lapack: no.\n"; }
  isAvailable = Factory.Query("Klu"); if(isAvailable){ std::cout << "Klu: yes.\n"; } else{ std::cout << "Klu: no.\n"; }
  isAvailable = Factory.Query("Umfpack"); if(isAvailable){ std::cout << "Umfpack: yes.\n"; } else{ std::cout << "Umfpack: no.\n"; }
  isAvailable = Factory.Query("Pardiso"); if(isAvailable){ std::cout << "Pardiso: yes.\n"; } else{ std::cout << "Pardiso: no.\n"; }
  isAvailable = Factory.Query("Taucs"); if(isAvailable){ std::cout << "Taucs: yes.\n"; } else{ std::cout << "Taucs: no.\n"; }
  isAvailable = Factory.Query("Superlu"); if(isAvailable){ std::cout << "Superlu: yes.\n"; } else{ std::cout << "Superlu: no.\n"; }
  isAvailable = Factory.Query("Superludist"); if(isAvailable){ std::cout << "Superludist: yes.\n"; } else{ std::cout << "Superludist: no.\n"; }
  isAvailable = Factory.Query("Mumps"); if(isAvailable){ std::cout << "Mumps: yes.\n"; } else{ std::cout << "Mumps: no.\n"; }
  isAvailable = Factory.Query("Dscpack"); if(isAvailable){ std::cout << "Dscpack: yes.\n"; } else{ std::cout << "Dscpack: no.\n"; }*/
  isAvailable = Factory.Query(amesos_type);
  if (!isAvailable) {
    if (MyPID==0) std::cout << "Amesos solver type not avaiable" << std::endl;
  }
  /*else {
    if (MyPID==0) std::cout << "Amesos " << amesos_type << " solver avaiable" << std::endl;
  }*/

  xx.PutScalar(0.0);

  Epetra_LinearProblem problem;
  problem.SetOperator(&A);
  problem.SetLHS(&xx);
  problem.SetRHS(&b);

  Solver = Factory.Create(amesos_type, problem);
  Solver->SetParameters(*Direct);
  Solver->SymbolicFactorization();
  Solver->NumericFactorization();
  Solver->Solve();

  delete Solver;
}

void plasticitySmallStrains::assemble_system(const Epetra_Vector & Du_, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

  F.PutScalar(0.0);
  K.PutScalar(0.0);

  stiffness_rhs_homogeneousForcing(Du_, K, F);

  Comm->Barrier();

  K.GlobalAssemble();
  K.FillComplete();
  F.GlobalAssemble();
}

void plasticitySmallStrains::stiffness_rhs_homogeneousForcing(const Epetra_Vector & Du_, Epetra_FECrsMatrix & K, Epetra_FEVector & F){

  int node, e_gid, error;
  int n_gauss_points = Mesh->n_gauss_cells;
  double gauss_weight;

  int *Indices_cells;
  Indices_cells = new int [3*Mesh->el_type];

  Epetra_SerialDenseVector du_el(3*Mesh->el_type);
  Epetra_SerialDenseVector sig_el(6);
  Epetra_SerialDenseMatrix m_tg_matrix(6,6);

  Epetra_SerialDenseVector Re(3*Mesh->el_type);
  Epetra_SerialDenseMatrix Ke(3*Mesh->el_type,3*Mesh->el_type);

  Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
  Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
  Epetra_SerialDenseMatrix B_times_TM(3*Mesh->el_type,6);

  for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];

      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
          for (int iddl=0; iddl<3; ++iddl){
              Indices_cells[3*inode+iddl] = 3*node+iddl;
              Re(3*inode+iddl) = 0.0;
              du_el(3*inode+iddl) = Du_[OverlapMap->LID(3*node+iddl)];
              for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                  for (int jddl=0; jddl<3; ++jddl){
                      Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
                  }
              }
          }
      }

      for (unsigned int gp=0; gp<n_gauss_points; ++gp){
          gauss_weight = Mesh->gauss_weight_cells(gp);
          for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
              dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
              dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
              dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);
          }
          compute_B_matrices(dx_shape_functions,matrix_B);

          for (unsigned int k=0; k<6; ++k) {
            sig_el(k) = (*sig)[k][GaussMap->LID(int(e_gid*n_gauss_points+gp))];
            for (unsigned int l=k; l<6; ++l) {
              m_tg_matrix(k,l) = (*tgm)[15-(6-k)*(6-k-1)/2+l][GaussMap->LID(int(e_gid*n_gauss_points+gp))];
              m_tg_matrix(l,k) = m_tg_matrix(k,l); // find a better way lol
            }
          }

          //if (MyPID==0) std::cout << "m_tg_matrix = " << m_tg_matrix << std::endl;

          error = Re.Multiply('T','N',-gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,sig_el,1.0);
          error = B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,m_tg_matrix,0.0);
          error = Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
      }

      for (unsigned int i=0; i<3*Mesh->el_type; ++i){
          error = F.SumIntoGlobalValues(1, &Indices_cells[i], &Re(i));
          for (unsigned int j=0; j<3*Mesh->el_type; ++j){
              error = K.SumIntoGlobalValues(1, &Indices_cells[i], 1, &Indices_cells[j], &Ke(i,j));
          }
      }
  }
  delete[] Indices_cells;
}

void plasticitySmallStrains::integrate_constitutive_problem(Epetra_Vector & Du_){

  int e_gid, node;

  int n_gauss_points = Mesh->n_gauss_cells;
  Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
  Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);

  double epcum_el;
  Epetra_SerialDenseVector du_el(3*Mesh->el_type);
  Epetra_SerialDenseVector sig_el(6);
  Epetra_SerialDenseVector deto_el(6);
  Epetra_SerialDenseMatrix m_tg_matrix(6,6);

  for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];
      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
          for (int iddl=0; iddl<3; ++iddl){
              du_el(3*inode+iddl) = Du_[OverlapMap->LID(3*node+iddl)];
          }
      }

    for (unsigned int gp=0; gp<n_gauss_points; ++gp){
      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
          dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
          dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);
      }

      compute_B_matrices(dx_shape_functions,matrix_B);
      deto_el.Multiply('N','N',1.0,matrix_B,du_el,0.0);

      for (unsigned int k=0; k<6; ++k) {
        sig_el(k) = (*sig_converged)[k][GaussMap->LID(int(e_gid*n_gauss_points+gp))];
      }
      epcum_el = (*epcum_converged)[GaussMap->LID(int(e_gid*n_gauss_points+gp))];

      constitutive_problem(e_lid, gp, deto_el, sig_el, epcum_el, m_tg_matrix);

      /*std::cout << "deto_el " << deto_el << std::endl;
      std::cout << "sig_el " << sig_el << std::endl;
      std::cout << "m_tg_matrix " << m_tg_matrix << std::endl;
      std::cout << "epcum_el = " << epcum_el << std::endl;*/

      for (unsigned int k=0; k<6; ++k) {
        (*sig)[k][GaussMap->LID(int(e_gid*n_gauss_points+gp))] = sig_el(k);
        (*eto)[k][GaussMap->LID(int(e_gid*n_gauss_points+gp))] += deto_el(k);
        for (unsigned int l=k; l<6; ++l) {
          (*tgm)[15-(6-k)*(6-k-1)/2+l][GaussMap->LID(int(e_gid*n_gauss_points+gp))] = m_tg_matrix(k,l);
        }
      }
      (*epcum)[GaussMap->LID(int(e_gid*n_gauss_points+gp))] = epcum_el;
    }
  }
}

void plasticitySmallStrains::elastic_predictor(Epetra_LinearProblem & problem_, AztecOO & solver_, Epetra_FECrsMatrix & K,
                                               Epetra_FEVector & rhs_, Epetra_Vector & lhs_, double & displacement_) {

  rhs_.PutScalar(0.0);
  CpuTime->ResetStartTime();
  assemble_system_LinearElasticity(K);
  Assemble_time = CpuTime->ElapsedTime();
  apply_dirichlet_conditions(K, rhs_, displacement_);

  lhs_.PutScalar(0.0);
  problem_.SetOperator(&K);
  problem_.SetLHS(&lhs_);
  problem_.SetRHS(&rhs_);
  solver_.SetProblem(problem_);
  solver_.SetParameters(*Krylov);

  CpuTime->ResetStartTime();
  solver_.Iterate(2000,tol);
  Aztec_time = CpuTime->ElapsedTime();
}

void plasticitySmallStrains::assemble_system_LinearElasticity(Epetra_FECrsMatrix & K){

  K.PutScalar(0.0);
  stiffness_homogeneousForcing_LinearElasticity(K);

  Comm->Barrier();
  int error;
  error=K.GlobalAssemble();
  error=K.FillComplete();
}

void plasticitySmallStrains::stiffness_homogeneousForcing_LinearElasticity(Epetra_FECrsMatrix & K){

  int node, e_gid, error;
  int n_gauss_points = Mesh->n_gauss_cells;
  double gauss_weight;

  int *Indices_cells;
  Indices_cells = new int [3*Mesh->el_type];

  Epetra_SerialDenseMatrix Ke(3*Mesh->el_type,3*Mesh->el_type);

  Epetra_SerialDenseMatrix dx_shape_functions(Mesh->el_type,3);
  Epetra_SerialDenseMatrix matrix_B(6,3*Mesh->el_type);
  Epetra_SerialDenseMatrix B_times_TM(3*Mesh->el_type,6);

  for (unsigned int e_lid=0; e_lid<Mesh->n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];

      for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
          for (int iddl=0; iddl<3; ++iddl){
              Indices_cells[3*inode+iddl] = 3*node+iddl;
              for (unsigned int jnode=0; jnode<Mesh->el_type; ++jnode){
                  for (int jddl=0; jddl<3; ++jddl){
                      Ke(3*inode+iddl,3*jnode+jddl) = 0.0;
                  }
              }
          }
      }

      for (unsigned int gp=0; gp<n_gauss_points; ++gp){
          gauss_weight = Mesh->gauss_weight_cells(gp);
          for (unsigned int inode=0; inode<Mesh->el_type; ++inode){
              dx_shape_functions(inode,0) = Mesh->DX_N_cells(gp+n_gauss_points*inode,e_lid);
              dx_shape_functions(inode,1) = Mesh->DY_N_cells(gp+n_gauss_points*inode,e_lid);
              dx_shape_functions(inode,2) = Mesh->DZ_N_cells(gp+n_gauss_points*inode,e_lid);
          }

          compute_B_matrices(dx_shape_functions,matrix_B);
          get_elasticity_tensor(e_lid, gp, ELASTICITY);

          error = B_times_TM.Multiply('T','N',gauss_weight*Mesh->detJac_cells(e_lid,gp),matrix_B,ELASTICITY,0.0);
          error = Ke.Multiply('N','N',1.0,B_times_TM,matrix_B,1.0);
      }

      for (unsigned int i=0; i<3*Mesh->el_type; ++i){
          for (unsigned int j=0; j<3*Mesh->el_type; ++j){
              error = K.SumIntoGlobalValues(1, &Indices_cells[i], 1, &Indices_cells[j], &Ke(i,j));
          }
      }
  }
  delete[] Indices_cells;
}

void plasticitySmallStrains::compute_B_matrices(Epetra_SerialDenseMatrix & dx_shape_functions, Epetra_SerialDenseMatrix & B){
    double factor = 1.0/std::sqrt(2.0);
    for (unsigned inode=0; inode<Mesh->el_type; ++inode){
        B(0,3*inode) = dx_shape_functions(inode,0);
        B(0,3*inode+1) = 0.0;
        B(0,3*inode+2) = 0.0;

        B(1,3*inode) = 0.0;
        B(1,3*inode+1) = dx_shape_functions(inode,1);
        B(1,3*inode+2) = 0.0;

        B(2,3*inode) = 0.0;
        B(2,3*inode+1) = 0.0;
        B(2,3*inode+2) = dx_shape_functions(inode,2);

        B(3,3*inode) = 0.0;
        B(3,3*inode+1) = factor*dx_shape_functions(inode,2);
        B(3,3*inode+2) = factor*dx_shape_functions(inode,1);

        B(4,3*inode) = factor*dx_shape_functions(inode,2);
        B(4,3*inode+1) = 0.0;
        B(4,3*inode+2) = factor*dx_shape_functions(inode,0);

        B(5,3*inode) = factor*dx_shape_functions(inode,1);
        B(5,3*inode+1) = factor*dx_shape_functions(inode,0);
        B(5,3*inode+2) = 0.0;
    }
}


void plasticitySmallStrains::constructScalarGaussMap(){
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
  Comm->Barrier();
  GaussMap = new Epetra_Map(-1, n_local_cells*n_gauss_cells, &local_gauss_points[0], 0, *Comm);
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

int plasticitySmallStrains::print_solution(Epetra_Vector & solution, std::string fileName){

    int NumTargetElements = 0;
    if (Comm->MyPID()==0){
        NumTargetElements = 3*Mesh->n_nodes;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(*StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(solution,ExportOnRoot,Insert);

    int error = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(),lhs_root,0,0,false);

    return error;
}

int plasticitySmallStrains::print_at_gauss_points(std::string filename){

  /*Epetra_Vector x_coord(*GaussMap);
  Epetra_Vector y_coord(*GaussMap);
  Epetra_Vector z_coord(*GaussMap);*/
  Epetra_MultiVector output(*GaussMap,4,true);

  int node;
  double xi, eta, zeta;
  Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type);
  Epetra_SerialDenseVector vector_X(3);
  Epetra_SerialDenseVector N(Mesh->el_type);

  int e_gid;
  int n_local_cells = Mesh->n_local_cells;
  int n_gauss_cells = Mesh->n_gauss_cells;

  for (unsigned int e_lid=0; e_lid<n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];
      for (int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
          matrix_X(0,inode) = Mesh->nodes_coord[3*node+0];
          matrix_X(1,inode) = Mesh->nodes_coord[3*node+1];
          matrix_X(2,inode) = Mesh->nodes_coord[3*node+2];
      }
      for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
          xi   = Mesh->xi_cells[gp];
          eta  = Mesh->eta_cells[gp];
          zeta = Mesh->zeta_cells[gp];

          switch (Mesh->el_type){
              case 4:
                  tetra4::shape_functions(N, xi, eta, zeta);
                  break;
              case 8:
                  hexa8::shape_functions(N, xi, eta, zeta);
                  break;
              case 10:
                  tetra10::shape_functions(N, xi, eta, zeta);
                  break;
              case 20:
                  hexa20::shape_functions(N, xi, eta, zeta);
                  break;
              default:
                  std::cout << "Non supported element type encountered" << std::endl;
                  //NotImplementedException();
          }

          vector_X.Multiply('N','N',1.0,matrix_X,N,0.0);
          output[0][GaussMap->LID(int(e_gid*n_gauss_cells+gp))] = vector_X(0);
          output[1][GaussMap->LID(int(e_gid*n_gauss_cells+gp))] = vector_X(1);
          output[2][GaussMap->LID(int(e_gid*n_gauss_cells+gp))] = vector_X(2);
          output[3][GaussMap->LID(int(e_gid*n_gauss_cells+gp))] = (*epcum)[GaussMap->LID(int(e_gid*n_gauss_cells+gp))];
       }
  }

  int error;
  int NumTargetElements = 0;
  if (Comm->MyPID()==0){
      NumTargetElements = Mesh->n_cells*n_gauss_cells;
  }
  Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
  Epetra_Export ExportOnRoot(*GaussMap,MapOnRoot);
  Epetra_MultiVector lhs_root(MapOnRoot,4,true);

  lhs_root.Export(output,ExportOnRoot,Insert);
  error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
  return 0;
}

int plasticitySmallStrains::print_mean_at_center(std::string filename){

  Epetra_Map CellsMap(-1,Mesh->n_local_cells,&Mesh->local_cells[0],0,*Comm);
  Epetra_MultiVector output(CellsMap,4,true);

  int node;
  double xi   = 0.0;
  double eta  = 0.0;
  double zeta = 0.0;
  Epetra_SerialDenseMatrix matrix_X(3,Mesh->el_type);
  Epetra_SerialDenseVector vector_X(3);
  Epetra_SerialDenseVector N(Mesh->el_type);

  int e_gid;
  int n_local_cells = Mesh->n_local_cells;
  int n_gauss_cells = Mesh->n_gauss_cells;
  double vol_cell;

  switch (Mesh->el_type){
      case 8:
          hexa8::shape_functions(N, xi, eta, zeta);
          break;
      case 20:
          hexa20::shape_functions(N, xi, eta, zeta);
          break;
      default:
          std::cout << "Non supported element type encountered in plasticitySmallStrains::print_mean_at_center" << std::endl;
          //NotImplementedException();
  }

  for (unsigned int e_lid=0; e_lid<n_local_cells; ++e_lid){
      e_gid = Mesh->local_cells[e_lid];
      vol_cell = Mesh->vol_cells(e_lid);
      for (int inode=0; inode<Mesh->el_type; ++inode){
          node = Mesh->cells_nodes[Mesh->el_type*e_gid+inode];
          matrix_X(0,inode) = Mesh->nodes_coord[3*node+0];
          matrix_X(1,inode) = Mesh->nodes_coord[3*node+1];
          matrix_X(2,inode) = Mesh->nodes_coord[3*node+2];
      }

      vector_X.Multiply('N','N',1.0,matrix_X,N,0.0);
      output[0][CellsMap.LID(e_gid)] = vector_X(0);
      output[1][CellsMap.LID(e_gid)] = vector_X(1);
      output[2][CellsMap.LID(e_gid)] = vector_X(2);

      for (unsigned int gp=0; gp<n_gauss_cells; ++gp){
          output[3][CellsMap.LID(e_gid)] += (*epcum)[GaussMap->LID(int(e_gid*n_gauss_cells+gp))]*Mesh->gauss_weight_cells(gp)*Mesh->detJac_cells(e_lid,gp);
      }
      output[3][CellsMap.LID(e_gid)] /= vol_cell;

  }

  int error;
  int NumTargetElements = 0;
  if (Comm->MyPID()==0){
      NumTargetElements = Mesh->n_cells;
  }
  Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
  Epetra_Export ExportOnRoot(CellsMap,MapOnRoot);
  Epetra_MultiVector lhs_root(MapOnRoot,4,true);

  lhs_root.Export(output,ExportOnRoot,Insert);
  error = EpetraExt::MultiVectorToMatrixMarketFile(filename.c_str(),lhs_root,0,0,false);
  return 0;
}

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N-1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
        *x = val;
    return xs;
}

template<typename T>
std::ostream &operator <<(std::ostream &os, const std::vector<T> &v) {
   using namespace std;
   copy(v.begin(), v.end(), ostream_iterator<T>(os, "\t"));
   return os;
}
