#include "Newton_Raphsonpp.hpp"

Newton_Raphson::Newton_Raphson(Finite_Element_Problem & Interface, Teuchos::ParameterList & Parameters){
    
    interface = &Interface;
    Comm = interface->Comm;
    Krylov = &Parameters.sublist("Krylov");
    MyPID = Comm->MyPID();
    
    setParameters(Parameters);
    
    x = new Epetra_Vector(*interface->StandardMap);
}

void Newton_Raphson::setParameters(Teuchos::ParameterList & Parameters){

    delta = Teuchos::getParameter<double>(Parameters.sublist("Newton"), "delta");
    iter_min = Teuchos::getParameter<int>(Parameters.sublist("Newton"), "iterMin");
    iter_max = Teuchos::getParameter<int>(Parameters.sublist("Newton"), "iterMax");
    nb_bis_max = Teuchos::getParameter<int>(Parameters.sublist("Newton"), "nbBisMax");
    norm_inf_tol = Teuchos::getParameter<double>(Parameters.sublist("Newton"), "NormFTol");
    norm_inf_max = Teuchos::getParameter<double>(Parameters.sublist("Newton"), "NormFMax");
    eps = Teuchos::getParameter<double>(Parameters.sublist("Newton"), "eps");
    success = Teuchos::getParameter<double>(Parameters.sublist("Newton"), "success_parameter");
    failure = Teuchos::getParameter<double>(Parameters.sublist("Newton"), "failure_parameter");
    bc_disp = Teuchos::getParameter<double>(Parameters.sublist("Newton"), "bc_disp");
    pressure_load = Teuchos::getParameter<double>(Parameters.sublist("Newton"), "pressure_load");
        
}

int Newton_Raphson::Solve_with_Stratimikos(Teuchos::RCP<Teuchos::ParameterList> solverBuilderSL){
    
    double assemble_time,displacement, norm_inf_rhs, time_init, time_max, krylov_res;
    int FLAG1, FLAG2, FLAG3, nb_bis, iter, krylov_its;
    time_init = delta;
    time_max = 1.0;
    time = 0.0;
    
    std::string solver_its = "krylov_its";
    std::string solver_res = "krylov_res";
    
    Epetra_Time Time(*Comm);
    
    Teuchos::RCP<Epetra_FECrsMatrix> epetra_A = Teuchos::rcp(new Epetra_FECrsMatrix(Copy,*interface->FEGraph));
    Teuchos::RCP<Epetra_FEVector> epetra_rhs = Teuchos::rcp(new Epetra_FEVector(*interface->StandardMap));
    Teuchos::RCP<Epetra_Vector> epetra_lhs = Teuchos::rcp(new Epetra_Vector(*interface->StandardMap));
    Teuchos::RCP<Epetra_Vector> epetra_y = Teuchos::rcp(new Epetra_Vector(*interface->StandardMap));
    Teuchos::RCP<Epetra_MultiVector> epetra_rhsv = epetra_rhs;
    
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
    linearSolverBuilder.setParameterList(solverBuilderSL);
    
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory =
    linearSolverBuilder.createLinearSolveStrategy("");
    Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > lows = lowsFactory->createOp();
    
    Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
    //lowsFactory->setOStream(out);
    //lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);
    
    Teuchos::RCP<const Thyra::LinearOpBase<double>> thyra_A = Thyra::epetraLinearOp(epetra_A);
    Teuchos::RCP<Thyra::VectorBase<double>> thyra_lhs = Thyra::create_Vector(epetra_lhs,thyra_A->domain());
    Teuchos::RCP<const Thyra::MultiVectorBase<double>> thyra_rhs = Thyra::create_MultiVector(epetra_rhsv,thyra_A->range());
    
    Thyra::SolveStatus<double> status;
        
    FLAG1=1;
    while (FLAG1==1){
        
        FLAG1=0;
        FLAG2=1;
        FLAG3=1;
        nb_bis = 0;
        time += delta;
        *epetra_y = *x;
        
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
            
            interface->pressure_load = time*pressure_load;
            iter = 0;
            while (FLAG3==1){
                FLAG3=0;
                iter++;
                if(iter>iter_max){
                    if (MyPID==0){std::cout << "Iteration limit exceeds.\n";}
                    return 1;
                }
                Time.ResetStartTime();
                interface->get_matrix_and_rhs(*x, *epetra_A, *epetra_rhs);
                assemble_time = Time.ElapsedTime();
                
                if (iter==1){
                    displacement = delta*bc_disp;
                }
                else{
                    displacement = 0.0;
                }
                interface->apply_dirichlet_conditions(*epetra_A, *epetra_rhs, displacement);
                
                if(iter>1){
                    epetra_rhs->NormInf(&norm_inf_rhs);
                    
                    if (MyPID==0){
                        if(iter>2){
                            std::cout << "\t\t\t" << iter << "\t" << norm_inf_rhs << "\t" << krylov_res << "\t\t" << assemble_time << "\t\t\t" << time*pressure_load/1000.0 << "\n";
                        }
                        else{
                            std::cout << "\n Time" << "\t" << "Timestep" << "\t" << "Iter" << "\t" << "NormInf" << "\t" << solver_res << "\t\t" << "assemble_time" << "\t\t" << "pressure_load [kPa]" << "\n";
                            std::cout << " " << time << "\t" << delta << "\t\t" << iter << "\t" << norm_inf_rhs << "\t\t" << krylov_res << "\t\t" << assemble_time << "\t\t\t" << time*pressure_load/1000.0 << "\n";
                        }
                    }
                    
                    if(norm_inf_rhs<norm_inf_tol){
                        if (iter<=iter_min){
                            delta = success*delta;
                        }
                        if (iter>=iter_max){
                            delta = delta/failure;
                        }
                        FLAG1=1;
                        break;
                    }
                    
                    if(norm_inf_rhs>norm_inf_max||iter==iter_max){
                        nb_bis++;
                        if (nb_bis<nb_bis_max){
                            delta /= failure;
                            time -= delta;
                            *x = *epetra_y;
                            if (MyPID==0){
                                std::cout << "Bisecting increment: " << nb_bis << "\n";
                            }
                        }
                        else{
                            if (MyPID==0){
                                std::cout << "Bisection number exceeds.\n";
                            }
                            return 2;
                        }
                        FLAG2=1;
                        FLAG3=1;
                        break;
                    }
                    
                }
                
                epetra_rhsv = epetra_rhs;
                epetra_lhs->PutScalar(0.0);
                
                Comm->Barrier();
                if (time==time_init){
                    Thyra::initializeOp<double>(*lowsFactory, thyra_A, lows.ptr());
                }
                else{
                    Thyra::uninitializeOp<double>(*lowsFactory, lows.ptr());
                    Thyra::initializeAndReuseOp<double>(*lowsFactory, thyra_A, lows.ptr());
                }
                status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *thyra_rhs, thyra_lhs.ptr());
                /*if (status.solveStatus){
                    if (Comm->MyPID()==0){
                        std::cout << "STRATIMIKOS FAILED TO CONVERGE.\n";
                        std::cout << status.message << "\n";
                    }
                }*/
                //krylov_its = -1; //not sure that thyra_lhs will be copied into epetra_lhs
                krylov_res = status.achievedTol;
                Comm->Barrier();
                x->Update(1.0,*epetra_lhs,1.0);
                
                FLAG3=1;
            }
        }
    }
    
    return 0;
}

int Newton_Raphson::Solve_with_Aztec(){
    
    Epetra_Time Time(*Comm);
    
    Epetra_LinearProblem problem;
    AztecOO solver;
    
    Epetra_FECrsMatrix stiffness(Copy,*interface->FEGraph);
    Epetra_FEVector rhs(*interface->StandardMap);
    Epetra_Vector lhs(*interface->StandardMap);
    Epetra_Vector y(*interface->StandardMap);
     
    double Assemble_time;
    double displacement;
    double norm_inf_rhs;
    double time_max = 1.0;
    double Krylov_res = 0;
    int FLAG1, FLAG2, FLAG3, nb_bis, iter;
    int Krylov_its = 0;
    
    time = 0.0;
    
    std::string solver_its = "GMRES_its";
    std::string solver_res = "GMRES_res";
    
    FLAG1=1;
    while (FLAG1==1){
        
        FLAG1=0;
        FLAG2=1;
        FLAG3=1;
        nb_bis = 0;
        time += delta;
        y = *x;
        
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
            
            interface->pressure_load = time*pressure_load;
            
            iter = 0;
            while (FLAG3==1){
                
                FLAG3=0;
                iter++;
                
                if(iter>iter_max){
                    if (MyPID==0){std::cout << "Iteration limit exceeds.\n";}
                    return 1;
                }
                
                Time.ResetStartTime();
                interface->get_matrix_and_rhs(*x, stiffness, rhs);
                Assemble_time = Time.ElapsedTime();
                
                if (iter==1){
                    displacement = delta*bc_disp;
                }
                else{
                    displacement = 0.0;
                }
            
                interface->apply_dirichlet_conditions(stiffness, rhs, displacement);
                
                if(iter>1){
                    
                    rhs.NormInf(&norm_inf_rhs);
                    
                    if (MyPID==0){
                        if(iter>2){
                            /*std::cout << std::setprecision(2);
                            std::cout << std::scientific;*/
                            std::cout << "\t\t\t" << iter << "\t" << norm_inf_rhs << "\t" << Krylov_its << "\t\t" << Krylov_res << "\t\t" << Assemble_time << "\t\t\t" << time*pressure_load/1000.0 << "\n";
                        }
                        else{
                            std::cout << "\n Time" << "\t" << "Timestep" << "\t" << "Iter" << "\t" << "NormInf" << "\t\t" << solver_its << "\t" << solver_res << "\t\t" << "assemble_time" << "\t\t" << "pressure_load [kPa]" << "\n";
                            /*std::cout << std::setprecision(2);
                            std::cout << std::scientific;*/
                            std::cout << " " << time << "\t" << delta << "\t\t" << iter << "\t" << norm_inf_rhs << "\t" << Krylov_its << "\t\t" << Krylov_res << "\t\t" << Assemble_time << "\t\t\t" << time*pressure_load/1000.0 << "\n";
                        }
                    }
                    
                    if(norm_inf_rhs<norm_inf_tol){
                        if (iter<=iter_min){
                            delta = success*delta;
                        }
                        if (iter>=iter_max){
                            delta = delta/failure;
                        }
                        FLAG1=1;
                        //std::string name="Newton_Solution_Time" + std::to_string(time) + ".mtx";
                        //print_newton_solution(name);
                        break;
                    }
                    
                    if(norm_inf_rhs>norm_inf_max||iter==iter_max){
                        nb_bis++;
                        if (nb_bis<nb_bis_max){
                            delta /= failure;
                            time -= delta;
                            *x = y;
                            if (MyPID==0){
                                std::cout << "Bisecting increment: " << nb_bis << "\n";
                            }
                        }
                        else{
                            if (MyPID==0){
                                std::cout << "Bisection number exceeds.\n";
                            }
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
                
                solver.Iterate(2000,1e-6);
                Krylov_its = solver.NumIters();
                Krylov_res = solver.TrueResidual();
                x->Update(1.0,lhs,1.0);
                
                FLAG3=1;
                
            }
        }
    }
    
    //print_newton_solution("Newton_solution.mtx");
    return 0;
}


int Newton_Raphson::print_newton_solution(std::string fileName){
    
    int NumTargetElements = 0;
    if (MyPID==0){
        NumTargetElements = 3*interface->Mesh->n_nodes;
    }
    Epetra_Map MapOnRoot(-1,NumTargetElements,0,*Comm);
    Epetra_Export ExportOnRoot(*interface->StandardMap,MapOnRoot);
    Epetra_MultiVector lhs_root(MapOnRoot,true);
    lhs_root.Export(*x,ExportOnRoot,Insert);
    
    int error = EpetraExt::MultiVectorToMatrixMarketFile(fileName.c_str(),lhs_root,0,0,false);
    
    return error;
}

void Newton_Raphson::Initialization(){
    x->PutScalar(0.0);
}

void Newton_Raphson::setInitialization(Epetra_Vector & init){
    *x = init;
}

Newton_Raphson::~Newton_Raphson(){
}
