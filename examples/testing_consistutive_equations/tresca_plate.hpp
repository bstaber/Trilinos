/*
Brian Staber (brian.staber@gmail.com)

1) von Mises plate under tension
2) Fenics and Zset give same results

*/

#ifndef TRESCA_PLATE_HPP
#define TRESCA_PLATE_HPP

#define _USE_MATH_DEFINES

#include "plasticitySmallStrains.hpp"

class tresca_plate : public plasticitySmallStrains
{
public:

  double E, nu, lambda, mu, kappa; // = 210000.0;
  double R0, H;
  double gasl, gasr, gala, gara, qs, qlsl, qlla, qrsr, qrra;
  bool IFPLAS, SMOOTH, LEFT, RIGHT;

  std::string BCType;

  Epetra_SerialDenseMatrix COMPLIANCE;

  tresca_plate(Epetra_Comm & comm, Teuchos::ParameterList & Parameters){

    if (Parameters.sublist("Behavior").isParameter("young")) {
      E  = Teuchos::getParameter<double>(Parameters.sublist("Behavior"), "young");
    }
    else{
      std::cout << "Young modulus not found, setting default value 210MPa" << std::endl;
      E = 210000.0;
    }

    if (Parameters.sublist("Behavior").isParameter("poisson")) {
      nu = Teuchos::getParameter<double>(Parameters.sublist("Behavior"), "poisson");
    }
    else{
      std::cout << "Poisson's coefficient not found, setting default value 0.3" << std::endl;
      nu = 0.3;
    }
    if (Parameters.sublist("Behavior").isParameter("yield")) {
      R0 = Teuchos::getParameter<double>(Parameters.sublist("Behavior"), "yield");
    }
    else{
      std::cout << "Yield stress not found, setting default value 250.0" << std::endl;
      R0 = 250.0;
    }
    if (Parameters.sublist("Behavior").isParameter("hardening")) {
      H = Teuchos::getParameter<double>(Parameters.sublist("Behavior"), "hardening");
    }
    else{
      std::cout << "Hardening modulus not found, setting default value 10000.0" << std::endl;
      H = 10000.0;
    }
    if (Parameters.sublist("BoundaryConditions").isParameter("type")) {
      BCType = Teuchos::getParameter<std::string>(Parameters.sublist("BoundaryConditions"), "type");
    }
    else {
      std::cout << "Boundary condition type not found, setting default value tension_totally_clamped" << std::endl;
      BCType = "tension_totally_clamepd";
    }

    plasticitySmallStrains::initialize(comm, Parameters);

    lambda = E*nu/((1.0+nu)*(1.0-2.0*nu));
    mu = E/(2.0*(1.0+nu));
    kappa = lambda + 2.0*mu/3.0;

    COMPLIANCE.Reshape(6,6);
  }

  ~tresca_plate(){
  }

  void constitutive_problem(const unsigned int & elid, const unsigned int & igp,
                            const Epetra_SerialDenseVector & DETO, Epetra_SerialDenseVector & SIG,
                            double & EPCUM, Epetra_SerialDenseMatrix & TGM){

    Epetra_SerialDenseVector SIGTR(6), DSIG(6), EYE(6);
    EYE(0) = 1.0; EYE(1) = 1.0; EYE(2) = 1.0; EYE(3) = 0.0; EYE(4) = 0.0; EYE(5) = 0.0;

    get_elasticity_tensor(elid, igp, ELASTICITY);
    get_compliance_tensor(elid, igp, COMPLIANCE);
    DSIG.Multiply('N','N',1.0,ELASTICITY,DETO,0.0);
    for (unsigned int k=0; k<6; ++k) SIGTR(k) = SIG(k) + DSIG(k);

    Epetra_SerialDenseVector EELTR(6), xtr(3), ytr(3), y(3);
    EELTR.Multiply('N','N',1.0,COMPLIANCE,SIGTR,0.0);
    sorted_eigenvalues(EELTR,xtr);
    ytr(0) = lambda*(xtr(0)+xtr(1)+xtr(2))+2.0*mu*xtr(0);
    ytr(1) = lambda*(xtr(0)+xtr(1)+xtr(2))+2.0*mu*xtr(1);
    ytr(2) = lambda*(xtr(0)+xtr(1)+xtr(2))+2.0*mu*xtr(2);

    //std::cout << "ytr = " << ytr << std::endl;
    //std::cout << "EELTR = " << EELTR << std::endl;
    //std::cout << "SIGTR = " << SIGTR << std::endl;
    //std::cout << "eigeeltr = " << xtr << std::endl;

    double yield = ytr(0)-ytr(2)-R0-H*EPCUM;

    /*std::cout << "yield = " << yield << std::endl;
    std::cout << "ytr(0) - ytr(2) = " << ytr(0)-ytr(2) << std::endl;
    std::cout << "radisu = " << R0+H*EPCUM << std::endl;
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "mu = " << mu << std::endl;
    std::cout << "ytr = " << ytr << std::endl;*/

    std::cout << "SIGTR = " << SIGTR << std::endl;
    std::cout << "EELTR = " << EELTR << std::endl;

    SIG = SIGTR;
    TGM = ELASTICITY;

    IFPLAS = false;
    SMOOTH = false;
    LEFT   = false;
    RIGHT  = false;

    double tolerance = 1.0e-12;

    if (yield>tolerance){
      IFPLAS = true;
      critical_values(ytr,EPCUM);

      if (qs<tolerance) {
        SMOOTH = true;
        double gamma_s = yield/(4.0*mu+H);
        y(0) = ytr(0) - 2.0*mu*gamma_s;
        y(1) = ytr(1);
        y(2) = ytr(2) + 2.0*mu*gamma_s;
        EPCUM += gamma_s;

        Epetra_SerialDenseVector E1(6), E2(6), E3(6);
        E1 = EigenProjection(EELTR,xtr(1),xtr(2),xtr(0));
        E2 = EigenProjection(EELTR,xtr(0),xtr(2),xtr(1));
        E3 = EigenProjection(EELTR,xtr(0),xtr(1),xtr(2));
        for (unsigned int k=0; k<6; ++k) SIG(k) = y(0)*E1(k) + y(1)*E2(k) + y(2)*E3(k);

        tgm_smooth(EELTR,xtr,y,E1,E2,E3,TGM);
      }

      if ((qlsl>=tolerance) & (qlla<tolerance)) {
        LEFT = true;
        if (SMOOTH==true) std::cout << "Can't be SMOOTH and LEFT at the same time" << std::endl;
        double gamma_l = (0.5*(ytr(0)+ytr(1))-ytr(2)-R0-H*EPCUM)/(3.0*mu+H);
        y(0) = 0.5*(ytr(0)+ytr(1)) - mu*gamma_l;
        y(1) = y(0);
        y(2) = ytr(2) + 2.0*mu*gamma_l;
        EPCUM += gamma_l;

        Epetra_SerialDenseVector E3(6), E12(6);
        E3 = EigenProjection(EELTR,xtr(0),xtr(1),xtr(2));
        for (unsigned int k=0; k<6; ++k) {
          E12(k) = EYE(k)-E3(k);
          SIG(k) = 0.5*(y(0)+y(1))*E12(k) + y(2)*E3(k);
        }

        tgm_left(EELTR,xtr,y,E12,E3,TGM);
      }

      if ((qrsr>=tolerance) & (qrra<tolerance)) {
        RIGHT = true;
        if ((SMOOTH==true) || (LEFT==true)) std::cout << "Can't be RIGHT and SMOOTH/LEFT at the same time: " << SMOOTH << "," << LEFT << "," << RIGHT << std::endl;
        double gamma_r = (ytr(0)-0.5*(ytr(1)+ytr(2))-R0-H*EPCUM)/(3.0*mu+H);
        y(0) = ytr(0) - 2.0*mu*gamma_r;
        y(1) = 0.5*(ytr(1)+ytr(2)) + mu*gamma_r;
        y(2) = y(1);
        std::cout << "Y_TR = " << ytr << std::endl;
        std::cout << "Y_NEW = " << y << std::endl;
        EPCUM += gamma_r;

        Epetra_SerialDenseVector E1(6), E23(6);
        E1 = EigenProjection(EELTR,xtr(1),xtr(2),xtr(0));
        for (unsigned int k=0; k<6; ++k) {
          E23(k) = EYE(k)-E1(k);
          SIG(k) = y(0)*E1(k) + 0.5*(y(1)+y(2))*E23(k);
        }

        /*std::cout << "E1 = " << E1 << std::endl;
        std::cout << "E23 = " << E23 << std::endl;*/

        tgm_right(EELTR,xtr,y,E1,E23,TGM);
      }
    }
    //std::cout << "States: " << LEFT << "\t" << SMOOTH << "\t" << RIGHT << std::endl;
  }

  void sorted_eigenvalues(Epetra_SerialDenseVector & X, Epetra_SerialDenseVector & eigx){
    double I1 = X(0)+X(1)+X(2);
    double I2 = X(0)*X(1)+X(0)*X(2)+X(1)*X(2)-0.5*X(3)*X(3)-0.5*X(4)*X(4)-0.5*X(5)*X(5);
    double I3 = X(0)*X(1)*X(2) - 0.5*X(1)*X(4)*X(4) - 0.5*X(2)*X(5)*X(5) - 0.5*X(0)*X(3)*X(3) + M_SQRT1_2*X(3)*X(4)*X(5);

    //std::cout << "INVARIANTS = " << I1 << "\t" << I2 << "\t" << I3 << std::endl;

    double Q = std::fmax(0.0,(1.0/9.0)*(I1*I1-3.0*I2));
    double R = (1.0/54.0)*(-2.0*I1*I1*I1+9.0*I1*I2-27.0*I3);
    double eta0 = 0.0;
    if (Q>0.0) eta0 = R/std::sqrt(Q*Q*Q);
    double eta = std::acos(std::fmin(std::fmax(eta0,-1.0),1.0))/3.0;

    double sqrtQ = std::sqrt(Q);
    eigx(0) = -2.0*sqrtQ*std::cos(eta+2.0*M_PI/3.0) + I1/3.0;
    eigx(1) = -2.0*sqrtQ*std::cos(eta-2.0*M_PI/3.0) + I1/3.0;
    eigx(2) = -2.0*sqrtQ*std::cos(eta) + I1/3.0;
  }

  void critical_values(Epetra_SerialDenseVector & v, double & EPCUM){
    gasl = (v(0)-v(1))/(2.0*mu);
    gasr = (v(1)-v(2))/(2.0*mu);
    gala = (v(0)+v(1)-2.0*v(2))/(6.0*mu);
    gara = (2.0*v(0)-v(1)-v(2))/(6.0*mu);

    double minga = std::fmin(gasl,gasr);
    double rad = R0+H*(EPCUM+minga);
    qs = v(0)-v(2)-rad-4.0*mu*minga;

    rad = R0+H*(EPCUM+gasl);
    qlsl = 0.5*(v(0)+v(1))-v(2)-rad-3.0*mu*gasl;
    rad = R0+H*(EPCUM+gala);
    qlla = 0.5*(v(0)+v(1))-v(2)-rad-3.0*mu*gala;

    rad = R0+H*(EPCUM+gasr);
    qrsr = v(0)-0.5*(v(1)+v(2))-rad-3.0*mu*gasr;
    rad = R0+H*(EPCUM+gara);
    qrra = v(0)-0.5*(v(1)+v(2))-rad-3.0*mu*gara;
  }

  void tgm_smooth(const Epetra_SerialDenseVector & X, const Epetra_SerialDenseVector & x, const Epetra_SerialDenseVector & y,
                  const Epetra_SerialDenseVector & E1, const Epetra_SerialDenseVector & E2, const Epetra_SerialDenseVector & E3,
                  Epetra_SerialDenseMatrix & TGM) {
    Epetra_SerialDenseVector EYE(6);
    EYE(0) = 1.0; EYE(1) = 1.0; EYE(2) = 1.0; EYE(3) = 0.0; EYE(4) = 0.0; EYE(5) = 0.0;
    Epetra_SerialDenseMatrix DXX(6,6), ikjl_EYEX(6,6), iljk_EYEX(6,6), ikjl_XEYE(6,6), iljk_XEYE(6,6);

    ikjl_EYEX = ikjl(EYE,X);
    ikjl_XEYE = ikjl(X,EYE);
    iljk_EYEX = iljk(EYE,X);
    iljk_XEYE = iljk(X,EYE);

    DXX = ikjl_EYEX;
    DXX += iljk_EYEX;
    DXX += ikjl_XEYE;
    DXX += iljk_XEYE;
    DXX.Scale(0.5);

    Epetra_SerialDenseMatrix E1E1(6,6), E2E2(6,6), E3E3(6,6);
    E1E1 = ijkl(E1,E1);
    E2E2 = ijkl(E2,E2);
    E3E3 = ijkl(E3,E3);

    Epetra_SerialDenseMatrix UNIT(6,6), ikjl_EYE(6,6), iljk_EYE(6,6);
    ikjl_EYE = ikjl(EYE,EYE);
    iljk_EYE = iljk(EYE,EYE);
    UNIT = ikjl_EYE;
    UNIT += iljk_EYE;
    UNIT.Scale(0.5);

    Epetra_SerialDenseMatrix DE1(6,6), DE2(6,6), DE3(6,6);
    for (unsigned int i=0; i<6; ++i) {
      for (unsigned int j=0; j<6; ++j) {
        DE1(i,j) = (DXX(i,j)-(x(1)+x(2))*UNIT(i,j)-(2.0*x(0)-x(1)-x(2))*E1E1(i,j)-(x(1)-x(2))*(E2E2(i,j)-E3E3(i,j)))/((x(0)-x(1))*(x(0)-x(2)));
        DE2(i,j) = (DXX(i,j)-(x(0)+x(2))*UNIT(i,j)-(2.0*x(1)-x(0)-x(2))*E2E2(i,j)-(x(0)-x(2))*(E1E1(i,j)-E3E3(i,j)))/((x(1)-x(0))*(x(1)-x(2)));
        DE3(i,j) = (DXX(i,j)-(x(0)+x(1))*UNIT(i,j)-(2.0*x(2)-x(0)-x(1))*E3E3(i,j)-(x(0)-x(1))*(E1E1(i,j)-E2E2(i,j)))/((x(2)-x(0))*(x(2)-x(1)));
      }
    }

    Epetra_SerialDenseVector TEMP(6), DLAMBDA(6);
    for (unsigned int k=0; k<6; ++k) TEMP(k) = 2.0*mu*E1(k)-2.0*mu*E3(k);
    DLAMBDA = TEMP;
    DLAMBDA.Scale(1.0/(4.0*mu+H));

    Epetra_SerialDenseMatrix ijkl_EYE(6,6), ijkl_TEMPDLAMBDA(6,6);
    ijkl_EYE = ijkl(EYE,EYE);
    ijkl_TEMPDLAMBDA = ijkl(TEMP,DLAMBDA);

    for (unsigned int i=0; i<6; ++i) {
      for (unsigned int j=0; j<6; ++j) {
        TGM(i,j) = y(0)*DE1(i,j) + y(1)*DE2(i,j) + y(2)*DE3(i,j) + 2.0*mu*E1E1(i,j) + 2.0*mu*E2E2(i,j) + 2.0*mu*E3E3(i,j) + lambda*ijkl_EYE(i,j) - ijkl_TEMPDLAMBDA(i,j);
      }
    }
  }

  void tgm_left(const Epetra_SerialDenseVector & X, const Epetra_SerialDenseVector & x, const Epetra_SerialDenseVector & y,
                const Epetra_SerialDenseVector & E12, const Epetra_SerialDenseVector & E3,
                Epetra_SerialDenseMatrix & TGM) {
    Epetra_SerialDenseVector EYE(6);
    EYE(0) = 1.0; EYE(1) = 1.0; EYE(2) = 1.0; EYE(3) = 0.0; EYE(4) = 0.0; EYE(5) = 0.0;
    Epetra_SerialDenseMatrix DXX(6,6), ikjl_EYEX(6,6), iljk_EYEX(6,6), ikjl_XEYE(6,6), iljk_XEYE(6,6);

    ikjl_EYEX = ikjl(EYE,X);
    ikjl_XEYE = ikjl(X,EYE);
    iljk_EYEX = iljk(EYE,X);
    iljk_XEYE = iljk(X,EYE);

    DXX = ikjl_EYEX;
    DXX += iljk_EYEX;
    DXX += ikjl_XEYE;
    DXX += iljk_XEYE;
    DXX.Scale(0.5);

    Epetra_SerialDenseMatrix E12E12(6,6), E3E3(6,6), E3E12(6,6), E12E3(6,6);
    E12E12 = ijkl(E12,E12);
    E3E3   = ijkl(E3,E3);
    E12E3  = ijkl(E12,E3);
    E3E12  = ijkl(E3,E12);

    Epetra_SerialDenseMatrix UNIT(6,6), ikjl_EYE(6,6), iljk_EYE(6,6);
    ikjl_EYE = ikjl(EYE,EYE);
    iljk_EYE = iljk(EYE,EYE);
    UNIT = ikjl_EYE;
    UNIT += iljk_EYE;
    UNIT.Scale(0.5);

    Epetra_SerialDenseMatrix DE3(6,6), ijkl_XE12(6,6), ijkl_E12X(6,6);
    ijkl_XE12 = ijkl(X,E12);
    ijkl_E12X = ijkl(E12,X);
    for (unsigned int i=0; i<6; ++i) {
      for (unsigned int j=0; j<6; ++j) {
        DE3(i,j) = (DXX(i,j)-(x(0)+x(1))*UNIT(i,j)-ijkl_XE12(i,j)-ijkl_E12X(i,j) + (x(0)+x(1))*E12E12(i,j) + (x(0)+x(1)-2.0*x(2))*E3E3(i,j) + x(2)*(E3E12(i,j)+E12E3(i,j)))/((x(2)-x(0))*(x(2)-x(1)));
      }
    }

    Epetra_SerialDenseVector TEMP(6), DLAMBDA(6);
    for (unsigned int k=0; k<6; ++k) TEMP(k) = mu*E12(k)-2.0*mu*E3(k);
    DLAMBDA = TEMP;
    DLAMBDA.Scale(1.0/(3.0*mu+H));

    Epetra_SerialDenseMatrix ijkl_EYE(6,6), ijkl_TEMPDLAMBDA(6,6);
    ijkl_EYE = ijkl(EYE,EYE);
    ijkl_TEMPDLAMBDA = ijkl(TEMP,DLAMBDA);

    for (unsigned int i=0; i<6; ++i) {
      for (unsigned int j=0; j<6; ++j) {
        TGM(i,j) = (y(2)-0.5*(y(0)+y(1)))*DE3(i,j) + mu*E12E12(i,j) + 2.0*mu*E3E3(i,j) + lambda*ijkl_EYE(i,j) - ijkl_TEMPDLAMBDA(i,j);
      }
    }
  }

  void tgm_right(const Epetra_SerialDenseVector & X, const Epetra_SerialDenseVector & x, const Epetra_SerialDenseVector & y,
                 const Epetra_SerialDenseVector & E1, const Epetra_SerialDenseVector & E23,
                 Epetra_SerialDenseMatrix & TGM) {
    Epetra_SerialDenseVector EYE(6);
    EYE(0) = 1.0; EYE(1) = 1.0; EYE(2) = 1.0; EYE(3) = 0.0; EYE(4) = 0.0; EYE(5) = 0.0;
    Epetra_SerialDenseMatrix DXX(6,6), ikjl_EYEX(6,6), iljk_EYEX(6,6), ikjl_XEYE(6,6), iljk_XEYE(6,6);

    ikjl_EYEX = ikjl(EYE,X);
    ikjl_XEYE = ikjl(X,EYE);
    iljk_EYEX = iljk(EYE,X);
    iljk_XEYE = iljk(X,EYE);

    DXX = ikjl_EYEX;
    DXX += iljk_EYEX;
    DXX += ikjl_XEYE;
    DXX += iljk_XEYE;
    DXX.Scale(0.5);

    Epetra_SerialDenseMatrix E23E23(6,6), E1E1(6,6), E1E23(6,6), E23E1(6,6);
    E23E23 = ijkl(E23,E23);
    E1E1   = ijkl(E1,E1);
    E23E1  = ijkl(E23,E1);
    E1E23  = ijkl(E1,E23);

    Epetra_SerialDenseMatrix UNIT(6,6), ikjl_EYE(6,6), iljk_EYE(6,6);
    ikjl_EYE = ikjl(EYE,EYE);
    iljk_EYE = iljk(EYE,EYE);
    UNIT = ikjl_EYE;
    UNIT += iljk_EYE;
    UNIT.Scale(0.5);

    Epetra_SerialDenseMatrix DE1(6,6), ijkl_XE23(6,6), ijkl_E23X(6,6);
    ijkl_XE23 = ijkl(X,E23);
    ijkl_E23X = ijkl(E23,X);
    for (unsigned int i=0; i<6; ++i) {
      for (unsigned int j=0; j<6; ++j) {
        DE1(i,j) = (DXX(i,j)-(x(1)+x(2))*UNIT(i,j)-ijkl_XE23(i,j)-ijkl_E23X(i,j) + (x(1)+x(2))*E23E23(i,j) + (x(1)+x(2)-2.0*x(0))*E1E1(i,j) + x(0)*(E1E23(i,j)+E23E1(i,j)))/((x(0)-x(1))*(x(0)-x(2)));
      }
    }

    Epetra_SerialDenseVector TEMP(6), DLAMBDA(6);
    for (unsigned int k=0; k<6; ++k) TEMP(k) = 2.0*mu*E1(k)-mu*E23(k);
    DLAMBDA = TEMP;
    DLAMBDA.Scale(1.0/(3.0*mu+H));

    Epetra_SerialDenseMatrix ijkl_EYE(6,6), ijkl_TEMPDLAMBDA(6,6);
    ijkl_EYE = ijkl(EYE,EYE);
    ijkl_TEMPDLAMBDA = ijkl(TEMP,DLAMBDA);

    for (unsigned int i=0; i<6; ++i) {
      for (unsigned int j=0; j<6; ++j) {
        TGM(i,j) = (y(0)-0.5*(y(1)+y(2)))*DE1(i,j) + mu*E23E23(i,j) + 2.0*mu*E1E1(i,j) + lambda*ijkl_EYE(i,j) - ijkl_TEMPDLAMBDA(i,j);
      }
    }
  }

  void get_elasticity_tensor(const unsigned int & e_lid, const unsigned int & gp, Epetra_SerialDenseMatrix & tgm){
    double c1 = lambda+2.0*mu;
    tgm(0,0) = c1;     tgm(0,1) = lambda; tgm(0,2) = lambda; tgm(0,3) = 0.0;     tgm(0,4) = 0.0;    tgm(0,5) = 0.0;
    tgm(1,0) = lambda; tgm(1,1) = c1;     tgm(1,2) = lambda; tgm(1,3) = 0.0;     tgm(1,4) = 0.0;    tgm(1,5) = 0.0;
    tgm(2,0) = lambda; tgm(2,1) = lambda; tgm(2,2) = c1;     tgm(2,3) = 0.0;     tgm(2,4) = 0.0;    tgm(2,5) = 0.0;
    tgm(3,0) = 0.0;    tgm(3,1) = 0.0;    tgm(3,2) = 0.0;    tgm(3,3) = 2.0*mu;  tgm(3,4) = 0.0;    tgm(3,5) = 0.0;
    tgm(4,0) = 0.0;    tgm(4,1) = 0.0;    tgm(4,2) = 0.0;    tgm(4,3) = 0.0;     tgm(4,4) = 2.0*mu; tgm(4,5) = 0.0;
    tgm(5,0) = 0.0;    tgm(5,1) = 0.0;    tgm(5,2) = 0.0;    tgm(5,3) = 0.0;     tgm(5,4) = 0.0;    tgm(5,5) = 2.0*mu;
  }

  void get_compliance_tensor(const unsigned int & e_id, const unsigned int & gp, Epetra_SerialDenseMatrix & itgm){
    double c1 = 1.0/E;
    double c2 = -nu*c1;
    double c3 = (1.0+nu)*c1;
    itgm(0,0) = c1;     itgm(0,1) = c2;     itgm(0,2) = c2;     itgm(0,3) = 0.0;     itgm(0,4) = 0.0;    itgm(0,5) = 0.0;
    itgm(1,0) = c2;     itgm(1,1) = c1;     itgm(1,2) = c2;     itgm(1,3) = 0.0;     itgm(1,4) = 0.0;    itgm(1,5) = 0.0;
    itgm(2,0) = c2;     itgm(2,1) = c2;     itgm(2,2) = c1;     itgm(2,3) = 0.0;     itgm(2,4) = 0.0;    itgm(2,5) = 0.0;
    itgm(3,0) = 0.0;    itgm(3,1) = 0.0;    itgm(3,2) = 0.0;    itgm(3,3) = c3;      itgm(3,4) = 0.0;    itgm(3,5) = 0.0;
    itgm(4,0) = 0.0;    itgm(4,1) = 0.0;    itgm(4,2) = 0.0;    itgm(4,3) = 0.0;     itgm(4,4) = c3;     itgm(4,5) = 0.0;
    itgm(5,0) = 0.0;    itgm(5,1) = 0.0;    itgm(5,2) = 0.0;    itgm(5,3) = 0.0;     itgm(5,4) = 0.0;    itgm(5,5) = c3;
  }

  void setup_dirichlet_conditions(){
    if (BCType=="tension_totally_clamped") setup_tension_totally_clamped();
    if (BCType=="tension_with_slip") setup_tension_with_slip();
  }

  void setup_tension_totally_clamped(){
      n_bc_dof = 0;
      double x,y,z;
      unsigned int node;
      for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
          node = Mesh->local_nodes[i];
          x = Mesh->nodes_coord[3*node+0];
          y = Mesh->nodes_coord[3*node+1];
          z = Mesh->nodes_coord[3*node+2];
          if(y==0.0){
              n_bc_dof+=3;
          }
          if(y==10.0){
              n_bc_dof+=1;
          }
      }

      int indbc = 0;
      dof_on_boundary = new int [n_bc_dof];
      for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
          node = Mesh->local_nodes[inode];
          x = Mesh->nodes_coord[3*node+0];
          y = Mesh->nodes_coord[3*node+1];
          z = Mesh->nodes_coord[3*node+2];
          if(y==0.0){
              dof_on_boundary[indbc+0] = 3*inode+0;
              dof_on_boundary[indbc+1] = 3*inode+1;
              dof_on_boundary[indbc+2] = 3*inode+2;
              indbc+=3;
          }
          if(y==10.0){
              dof_on_boundary[indbc] = 3*inode+1;
              indbc+=1;
          }
      }
  }

  void setup_tension_with_slip(){
      n_bc_dof = 0;
      double x,y,z;
      unsigned int node;
      for (unsigned int i=0; i<Mesh->n_local_nodes_without_ghosts; ++i){
          node = Mesh->local_nodes[i];
          x = Mesh->nodes_coord[3*node+0];
          y = Mesh->nodes_coord[3*node+1];
          z = Mesh->nodes_coord[3*node+2];
          if (y==0.0||y==10.0){
            if (x==0.0&&y==0.0&&z==0.0) {
              n_bc_dof+=1;
            }
            n_bc_dof+=1;
            if ((x==0.0&&y==0.0&&z==0.0)||(x==4.0&&y==0.0&&z==0.0)) {
              n_bc_dof+=1;
            }
          }
      }

      int indbc = 0;
      dof_on_boundary = new int [n_bc_dof];
      for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
          node = Mesh->local_nodes[inode];
          x = Mesh->nodes_coord[3*node+0];
          y = Mesh->nodes_coord[3*node+1];
          z = Mesh->nodes_coord[3*node+2];
          if (y==0.0||y==10.0){
            if (x==0.0&&y==0.0&&z==0.0) {
              dof_on_boundary[indbc] = 3*inode+0;
              indbc+=1;
            }
            dof_on_boundary[indbc] = 3*inode+1;
            indbc+=1;
            if ((x==0.0&&y==0.0&&z==0.0)||(x==4.0&&y==0.0&&z==0.0)) {
              dof_on_boundary[indbc] = 3*inode+2;
              indbc+=1;
            }
          }
      }

      /*for (unsigned int i=0; i<n_bc_dof; ++i) {
        std::cout << dof_on_boundary[i] << std::endl;
      }*/
  }

  void apply_dirichlet_conditions(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double factor){
    if (BCType=="tension_totally_clamped") apply_tension_totally_clamped(K,F,factor);
    if (BCType=="tension_with_slip") apply_tension_with_slip(K,F,factor);
  }

  void apply_tension_totally_clamped(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double factor){
      Epetra_MultiVector v(*StandardMap,true);
      v.PutScalar(0.0);

      int node;
      double x,y,z;
      for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
          node = Mesh->local_nodes[inode];
          x = Mesh->nodes_coord[3*node+0];
          y = Mesh->nodes_coord[3*node+1];
          z = Mesh->nodes_coord[3*node+2];
          if (y==0.0){
              v[0][StandardMap->LID(3*node+0)] = 0.0;
              v[0][StandardMap->LID(3*node+1)] = 0.0;
              v[0][StandardMap->LID(3*node+2)] = 0.0;
          }
          if (y==10.0){
              v[0][StandardMap->LID(3*node+1)] = factor*y;
          }
      }

      Epetra_MultiVector rhs_dir(*StandardMap,true);
      K.Apply(v,rhs_dir);
      F.Update(-1.0,rhs_dir,1.0);

      for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
          node = Mesh->local_nodes[inode];
          x = Mesh->nodes_coord[3*node+0];
          y = Mesh->nodes_coord[3*node+1];
          z = Mesh->nodes_coord[3*node+2];
          if (y==0.0){
              F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
              F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
              F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
          }
          if (y==10.0){
              F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
          }
      }
      ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
  }

  void apply_tension_with_slip(Epetra_FECrsMatrix & K, Epetra_FEVector & F, double factor){
    Epetra_MultiVector v(*StandardMap,true);
    v.PutScalar(0.0);

    int node;
    double x,y,z;
    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        x = Mesh->nodes_coord[3*node+0];
        y = Mesh->nodes_coord[3*node+1];
        z = Mesh->nodes_coord[3*node+2];
        if (y==0.0||y==10.0){
          if (x==0.0&&y==0.0&&z==0.0) {
            v[0][StandardMap->LID(3*node+0)] = 0.0;
          }
          v[0][StandardMap->LID(3*node+1)] = factor*y;
          if ((x==0.0&&y==0.0&&z==0.0)||(x==4.0&&y==0.0&&z==0.0)) {
            v[0][StandardMap->LID(3*node+2)] = 0.0;
          }
        }
    }

    Epetra_MultiVector rhs_dir(*StandardMap,true);
    K.Apply(v,rhs_dir);
    F.Update(-1.0,rhs_dir,1.0);

    for (unsigned int inode=0; inode<Mesh->n_local_nodes_without_ghosts; ++inode){
        node = Mesh->local_nodes[inode];
        x = Mesh->nodes_coord[3*node+0];
        y = Mesh->nodes_coord[3*node+1];
        z = Mesh->nodes_coord[3*node+2];
        if (y==0.0||y==10.0){
          if (x==0.0&&y==0.0&&z==0.0) {
            F[0][StandardMap->LID(3*node+0)] = v[0][StandardMap->LID(3*node+0)];
          }
          F[0][StandardMap->LID(3*node+1)] = v[0][StandardMap->LID(3*node+1)];
          if ((x==0.0&&y==0.0&&z==0.0)||(x==4.0&&y==0.0&&z==0.0)){
            F[0][StandardMap->LID(3*node+2)] = v[0][StandardMap->LID(3*node+2)];
          }
        }
    }
    ML_Epetra::Apply_OAZToMatrix(dof_on_boundary,n_bc_dof,K);
  }

  Epetra_SerialDenseVector prod(const Epetra_SerialDenseVector & A, const Epetra_SerialDenseVector & B) {
    Epetra_SerialDenseVector ret(6);
    ret(0) = A(0)*B(0) + (A(4)*B(4))/2.0 + (A(5)*B(5))/2.0;
    ret(1) = A(1)*B(1) + (A(3)*B(3))/2.0 + (A(5)*B(5))/2.0;
    ret(2) = A(2)*B(2) + (A(3)*B(3))/2.0 + (A(4)*B(4))/2.0;
    ret(3) = M_SQRT2*(A(5)*B(4)/2.0 + M_SQRT2*A(1)*B(3)/2.0 + M_SQRT2*A(3)*B(2)/2.0);
    ret(4) = M_SQRT2*(A(5)*B(3)/2.0 + M_SQRT2*A(0)*B(4)/2.0 + M_SQRT2*A(4)*B(2)/2.0);
    ret(5) = M_SQRT2*(A(4)*B(3)/2.0 + M_SQRT2*A(0)*B(5)/2.0 + M_SQRT2*A(5)*B(1)/2.0);
    return ret;
  }

  Epetra_SerialDenseVector EigenProjection(const Epetra_SerialDenseVector & A, double & a, double & b, double & c) {
    Epetra_SerialDenseVector E(6);
    E(0) = A(4)*A(4)/2.0 + A(5)*A(5)/2.0 + (A(0)-a)*(A(0)-b);
    E(1) = A(3)*A(3)/2.0 + A(5)*A(5)/2.0 + (A(1)-a)*(A(1)-b);
    E(2) = A(3)*A(3)/2.0 + A(4)*A(4)/2.0 + (A(2)-a)*(A(2)-b);
    E(3) = M_SQRT2*(A(4)*A(5)/2.0 + M_SQRT2*A(3)*(A(1)-a)/2.0 + M_SQRT2*A(3)*(A(2)-b)/2.0);
    E(4) = M_SQRT2*(A(3)*A(5)/2.0 + M_SQRT2*A(4)*(A(0)-a)/2.0 + M_SQRT2*A(4)*(A(2)-b)/2.0);
    E(5) = M_SQRT2*(A(3)*A(4)/2.0 + M_SQRT2*A(5)*(A(0)-a)/2.0 + M_SQRT2*A(5)*(A(1)-b)/2.0);

    double d = (c-a)*(c-b);
    double f = 1.0/d;
    E.Scale(f);
    return E;
  }

  Epetra_SerialDenseMatrix ikjl(const Epetra_SerialDenseVector & A, const Epetra_SerialDenseVector & B) {
    Epetra_SerialDenseMatrix ret(6,6);
    ret(0,0)=A(0)*B(0);               ret(0,1)=(A(5)*B(5))/2.0;         ret(0,2)=(A(4)*B(4))/2.0;          ret(0,3)=(M_SQRT2*A(5)*B(4))/2.0;  ret(0,4)=A(0)*B(4);               ret(0,5)=A(0)*B(5);
    ret(1,0)=(A(5)*B(5))/2.0;         ret(1,1)=A(1)*B(1);               ret(1,2)=(A(3)*B(3))/2.0;          ret(1,3)=A(1)*B(3);                ret(1,4)=(M_SQRT2*A(5)*B(3))/2.0; ret(1,5)=A(5)*B(1);
    ret(2,0)=(A(4)*B(4))/2.0;         ret(2,1)=(A(3)*B(3))/2.0;         ret(2,2)=A(2)*B(2);                ret(2,3)=A(3)*B(2);                ret(2,4)=A(4)*B(2);               ret(2,5)=(M_SQRT2*A(4)*B(3))/2.0;
    ret(3,0)=(M_SQRT2*A(5)*B(4))/2.0; ret(3,1)=A(1)*B(3);               ret(3,2)=A(3)*B(2);                ret(3,3)=2.0*A(1)*B(2);              ret(3,4)=M_SQRT2*A(5)*B(2);       ret(3,5)=A(5)*B(3);
    ret(4,0)=A(0)*B(4);               ret(4,1)=(M_SQRT2*A(5)*B(3))/2.0; ret(4,2)=A(4)*B(2);                ret(4,3)=M_SQRT2*A(5)*B(2);        ret(4,4)=2.0*A(0)*B(2);           ret(4,5)=M_SQRT2*A(0)*B(3);
    ret(5,0)=A(0)*B(5);               ret(5,1)=A(5)*B(1);               ret(5,2)=(M_SQRT2*A(4)*B(3))/2.0;  ret(5,3)=A(5)*B(3);                ret(5,4)=M_SQRT2*A(0)*B(3);       ret(5,5)=2.0*A(0)*B(1);
    return ret;
  }

  Epetra_SerialDenseMatrix iljk(const Epetra_SerialDenseVector & A, const Epetra_SerialDenseVector & B) {
    Epetra_SerialDenseMatrix ret(6,6);
    ret(0,0)=A(0)*B(0);               ret(0,1)=(A(5)*B(5))/2.0;         ret(0,2)=(A(4)*B(4))/2.0;         ret(0,3)=(M_SQRT2*A(4)*B(5))/2.0; ret(0,4)=A(4)*B(0);               ret(0,5)=A(5)*B(0);
    ret(1,0)=(A(5)*B(5))/2.0;         ret(1,1)=A(1)*B(1);               ret(1,2)=(A(3)*B(3))/2.0;         ret(1,3)=A(3)*B(1);               ret(1,4)=(M_SQRT2*A(3)*B(5))/2.0; ret(1,5)=A(1)*B(5);
    ret(2,0)=(A(4)*B(4))/2.0;         ret(2,1)=(A(3)*B(3))/2.0;         ret(2,2)=A(2)*B(2);               ret(2,3)=A(2)*B(3);               ret(2,4)=A(2)*B(4);               ret(2,5)=(M_SQRT2*A(3)*B(4))/2.0;
    ret(3,0)=(M_SQRT2*A(5)*B(4))/2.0; ret(3,1)=A(1)*B(3);               ret(3,2)=A(3)*B(2);               ret(3,3)=A(3)*B(3);               ret(3,4)=A(3)*B(4);               ret(3,5)=M_SQRT2*A(1)*B(4);
    ret(4,0)=A(0)*B(4);               ret(4,1)=(M_SQRT2*A(5)*B(3))/2.0; ret(4,2)=A(4)*B(2);               ret(4,3)=A(4)*B(3);               ret(4,4)=A(4)*B(4);               ret(4,5)=A(5)*B(4);
    ret(5,0)=A(0)*B(5);               ret(5,1)=A(5)*B(1);               ret(5,2)=(M_SQRT2*A(4)*B(3))/2.0; ret(5,3)=M_SQRT2*A(4)*B(1);       ret(5,4)=A(4)*B(5);               ret(5,5)=A(5)*B(5);
    return ret;
  }

  Epetra_SerialDenseMatrix ijkl(const Epetra_SerialDenseVector & A, const Epetra_SerialDenseVector & B) {
    Epetra_SerialDenseMatrix ret(6,6);
    for (unsigned int i=0; i<6; ++i) {
      for (unsigned int j=0; j<6; ++j) {
        ret(i,j) = A(i)*B(j);
      }
    }
    return ret;
  }

};
#endif
