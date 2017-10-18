#include  <math.h>
#include "Compressible_Mooney_Transverse_Isotropic.hpp"
#include "Newton_Raphsonpp.hpp"
#include "distributenrldata.hpp"

class objectiveFunction
{
private:
    
    Teuchos::ParameterList _paramList;
    Epetra_Comm * comm;
    Teuchos::RCP<Newton_Raphson> newton;
    Teuchos::RCP<TIMooney> interface;
    
    
public:
        
    objectiveFunction(Epetra_Comm & Comm, Teuchos::ParameterList & paramList){
        
        comm = &Comm;
        _paramList = paramList;
        interface = Teuchos::rcp(new TIMooney(Comm,paramList));
        newton = Teuchos::rcp(new Newton_Raphson(*interface,paramList));
        
        Teuchos::RCP<distributenrldata> nrldata = Teuchos::rcp(new distributenrldata(*interface->Mesh));
    }
    
    ~objectiveFunction(){
    }
    
};
