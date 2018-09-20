#ifndef DAMAGEFIELD_HPP
#define DAMAGEFIELD_HPP

include "linearFiniteElementProblem.hpp"

class damageField : public linearFiniteElementProblem
{
public:
    mesh * Mesh;
    
    damageField(mesh & Mesh);
    ~damageField();

    void create_FECrsGraph();
};

#endif
