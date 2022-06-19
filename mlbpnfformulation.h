#ifndef __MLBPNF_FORMULATION_H__
#define __MLBPNF_FORMULATION_H__


#include "problems.h"
#include "mipsolver.h"

template<typename> struct Instance;
template<typename> struct Solution;

class MLBPNFFormulation : public MIPFormulation<MLBP>
{
public:
	virtual void createDecisionVariables(IloEnv env, const Instance<MLBP>& inst);
	virtual void addConstraints(IloEnv env, IloModel model, const Instance<MLBP>& inst);
	virtual void addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBP>& inst);
	virtual void extractSolution(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol);
private:

	// binary decision variables x_{kij}: item/bin of index i of level k is inserted into bin j of level k + 1 (=1) or not (=0)
	IloArray<IloArray<IloNumVarArray>> x;

	// binary decision variables y_{ki}: item/bin of index i of level k is used (=1) or not (=0)
	IloArray<IloNumVarArray> y;

	// integer decision variables f_{kij}: flow between item/bin of index i of level k to bin of index j of level k + 1 (1 flow means 1 item)
	IloArray<IloArray<IloNumVarArray>> f;
};


#endif // __MLBPNF_FORMULATION_H__
