#ifndef __MLBPTWNF_FORMULATION_H__
#define __MLBPTWNF_FORMULATION_H__


#include "problems.h"
#include "mipsolver.h"

template<typename> struct Instance;
template<typename> struct Solution;

class MLBPTWNFFormulation : public MIPFormulation<MLBPTW>
{
public:
	virtual void createDecisionVariables(IloEnv env, const Instance<MLBPTW>& inst);
	virtual void addConstraints(IloEnv env, IloModel model, const Instance<MLBPTW>& inst);
	virtual void addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBPTW>& inst);
	virtual void extractSolution(IloCplex cplex, const Instance<MLBPTW>& inst, Solution<MLBPTW>& sol);
private:

	// binary decision variables x_{kij}: item/bin of index i of level k is inserted into bin j of level k + 1 (=1) or not (=0)
	IloArray<IloArray<IloNumVarArray>> x;

	// binary decision variables y_{ki}: item/bin of index i of level k is used (=1) or not (=0)
	IloArray<IloNumVarArray> y;

	// integer decision variables f_{kij}: flow between item/bin of index i of level k to bin of index j of level k + 1 (1 flow means 1 item)
	IloArray<IloArray<IloNumVarArray>> f;

	// integer decision variables u_{i}: earliest packing time of item i
	IloArray<IloNumVar> u;

	// binary decision variables ib_{kij}: item i is assigned to bin j at level k
	IloArray<IloArray<IloNumVarArray>> ib;

};


#endif // __MLBPTWNF_FORMULATION_H__
