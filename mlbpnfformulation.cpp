#include "mlbpnfformulation.h"

#include "instance.h"
#include "solution.h"
#include "users.h"


void MLBPNFFormulation::createDecisionVariables(IloEnv env, const Instance<MLBP>& inst)
{

	// decision variables x_{kij}
	x = IloArray<IloArray<IloNumVarArray>>(env, inst.m);

	// flow variables f_{kij}
	f = IloArray<IloArray<IloNumVarArray>>(env, inst.m);

	// counters of how many decision variables have been added for debugging
	int count = 0;

	for (int k = 0; k < inst.m; k++) {
		x[k] = IloArray<IloNumVarArray>(env, inst.n[k]);
		f[k] = IloArray<IloNumVarArray>(env, inst.n[k]);
		for (int i : inst.B[k]) {
			x[k][i] = IloNumVarArray(env, inst.n[k + 1], 0, 1, ILOBOOL);
			// Idea: speedup if ILOBOOL for flow on level 0?
			f[k][i] = IloNumVarArray(env, inst.n[k + 1], 0, inst.n[0], ILOINT);
			count += inst.n[k + 1];
		}
	}
	MLB_OUT(TRACE) << "added " << count << " x_{kij} and f_{kij} variables." << std::endl;

	// decision variables y_{ki}
	y = IloArray<IloNumVarArray>(env, inst.m + 1);

	// counters of how many decision variables have been added for debugging
	int ys = 0;

	for (int k = 0; k <= inst.m; k++) {
		y[k] = IloNumVarArray(env, inst.n[k], 0, 1, ILOBOOL);
		ys += inst.n[k];
	}

	MLB_OUT(TRACE) << "added " << ys << " y_{ki} variables" << std::endl;
}

void MLBPNFFormulation::addConstraints(IloEnv env, IloModel model, const Instance<MLBP>& inst)
{
	// there can only be flow between 2 item/bins if the lower level bin is assigned to the higher level bin
	int count = 0;
	for (int k = 0; k < inst.m; k++) {
		for (int i : inst.B[k]) {
			for (int j : inst.B[k + 1]) {
				model.add(f[k][i][j] <= x[k][i][j] * inst.n[0]);
				count++;
			}
		}
	}
	MLB_OUT(TRACE) << "added " << count << " constraints such that there is only flow between bins if they are assigned to eachother." << std::endl;

	// there needs to be flow between two item/bins when they are assigned to eachother
	count = 0;
	for (int k = 0; k < inst.m; k++) {
		for (int i : inst.B[k]) {
			for (int j : inst.B[k + 1]) {
				model.add(x[k][i][j] <= f[k][i][j]);
				count++;
			}
		}
	}
	MLB_OUT(TRACE) << "added " << count << " constraints such that there is flow between two item/bins if they are assigned to eachtother." << std::endl;
	

	// the flow between each item and all the bins of level 1 needs to be 1
	for (int i : inst.B[0]) {
		IloExpr sum(env);
		for (int j : inst.B[1]) {
			sum += f[0][i][j];
		}
		model.add(sum == 1);
		sum.end();
	}
	MLB_OUT(TRACE) << "added " << inst.n[0] << " constraints such that each item provides 1 flow to the network." << std::endl;

	// the flow between the top level bins and the layer below is equal to the amount of items
	IloExpr sum(env);
	for (int i : inst.B[inst.m - 1]) {
		for (int j : inst.B[inst.m]) {
			sum += f[inst.m - 1][i][j];
		}
	}
	model.add(sum == inst.n[0]);
	sum.end();
	MLB_OUT(TRACE) << "added 1 constraint such that the flow to the top level bins equals the amount of items." << std::endl;

	// each bin apart from the ones on the top level have the same in- as outflow
	count = 0;
	for (int k = 1; k < inst.m; k++) {
		for (int i : inst.B[k]) {
			IloExpr sum(env);
			for (int lower : inst.B[k - 1]) {
				sum += f[k - 1][lower][i];
			}
			for (int upper : inst.B[k + 1]) {
				sum -= f[k][i][upper];
			}
			model.add(sum == 0);
			sum.end();
			count++;
		}
	}
	MLB_OUT(TRACE) << "added " << count << " constraints such each inner bins has the same in- as outflow." << std::endl;
	
	/************************************************************************/
	/** Constraints from basic Multi-Level Bin Packing Problem formulation **/
	/************************************************************************/

	// a bin can only be assigned to another bin if it is used
	count = 0;
	for (int k = 1; k < inst.m; k++) {
		for (int i : inst.B[k]) {
			for (int j : inst.B[k + 1]) {
				model.add(x[k][i][j] <= y[k][i]);
				count++;
			}
		}
	}
	MLB_OUT(TRACE) << "added " << count << " constraints such that only bins that are used are assigned to another bin" << std::endl;

	// each item must be inserted into exactly one bin of level 1
	for (int i : inst.B[0]) {
		IloExpr sum(env);
		for (int j : inst.B[1]) {
			sum += x[0][i][j];
		}
		model.add(sum == 1);
		sum.end();
	}
	MLB_OUT(TRACE) << "added " << inst.n[0] << " constraints such that each item is inserted in to exactly 1 bin of level 1" << std::endl;

	// each bin must be inserted into exactly one bin if it is used
	count = 0;
	for (int k = 1; k < inst.m; k++) {
		count += inst.n[k];
		for (int i : inst.B[k]) {
			IloExpr sum(env);
			for (int j : inst.B[k + 1]) {
				sum += x[k][i][j];
			}
			model.add(sum == y[k][i]);
			// model.add((sum == 1 && y[k][i] == 1) || (sum == 0 && y[k][i] == 0));
			sum.end();
		}
	}
	MLB_OUT(TRACE) << "added " << count << " constraints such that each bin must be inserted into exactly one bin if it is used" << std::endl;

	// the capacity of each used bin must not be exceeded
	count = 0;
	for (int k : inst.M) {
		for (int j : inst.B[k]) { // index of bin of which to check capacity
			IloExpr sum(env);
			for (int i : inst.B[k - 1]) { // index of the item/bin that was put into the bin of which to check capacity
				sum += x[k - 1][i][j] * inst.s[k - 1][i];
			}
			model.add(sum <= y[k][j] * inst.w[k][j]);
			count++;
			sum.end();
		}
	}


	MLB_OUT(TRACE) << "added " << count << " constraints such that the capacity of each used bin must not be exceeded" << std::endl;
}

void MLBPNFFormulation::addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBP>& inst)
{
	IloExpr sum(env);
	for (int k : inst.M) {
		for (int i : inst.B[k]) {
			sum += y[k][i] * inst.c[k][i];
		}
	}
	model.add(IloMinimize(env, sum));
	sum.end();
}

void MLBPNFFormulation::extractSolution(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol)
{
	sol.total_cost = 0;
	for (int k : inst.M) {
		for (int i : inst.B[k]) {
			if (cplex.getValue(y[k][i]) > 0.5)
				sol.total_cost += inst.c[k][i];
		}
	}

	for (int k = 0; k < inst.m; k++) {
		sol.item_to_bins[k].assign(inst.n[k], -1);

		for (int i : inst.B[k])
			for (int j : inst.B[k + 1])
				if (cplex.getValue(x[k][i][j]) > 0.5)
					sol.item_to_bins[k][i] = j;

	}
}

