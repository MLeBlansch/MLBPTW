#include "mlbpformulation.h"

#include "instance.h"
#include "solution.h"
#include "users.h"


void MLBPFormulation::createDecisionVariables(IloEnv env, const Instance<MLBP>& inst)
{

	// decision variables x_{kij}
	x = IloArray<IloArray<IloNumVarArray>>(env, inst.m);	

	// counters of how many decision variables have been added for debugging
	int xs = 0;

	for (int k = 0; k < inst.m; k++) {
		x[k] = IloArray<IloNumVarArray>(env, inst.n[k]);
		for (int i : inst.B[k]) {
			x[k][i] = IloNumVarArray(env, inst.n[k+1], 0, 1, ILOBOOL);
			xs += inst.n[k + 1];
		}
	}
	MLB_OUT(TRACE) << "added " << xs << " x_{kij} variables" << std::endl;

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

void MLBPFormulation::addConstraints(IloEnv env, IloModel model, const Instance<MLBP>& inst)
{
	// a bin can only be assigned to another bin if it is used
	int count = 0;
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

void MLBPFormulation::addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBP>& inst)
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

void MLBPFormulation::extractSolution(IloCplex cplex, const Instance<MLBP>& inst, Solution<MLBP>& sol)
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

