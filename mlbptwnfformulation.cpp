#include "mlbptwnfformulation.h"

#include "instance.h"
#include "solution.h"
#include "users.h"

void MLBPTWNFFormulation::createDecisionVariables(IloEnv env, const Instance<MLBPTW>& inst)
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
			f[k][i] = IloNumVarArray(env, inst.n[k + 1], 0, inst.n[0], ILOINT);
			count += inst.n[k + 1];
		}
	}
	MLB_OUT(TRACE) << "added " << count << " x_{kij} and f_{kij} variables." << std::endl;

	// decision variables ib_{kij}
	ib = IloArray<IloArray<IloNumVarArray>>(env, inst.m + 1);

	count = 0;
	for (int k = 0; k <= inst.m; k++) {
		ib[k] = IloArray<IloNumVarArray>(env, inst.n[0]);
		for (int i : inst.B[0]) {
			ib[k][i] = IloNumVarArray(env, inst.n[k], 0, 1, ILOBOOL);
			count += inst.n[k];
		}
	}

	MLB_OUT(TRACE) << "added " << count << " ib_{kij} variables" << std::endl;



	// decision variables y_{ki}
	y = IloArray<IloNumVarArray>(env, inst.m + 1);

	// counters of how many decision variables have been added for debugging
	count = 0;
	for (int k = 0; k <= inst.m; k++) {
		y[k] = IloNumVarArray(env, inst.n[k], 0, 1, ILOBOOL);
		count += inst.n[k];
	}

	MLB_OUT(TRACE) << "added " << count << " y_{ki} variables" << std::endl;

	// decision variables u_{i}
	u = IloArray<IloNumVar>(env, inst.n[0]);

	for (int i : inst.B[0]) {
		u[i] = IloNumVar(env, inst.e[i], inst.l[i], ILOINT);
	}

	MLB_OUT(TRACE) << "added " << inst.n[0] << " u_{i} variables" << std::endl;

}

void MLBPTWNFFormulation::addConstraints(IloEnv env, IloModel model, const Instance<MLBPTW>& inst)
{
	// if item i is packed in bin j at level k, and bin j is assigned to bin l at level k + 1, then item i is assigned to bin j at level k + 1
	int count = 0;
	for (int k = 1; k < inst.m; k++) {
		for (int i : inst.B[0]) {
			for (int j : inst.B[k]) {
				for (int l : inst.B[k + 1]) {
					model.add(ib[k][i][j] + x[k][j][l] <= 1 + ib[k + 1][i][l]);
					count++;
				}
			}
		}
	}
	MLB_OUT(TRACE) << "added " << count << " constraints for transitivity between x_{kij} and ib_{kij}" << std::endl;


	// if two items with overlapping time windows are packed into the same top level bin the earliest packing time is as big as the latest starting time between those items
	count = 0;
	for (int a : inst.B[0]) {
		for (int b = a + 1; b < inst.n[0]; b++) {
			for (int top : inst.B[inst.m]) {
				if (inst.l[a] < inst.e[b] || inst.e[a] > inst.l[b]) {
					// Non-overlapping time windows cannot be in the same bin
					model.add(IloIfThen(env, ib[inst.m][a][top] >= 0.5, ib[inst.m][b][top] <= 0.5));
				}
				else if (inst.e[a] > inst.e[b]) {
					model.add(IloIfThen(env, ib[inst.m][a][top] >= 0.5 && ib[inst.m][b][top] >= 0.5, u[b] >= u[a]));
				}
				else if (inst.e[b] > inst.e[a]) {
					model.add(IloIfThen(env, ib[inst.m][a][top] >= 0.5 && ib[inst.m][b][top] >= 0.5, u[a] >= u[b]));
				}
				count++;
			}
		}
	}
	MLB_OUT(TRACE) << "added " << count << " constraints to enforce only items with overlapping time windows can be packed together and their earliest packing time coincides" << std::endl;

	// each item can only be assigned to 1 bin at each level
	for (int k = 0; k <= inst.m; k++) {
		for (int i : inst.B[0]) {
			IloExpr sum(env);
			for (int j : inst.B[k]) {
				sum += ib[k][i][j];
			}
			model.add(sum == 1);
			sum.end();
		}
	}
	MLB_OUT(TRACE) << "added " << ((inst.m + 1) * inst.n[0]) << " constraints such that each bin can be assigned at most to 1 bin at each level" << std::endl;

	// each ib at level 0 is assigned to the same as x
	for (int i : inst.B[0]) {
		for (int j : inst.B[1]) {
			model.add(ib[1][i][j] >= x[0][i][j]);
		}
	}
	MLB_OUT(TRACE) << "added " << inst.n[0] * inst.n[1] << " constraints to enforce ib to be the same as x at level 0" << std::endl;


	/*******************************************************************************/
	/** Constraints from Multi-Level Bin Packing Problem Network Flow formulation **/
	/*******************************************************************************/

	// there can only be flow between 2 item/bins if the lower level bin is assigned to the higher level bin
	count = 0;
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

void MLBPTWNFFormulation::addObjectiveFunction(IloEnv env, IloModel model, const Instance<MLBPTW>& inst)
{
	IloExpr sum(env);
	for (int k : inst.M) {
		for (int j : inst.B[k]) {
			sum += y[k][j] * inst.c[k][j];
		}
	}
	for (int i : inst.B[0]) {
		sum += inst.p * (u[i] - inst.e[i]);
	}
	model.add(IloMinimize(env, sum));
	sum.end();
}

void MLBPTWNFFormulation::extractSolution(IloCplex cplex, const Instance<MLBPTW>& inst, Solution<MLBPTW>& sol)
{
	sol.total_cost = 0;
	for (int k : inst.M) {
		for (int j : inst.B[k]) {
			if (cplex.getValue(y[k][j]) > 0.5)
				sol.total_cost += inst.c[k][j];
		}
	}
	for (int i : inst.B[0]) {
		sol.total_cost += inst.p * (cplex.getValue(u[i]) - inst.e[i]);
	}

	for (int k = 0; k < inst.m; k++) {
		sol.item_to_bins[k].assign(inst.n[k], -1);

		for (int i : inst.B[k])
			for (int j : inst.B[k + 1]) {
				if (cplex.getValue(x[k][i][j]) > 0.5)
					sol.item_to_bins[k][i] = j;
			}

	}
}
