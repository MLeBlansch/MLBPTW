#include "solution.h"

#include "lib/util.h"
#include "users.h"

#include <algorithm>
#include <numeric>
#include <cassert>
#include <climits>



/*****************************************************************************************/
/** Bin Packing Problem ******************************************************************/
/*****************************************************************************************/
Solution<BP>::Solution(const Instance<BP>& inst) : item_to_bins(inst.s.size())
{
	std::iota(item_to_bins.begin(), item_to_bins.end(), 0);
	total_bins = (int)inst.s.size();
}

std::ostream& operator<<(std::ostream& os, const Solution<BP>& sol)
{
	os << sol.item_to_bins;
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem ******************************************************/
/*****************************************************************************************/
Solution<MLBP>::Solution(const Instance<MLBP>& inst) : db(-1), item_to_bins(inst.m)
{
	item_to_bins.assign(inst.m, std::vector<int>(inst.n[0]));
	for (std::vector<int> lvl : item_to_bins) {
		std::iota(lvl.begin(), lvl.end(), 0);
	}

	total_cost = INT_MAX;
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBP>& sol) {
	os << "[" << std::endl;
	for (std::vector<int> lvl : sol.item_to_bins) {
		os << lvl << std::endl;
	}
	os << "]";
	return os;
}




/*****************************************************************************************/
/** Class Constrained Multi-Level Bin Packing Problem ************************************/
/*****************************************************************************************/
Solution<CCMLBP>::Solution(const Instance<CCMLBP>& inst) : Solution<MLBP>(inst)
{
}

std::ostream& operator<<(std::ostream& os, const Solution<CCMLBP>& sol) {
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Conflict Constraints ****************************/
/*****************************************************************************************/
Solution<MLBPCC>::Solution(const Instance<MLBPCC>& inst) : Solution<MLBP>(inst)
{
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBPCC>& sol) {
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Partial Orders **********************************/
/*****************************************************************************************/
Solution<MLBPPO>::Solution(const Instance<MLBPPO>& inst) : Solution<MLBP>(inst)
{
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBPPO>& sol) {
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Time Windows ************************************/
/*****************************************************************************************/
Solution<MLBPTW>::Solution(const Instance<MLBPTW>& inst) : Solution<MLBP>(inst)
{
	item_to_bins.assign(inst.m, std::vector<int>(inst.n[0]));
	for (std::vector<int> lvl : item_to_bins) {
		std::iota(lvl.begin(), lvl.end(), 0);
	}

	total_cost = INT_MAX;
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBPTW>& sol) {
	os << "[" << std::endl;
	for (std::vector<int> lvl : sol.item_to_bins) {
		os << lvl << std::endl;
	}
	os << "]";
	return os;
}




/*****************************************************************************************/
/** Multi-Level Bin Packing Problem with Fragmentation Constraints ***********************/
/*****************************************************************************************/
Solution<MLBPFC>::Solution(const Instance<MLBPFC>& inst) : Solution<MLBP>(inst)
{
}

std::ostream& operator<<(std::ostream& os, const Solution<MLBPFC>& sol) {
	return os;
}


