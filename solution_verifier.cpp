#include "solution_verifier.h"
#include "instance.h"
#include "solution.h"
#include "users.h"

#include <sstream>
#include <algorithm>
#include <set>
#include <climits>



/*************************************************************************************************/
/* BP ********************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<BP>::verify(const Instance<BP>& inst, const Solution<BP>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	std::vector<int> levels(sol.total_bins, 0);  // size of each bin
	for (int i : inst.I) {
		int bin = sol.item_to_bins[i];
		if (bin >= sol.total_bins) {
			if (error_msg) {
				std::stringstream ss;
				ss << "Detected bin with id " << bin << " but solutions claims that there are only " << sol.total_bins << ".";
				error_msg->push_back(ss.str());
			}
			ret = false;
		} else if (bin < 0) {
			if (error_msg) {
				std::stringstream ss;
				ss << "Item " << i << " is not assigned to any bin.";
				error_msg->push_back(ss.str());
			}
			ret = false;
		} else {
			levels[bin] += inst.s[i];
		}
	}
	for (int bin = 0; bin < sol.total_bins; bin++) {
		if (levels[bin] > inst.smax) {
			if (error_msg) {
				std::stringstream ss;
				ss << "Bin " << bin << " of size " << levels[bin] << " exceeds maximum capaicty (" << inst.smax << ").";
				error_msg->push_back(ss.str());
			}
			ret = false;
		}
	}
	return ret;
}


/*************************************************************************************************/
/* MLBP ******************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBP>::verify(const Instance<MLBP>& inst, const Solution<MLBP>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;

	std::vector<std::vector<int>> usedBins(inst.m);
	for (int k = 0; k < inst.m; k++) {
		usedBins[k] = std::vector<int>(inst.n[k+1], -1);
	}

	for (int k = 0; k < inst.m; k++) {
		std::vector<int> capacity(inst.n[k + 1], 0);
		for (int i : inst.B[k]) {
			int bin = sol.item_to_bins[k][i];
			if (bin >= inst.n[k + 1]) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Detected bin with id " << bin << " but there are only " << inst.n[k + 1] << " bins at level " << k << ".";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (k == 0 && bin < 0) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Item " << i << " is not assigned to any bin.";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (k > 0 && bin < 0 && usedBins[k-1][i] > 0) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Bin with id " << i << " of level " << k << " is used but not assigned to any bin.";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (k > 0 && bin > 0 && usedBins[k-1][i] < 0) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Bin with id " << i << " of level " << k << " is assigned to a bin but not used.";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (bin >= 0) {
				capacity[bin] += inst.s[k][i];
				usedBins[k][bin] = 1;
			}
		}
		for (int bin = 0; bin < inst.n[k + 1]; bin++) {
			if (capacity[bin] > inst.c[k+1][bin]) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Bin " << bin << " of level " << k << " with capacity " << inst.c[k][bin] << " exceeds maximum capacity (" << capacity[bin] << ").";
					MLB_OUT(TRACE) << ss.str();
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
		}
	}
	usedBins.clear();
	return ret;
}





/*************************************************************************************************/
/* CCMLBP ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<CCMLBP>::verify(const Instance<CCMLBP>& inst, const Solution<CCMLBP>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	return ret;
}




/*************************************************************************************************/
/* MLBPCC ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBPCC>::verify(const Instance<MLBPCC>& inst, const Solution<MLBPCC>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	return ret;
}




/*************************************************************************************************/
/* MLBPPO ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBPPO>::verify(const Instance<MLBPPO>& inst, const Solution<MLBPPO>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	return ret;
}




/*************************************************************************************************/
/* MLBPTW ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBPTW>::verify(const Instance<MLBPTW>& inst, const Solution<MLBPTW>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;

	std::vector<std::vector<int>> usedBins(inst.m);
	std::vector<std::vector<int>> e(inst.m + 1);
	std::vector<std::vector<int>> l(inst.m + 1);
	e[0] = inst.e;
	l[0] = inst.l;

	for (int k = 0; k < inst.m; k++) {
		usedBins[k] = std::vector<int>(inst.n[k + 1], -1);
		
		e[k+1] = std::vector<int>(inst.n[k+1], 0);
		l[k+1] = std::vector<int>(inst.n[k+1], INT_MAX);
	}

	for (int k = 0; k < inst.m; k++) {
		std::vector<int> capacity(inst.n[k + 1], 0);
		for (int i : inst.B[k]) {
			int bin = sol.item_to_bins[k][i];
			if (bin >= inst.n[k + 1]) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Detected bin with id " << bin << " but there are only " << inst.n[k + 1] << " bins at level " << k << ".";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (k == 0 && bin < 0) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Item " << i << " is not assigned to any bin.";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (k > 0 && bin < 0 && usedBins[k - 1][i] > 0) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Bin with id " << i << " of level " << k << " is used but not assigned to any bin.";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (k > 0 && bin > 0 && usedBins[k - 1][i] < 0) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Bin with id " << i << " of level " << k << " is assigned to a bin but not used.";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (bin >= 0 && (e[k][i] > l[k+1][bin] || l[k][i] < e[k+1][bin])) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Item/bin " << i << " of level " << k << " with interval [" << e[k][i] << "," << l[k][i] << "] is assigned to bin " << bin << " with interval[" << e[k+1][bin] << ", " << l[k+1][bin] << "] at level " << k << " .";
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
			else if (bin >= 0) {
				capacity[bin] += inst.s[k][i];
				usedBins[k][bin] = 1;
				e[k + 1][bin] = std::max(e[k + 1][bin], e[k][i]);
				l[k + 1][bin] = std::min(l[k + 1][bin], l[k][i]);
			}
		}
		for (int bin = 0; bin < inst.n[k + 1]; bin++) {
			if (capacity[bin] > inst.c[k + 1][bin]) {
				if (error_msg) {
					std::stringstream ss;
					ss << "Bin " << bin << " of level " << k << " with capacity " << inst.c[k][bin] << " exceeds maximum capacity (" << capacity[bin] << ").";
					MLB_OUT(TRACE) << ss.str();
					error_msg->push_back(ss.str());
				}
				ret = false;
			}
		}
	}
	usedBins.clear();
	return ret;
}




/*************************************************************************************************/
/* MLBPFC ****************************************************************************************/
/*************************************************************************************************/
bool SolutionVerifier<MLBPFC>::verify(const Instance<MLBPFC>& inst, const Solution<MLBPFC>& sol, std::vector<std::string>* error_msg)
{
	bool ret = true;
	return ret;
}
