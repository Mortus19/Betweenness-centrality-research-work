#include <bits/stdc++.h>
#include <windows.h>
#include <psapi.h>
using namespace std;

#define all(x) x.begin(), x.end()
#define sz(x) (int)x.size()
#define debug(x) cerr << __LINE__ << ": " << #x << " = " << (x) << endl

namespace Statistic
{
double getMaxError(vector<double> &exactBC, vector<double> &approximateBC);
double getMaeError(vector<double> &exactBC, vector<double> &approximateBC);
double getRmseError(vector<double> &exactBC, vector<double> &approximateBC);
double getMaxRelativeError(vector<double>& exactBC, vector<double>& approximateBC, double threshold = 1e-3);
double getPrecisionAtK(vector<double>& exactBC, vector<double>& approximateBC,int k);
double getPeakRAM();
void getAllStatistic(string filename, string nameAlgo, int threads, vector<double> &exactBC, vector<double> &approximateBC, double RAM, double time, int k);
void getAllStatistic(string filename, string nameAlgo, int threads, double RAM, double time);
}

