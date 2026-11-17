#include "StatisticTools.h"
namespace Statistic
{

double getMaxError(vector<double> &exactBC, vector<double> &approximateBC) {
    double res = 0;
    for (int i = 0; i < (int)exactBC.size(); i++) {
        res = max(res, abs(exactBC[i] - approximateBC[i]));
    }
    return res;
}

double getMaeError(vector<double> &exactBC, vector<double> &approximateBC) {
    double res = 0;
    for (int i = 0; i < (int)exactBC.size(); i++) {
        res += abs(exactBC[i] - approximateBC[i]);
    }
    res /= (double)(exactBC.size());
    return res;
}

double getRmseError(vector<double> &exactBC, vector<double> &approximateBC) {
    double res = 0;
    for (int i = 0; i < (int)exactBC.size(); i++) {
        double delta = abs(exactBC[i] - approximateBC[i]);
        res += delta * delta;
    }
    res /= (double)(exactBC.size());
    res = sqrt(res);
    return res;
}

double getMaxRelativeError(vector<double>& exactBC, vector<double>& approximateBC, double threshold) {
    double res = 0;
    int fnd = 0;
    for (int i = 0; i < sz(exactBC); i++) {
        if (exactBC[i] < threshold || exactBC[i] == 0) continue;
        res = max(res, abs(exactBC[i] - approximateBC[i]) / exactBC[i]);
        fnd = 1;
    }
    if (fnd == 0) return numeric_limits<double>::quiet_NaN();
    return res;
}

double getPrecisionAtK(vector<double>& exactBC, vector<double>& approximateBC,int k) {
    int n = sz(exactBC);
    k = min(k, n);
    vector<int> a(n), b(n);
    iota(all(a), 0);
    iota(all(b), 0);
    auto cmp = [&](int x, int y, const vector<double>& bc) {
        if (bc[x] != bc[y]) return bc[x] > bc[y];
        return x < y;
    };
    sort(all(a), [&](int x, int y) {
        return cmp(x, y, exactBC);
    });
    sort(all(b), [&](int x, int y) {
        return cmp(x, y, approximateBC);
    });
    vector<bool> used(n, false);
    for (int i = 0; i < k; i++) {
        used[a[i]] = true;
    }
    int cnt = 0;
    for (int i = 0; i < k; i++) {
        cnt += used[b[i]];
    }
    return (double)cnt / k;
}

double getPeakRAM() {
    PROCESS_MEMORY_COUNTERS pmc{};

    GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc));

    return pmc.PeakWorkingSetSize / (1024.0 * 1024.0);
}

void getAllStatistic(string filename, vector<double> &exactBC, vector<double> &approximateBC, double RAM, double time, int k)
{
    cout << fixed << setprecision(6);

    cout << "Dataset             : " << filesystem::path(filename).filename().string() << '\n';
    cout << "\n========== TBC STATISTICS ==========\n";
    cout << "Samples (K)         : " << k << '\n';
    cout << "Execution time (s)  : " << time << '\n';
    cout << "Peak RAM (MiB)      : " << RAM << '\n';

    cout << "\n--------------- ERRORS -------------\n";

    cout << "Max absolute error  : " << getMaxError(exactBC, approximateBC) << '\n';

    cout << "MAE                 : " << getMaeError(exactBC, approximateBC) << '\n';

    cout << "RMSE                : "
         << getRmseError(exactBC, approximateBC) << '\n';

    cout << "Max relative error  : "
         << getMaxRelativeError(exactBC, approximateBC, 1e-3) * 100.0
         << "%\n";

    cout << "\n-------------- TOP-K ---------------\n";

    for (int topK : {10, 25, 50, 100, 250, 500, 1000, 2000}) {
        cout << "Top K precision: " << left << setw(4) << topK
             << right << "       : "
             << getPrecisionAtK(exactBC, approximateBC, topK)
             << '\n';
    }

    cout << "====================================\n";
}

}