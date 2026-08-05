#include <Rcpp.h>
using namespace Rcpp;

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

// [[Rcpp::export]]
double peak_rss_mb() {
#if defined(_WIN32)
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc))) {
        return (double) pmc.PeakWorkingSetSize / (1024.0 * 1024.0);
    }
    return NA_REAL;
#elif defined(__unix__) || defined(__APPLE__)
    struct rusage ru;
    getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    return (double) ru.ru_maxrss / (1024.0 * 1024.0);
#else
    return (double) ru.ru_maxrss / 1024.0;
#endif
#else
    return NA_REAL;
#endif
}
