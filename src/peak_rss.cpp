#include <Rcpp.h>
using namespace Rcpp;

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>
#elif defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

// Goal:   Ask the operating system for the most memory (RAM) this process has
//         ever used. It is a running maximum, not the amount in use right now,
//         so freeing memory does not lower it. Takes no arguments.
//
// Output: One number: maximum memory used, in MB. NA if the platform is not
//         supported, or if the Windows query fails.
//
// The three branches below do the same thing with different OS calls; the unit
// each one returns differs (Windows and macOS give bytes, Linux gives kB),
// which is why the divisors are not identical.
//
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
