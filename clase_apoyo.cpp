//
// Created by mariwogr on 14/11/21.
//


#include <chrono>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <vector>
#include <limits>

int main() {
    using namespace std::chrono;
    using clk = high_resolution_clock;

    auto start = clk::now();

    constexpr long nsteps = 10'000'000;
    double deltax = 1.0 / nsteps;

    int n=0;
#pragma omp parallel
    {
        n = omp_get_num_threads();
    }

    double sum[n];
    for (int c=0 ; c<n; c++)
        sum[c]=0;

#pragma omp parallel
    {
        long i;
        int t = omp_get_thread_num();

        for (i = t; i < nsteps; i=i+n) {

            double x = (static_cast<double>(i) + 0.5) * deltax;
            sum[t] += 1.0/ (1.0 + x * x);

        }
#pragma omp barrier
    }

    double pi_value = 0.0;

    for (int i = 0 ; i<n; i++){
        //std::cout << i << " " << sum[i] << std::endl;
        pi_value += 4.0 * deltax * sum[i];
    }

    auto stop = clk::now();
    auto diff = duration_cast<microseconds>(stop - start);

    std::cout << "pi_value = " << std::setprecision(10) << pi_value << "\n";
    std::cout << "Time = " << diff.count() << " us\n";
}