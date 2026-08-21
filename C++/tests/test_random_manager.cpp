// test_random_manager.cpp — unit tests for inverse CDFs and LHS stratification

#include "random_manager.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>

static int failures = 0;

#define CHECK_NEAR(actual, expected, tol, label) \
    do { \
        if (std::fabs((actual) - (expected)) > (tol)) { \
            printf("FAIL: %s = %.12f (expected %.12f, tol %.2e)\n", \
                   label, (double)(actual), (double)(expected), (double)(tol)); \
            failures++; \
        } \
    } while(0)

int main() {
    // ---- normal_ppf ----
    CHECK_NEAR(normal_ppf(0.5),  0.0,        1e-9,  "normal_ppf(0.5)");
    CHECK_NEAR(normal_ppf(0.025), -1.959963985, 1e-6, "normal_ppf(0.025)");
    CHECK_NEAR(normal_ppf(0.975),  1.959963985, 1e-6, "normal_ppf(0.975)");
    CHECK_NEAR(normal_ppf(0.001), -3.090232306, 1e-6, "normal_ppf(0.001)");
    CHECK_NEAR(normal_ppf(0.999),  3.090232306, 1e-6, "normal_ppf(0.999)");
    CHECK_NEAR(normal_ppf(0.25),  -0.674489750, 1e-6, "normal_ppf(0.25)");
    CHECK_NEAR(normal_ppf(0.75),   0.674489750, 1e-6, "normal_ppf(0.75)");

    // ---- lognormal_ppf ----
    // scipy: lognorm.ppf(0.5, s=0.1, scale=exp(0.3)) = exp(0.3) = 1.34986
    double ln = lognormal_ppf(0.5, 0.3, 0.1);
    CHECK_NEAR(ln, std::exp(0.3), 1e-6, "lognormal_ppf(0.5,0.3,0.1)");

    // ---- uniform_ppf ----
    CHECK_NEAR(uniform_ppf(0.3, 2.0, 5.0), 2.9, 1e-12, "uniform_ppf(0.3,2,5)");

    // ---- log10normal_ppf ----
    // 10^normal_ppf(0.5, 7.0, 1.4) = 10^7.0 = 1e7
    double l1n = log10normal_ppf(0.5, 7.0, 1.4);
    CHECK_NEAR(l1n, 1e7, 1e-3, "log10normal_ppf(0.5,7,1.4)");

    // ---- beta_ppf ----
    // scipy.stats.beta.ppf(0.5, 5, 2) = 0.7355500167043401
    CHECK_NEAR(beta_ppf(0.5, 5.0, 2.0), 0.735550016704340, 1e-8, "beta_ppf(0.5,5,2)");
    // scipy.stats.beta.ppf(0.1, 5, 2) = 0.4896836934485084
    CHECK_NEAR(beta_ppf(0.1, 5.0, 2.0), 0.489683693448508, 1e-8, "beta_ppf(0.1,5,2)");
    // beta with loc/scale: beta_ppf(0.5,5,2)*3 + 2 = 4.20665005011302
    double be = beta_ppf(0.5, 5.0, 2.0) * 3.0 + 2.0;
    CHECK_NEAR(be, 4.20665005011302, 1e-6, "beta with loc/scale");

    // ---- binomial_ppf ----
    // Binomial(10, 0.3): P(X<=2) = 0.38278, P(X<=3) = 0.64961
    // binomial_ppf(0.4, 10, 0.3) = 3 (first k where cdf >= 0.4)
    int bk = binomial_ppf(0.4, 10, 0.3);
    if (bk != 3) { printf("FAIL: binomial_ppf(0.4,10,0.3) = %d (expected 3)\n", bk); failures++; }
    // Edge cases
    if (binomial_ppf(0.0, 10, 0.3) != 0) { printf("FAIL: binomial_ppf(0,...) != 0\n"); failures++; }
    if (binomial_ppf(1.0, 10, 0.3) != 10) { printf("FAIL: binomial_ppf(1,...) != 10\n"); failures++; }

    // ---- LHS stratification ----
    // Draw 10000 from 'uniform' buffer, check stratification
    {
        RandomNumberManager mgr(10000);
        std::vector<int> bins(10, 0);
        for (int i = 0; i < 10000; i++) {
            double v = mgr.get("uniform");
            int b = std::min(static_cast<int>(v * 10), 9);
            bins[b]++;
        }
        for (int b = 0; b < 10; b++) {
            if (bins[b] < 800 || bins[b] > 1200) {
                printf("FAIL: LHS bin %d has %d (expected ~1000)\n", b, bins[b]);
                failures++;
            }
        }
    }

    // ---- LHS stratification for normal buffer ----
    // With LHS over 10000 draws, the normal buffer should also have
    // even stratification — check the transformed values fall in 10
    // roughly equal deciles of the normal distribution.
    {
        RandomNumberManager mgr(10000);
        std::vector<int> bins(10, 0);
        for (int i = 0; i < 10000; i++) {
            double v = random_normal_lhs(mgr, 0.0, 1.0);
            // Normal CDF of the decile boundaries
            // -inf, -1.2816, -0.8416, -0.5244, -0.2533, 0, 0.2533, 0.5244, 0.8416, 1.2816, +inf
            static const double bounds[] = {
                -1e9, -1.2816, -0.8416, -0.5244, -0.2533,
                0.0, 0.2533, 0.5244, 0.8416, 1.2816, 1e9
            };
            for (int b = 0; b < 10; b++) {
                if (v >= bounds[b] && v < bounds[b+1]) {
                    bins[b]++;
                    break;
                }
            }
        }
        for (int b = 0; b < 10; b++) {
            if (bins[b] < 700 || bins[b] > 1300) {
                printf("FAIL: LHS-normal bin %d has %d (expected ~1000)\n", b, bins[b]);
                failures++;
            }
        }
    }

    // ---- Block sizing uses the supplied N, not the 1,000,000 default ----
    // With a tiny batch size (N=5), the LHS stratification is over exactly
    // 5 equal-width bins. If the default (1,000,000) batch were used
    // instead, 5 draws landing one-per-bin across [0,1/5),[1/5,2/5),...
    // would be astronomically unlikely (chance ~5!/5^5 ≈ 0.038, and even
    // then the values wouldn't concentrate near each bin's LHS stratum).
    {
        const int N = 5;
        RandomNumberManager mgr(N);
        std::vector<int> seen(N, 0);
        for (int i = 0; i < N; i++) {
            double v = mgr.get("uniform");
            int b = std::min(static_cast<int>(v * N), N - 1);
            seen[b]++;
        }
        for (int b = 0; b < N; b++) {
            if (seen[b] != 1) {
                printf("FAIL: block-sizing bin %d has %d draws (expected exactly 1 "
                       "for batch_size=N=%d)\n", b, seen[b], N);
                failures++;
            }
        }
    }

    // ---- Independent per-distribution buffers ----
    // Fully exhausting the 'uniform' buffer must not disturb the
    // independently-constructed 'beta' buffer's own stratification.
    {
        const int N = 4;
        RandomNumberManager mgr(N);
        for (int i = 0; i < N; i++) mgr.get("uniform"); // exhaust 'uniform'

        std::vector<int> beta_bins(N, 0);
        for (int i = 0; i < N; i++) {
            double v = mgr.get("beta"); // raw LHS draw in [0,1), untouched buffer
            int b = std::min(static_cast<int>(v * N), N - 1);
            beta_bins[b]++;
        }
        for (int b = 0; b < N; b++) {
            if (beta_bins[b] != 1) {
                printf("FAIL: independent-buffer 'beta' bin %d has %d draws "
                       "(expected 1; exhausting 'uniform' should not affect it)\n",
                       b, beta_bins[b]);
                failures++;
            }
        }
    }

    // ---- Refill after a distribution buffer is exhausted ----
    // Drawing batch_size+1 values from the same distribution forces a
    // refill; the post-refill batch must again be a valid, fully
    // stratified batch of size N (not empty, truncated, or garbage).
    {
        const int N = 6;
        RandomNumberManager mgr(N);
        for (int i = 0; i < N; i++) mgr.get("normal"); // exhaust first batch

        std::vector<int> post_refill_bins(N, 0);
        for (int i = 0; i < N; i++) {
            double v = mgr.get("normal"); // triggers refill on first call, then reads it
            if (v < 0.0 || v > 1.0) {
                printf("FAIL: post-refill draw %d out of [0,1] range: %.6f\n", i, v);
                failures++;
                continue;
            }
            int b = std::min(static_cast<int>(v * N), N - 1);
            post_refill_bins[b]++;
        }
        for (int b = 0; b < N; b++) {
            if (post_refill_bins[b] != 1) {
                printf("FAIL: post-refill bin %d has %d draws (expected exactly 1)\n",
                       b, post_refill_bins[b]);
                failures++;
            }
        }
    }

    if (failures == 0) {
        printf("PASS: all random_manager tests (%d checks)\n", 18);
        return 0;
    } else {
        printf("FAIL: %d test(s) failed\n", failures);
        return 1;
    }
}