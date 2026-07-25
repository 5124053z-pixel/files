/*
 * coupling_scaling.c
 *
 * For a sequence of bit-length windows, draw random integers n of that
 * bit length, and for pairs (n, n+1):
 *   1. compute total stopping time steps(n), steps(n+1) (Collatz, GMP)
 *   2. if they agree, walk both trajectories in lock-step (raw 3x+1 map,
 *      no special-casing at 1) and record the first t < s at which the
 *      values literally coincide ("early merge" / coupling time tau)
 *
 * Output (per window): bit length, samples, #agree, #early-merge,
 * and the full list of merge times (so the decay rate gamma can be
 * fit in post-processing, e.g. python, per window).
 *
 * This directly tests whether the decay rate gamma(N) of
 * P(tau_couple > t | merged) drifts with scale (window bit length),
 * which is the open question motivating this experiment.
 *
 * Build:
 *   gcc -O3 -fopenmp -o coupling_scaling coupling_scaling.c -lgmp
 *
 * Run:
 *   ./coupling_scaling <samples_per_window> <out_prefix>
 *
 * Example:
 *   ./coupling_scaling 20000 run1
 * produces run1_bits<L>.csv for each window bit length L, one line
 * per sample: "agree,merge_time" (merge_time = -1 if no early merge
 * found, i.e. they only coincide trivially at t=s, or -2 if steps
 * disagree).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gmp.h>
#include <omp.h>

/* Bit-length windows to test. Extend as desired -- GMP makes this
 * scale-free, cost grows with number of steps (~7*bits on average),
 * not with a fixed machine-word limit. */
static const long BIT_WINDOWS[] = {
    32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536
};
static const int NUM_WINDOWS = sizeof(BIT_WINDOWS) / sizeof(BIT_WINDOWS[0]);

/* Cap on total steps we're willing to walk for early-merge search,
 * as a safety valve against pathological trajectories. Generous
 * multiple of the expected ~7*bits average. */
#define STEP_CAP_MULT 200

/* total stopping time: steps until x == 1. Modifies x in place. */
static long total_stopping_time(mpz_t x) {
    long t = 0;
    while (mpz_cmp_ui(x, 1) != 0) {
        if (mpz_even_p(x)) {
            mpz_fdiv_q_2exp(x, x, 1);
        } else {
            mpz_mul_ui(x, x, 3);
            mpz_add_ui(x, x, 1);
        }
        t++;
        if (t > 100000000L) { /* pathological safety valve */
            return -1;
        }
    }
    return t;
}

/*
 * Search for the first t in [1, s) such that T^t(n) == T^t(n+1),
 * applying the raw map (no special-casing at 1). Returns t if found,
 * -1 if not found within [1, s).
 */
static long early_merge_time(const mpz_t n, long s) {
    mpz_t x, y;
    mpz_init_set(x, n);
    mpz_init(y);
    mpz_add_ui(y, n, 1);

    long result = -1;
    for (long t = 1; t < s; t++) {
        if (mpz_even_p(x)) mpz_fdiv_q_2exp(x, x, 1);
        else { mpz_mul_ui(x, x, 3); mpz_add_ui(x, x, 1); }

        if (mpz_even_p(y)) mpz_fdiv_q_2exp(y, y, 1);
        else { mpz_mul_ui(y, y, 3); mpz_add_ui(y, y, 1); }

        if (mpz_cmp(x, y) == 0) { result = t; break; }
    }
    mpz_clear(x);
    mpz_clear(y);
    return result;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <samples_per_window> <out_prefix>\n", argv[0]);
        return 1;
    }
    long M = atol(argv[1]);
    const char *prefix = argv[2];

    for (int w = 0; w < NUM_WINDOWS; w++) {
        long bits = BIT_WINDOWS[w];
        char fname[512];
        snprintf(fname, sizeof(fname), "%s_bits%ld.csv", prefix, bits);
        FILE *out = fopen(fname, "w");
        if (!out) { perror("fopen"); return 1; }
        fprintf(out, "agree,merge_time\n");

        long agree_count = 0, merge_count = 0;
        double t_start = omp_get_wtime();

        #pragma omp parallel
        {
            gmp_randstate_t rs;
            gmp_randinit_mt(rs);
            unsigned long seed = (unsigned long)time(NULL) ^ (unsigned long)omp_get_thread_num() * 2654435761UL;
            gmp_randseed_ui(rs, seed);

            mpz_t n, x, y;
            mpz_inits(n, x, y, NULL);

            #pragma omp for schedule(dynamic, 64) reduction(+:agree_count, merge_count)
            for (long i = 0; i < M; i++) {
                /* random n with exactly 'bits' bits: top bit set,
                 * remaining bits-1 random. */
                mpz_urandomb(n, rs, bits - 1);
                mpz_setbit(n, bits - 1); /* ensure exact bit length */
                if (mpz_even_p(n)) mpz_add_ui(n, n, 1); /* prefer odd start, arbitrary choice */

                mpz_set(x, n);
                long sX = total_stopping_time(x);

                mpz_add_ui(y, n, 1);
                long sY = total_stopping_time(y);

                long merge_t = -2; /* sentinel: steps disagree */
                if (sX >= 0 && sY >= 0 && sX == sY) {
                    agree_count++;
                    merge_t = early_merge_time(n, sX);
                    if (merge_t >= 0) merge_count++;
                }

                #pragma omp critical(io)
                {
                    fprintf(out, "%d,%ld\n", (sX == sY && sX >= 0) ? 1 : 0, merge_t);
                }
            }
            mpz_clears(n, x, y, NULL);
            gmp_randclear(rs);
        }

        double elapsed = omp_get_wtime() - t_start;
        fclose(out);

        fprintf(stderr,
            "bits=%6ld  samples=%6ld  agree=%6ld (%.4f)  early_merge=%6ld (of agree: %.4f)  time=%.1fs  -> %s\n",
            bits, M, agree_count, (double)agree_count / M,
            merge_count, agree_count ? (double)merge_count / agree_count : 0.0,
            elapsed, fname);
    }

    return 0;
}
