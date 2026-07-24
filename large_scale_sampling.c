/*
 * large_scale_sampling.c
 * =======================
 * ランダムな大きい整数 n (指定したビット数)について、
 * steps(n) == steps(n+1) となる頻度をサンプリングで測定する。
 *
 * 目的: 「隣接整数の一致率が約50%」という既知の(だが未証明の)事実が、
 * nを大きくするにつれてどこまでドリフトし続けるかを、
 * 全数チェックではなくランダムサンプリングで大きい桁数まで確認する。
 *
 * ビルド:
 *   gcc -O3 -march=native -fopenmp large_scale_sampling.c -o large_scale_sampling -lgmp
 *
 * 実行:
 *   ./large_scale_sampling <bit_min> <bit_max> <num_samples>
 *   例: ./large_scale_sampling 10 10000 2000
 *       -> ビット数を10,20,50,100,200,500,1000,2000,5000,10000...と
 *          対数的に増やしながら、各段階で2000サンプルずつ測定する
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <omp.h>
#include <time.h>

static long collatz_steps(mpz_t n_in) {
    mpz_t n;
    mpz_init_set(n, n_in);
    long steps = 0;
    while (mpz_cmp_ui(n, 1) != 0) {
        if (mpz_even_p(n)) {
            mpz_fdiv_q_2exp(n, n, 1);
        } else {
            mpz_mul_ui(n, n, 3);
            mpz_add_ui(n, n, 1);
        }
        steps++;
    }
    mpz_clear(n);
    return steps;
}

// 指定したビット数のランダムな奇数を1つ生成 (再現性のため、呼び出し毎にseedを変える)
static void random_n(mpz_t out, gmp_randstate_t state, unsigned long bits) {
    mpz_urandomb(out, state, bits);
    // 最上位ビットを立てて、確実に指定ビット数にする
    mpz_setbit(out, bits - 1);
}

int main(int argc, char **argv) {
    if (argc < 4) {
        fprintf(stderr, "usage: %s <bit_min> <bit_max> <num_samples_per_level>\n", argv[0]);
        return 1;
    }
    unsigned long bit_min = strtoul(argv[1], NULL, 10);
    unsigned long bit_max = strtoul(argv[2], NULL, 10);
    long num_samples = atol(argv[3]);

    printf("threads=%d\n", omp_get_max_threads());
    fflush(stdout);

    // ビット数を対数的に(だいたい2倍刻みで)増やしていくリストを作る
    unsigned long bit_levels[128];
    int num_levels = 0;
    unsigned long b = bit_min;
    while (b <= bit_max && num_levels < 128) {
        bit_levels[num_levels++] = b;
        b = (unsigned long)(b * 1.7) + 1;
    }

    for (int li = 0; li < num_levels; li++) {
        unsigned long bits = bit_levels[li];
        long match = 0;

        double t0 = omp_get_wtime();

        #pragma omp parallel reduction(+:match)
        {
            gmp_randstate_t state;
            gmp_randinit_mt(state);
            unsigned long seed = (unsigned long)time(NULL) * 1000003UL
                                  + (unsigned long)omp_get_thread_num() * 7919UL
                                  + bits;
            gmp_randseed_ui(state, seed);

            mpz_t n, n1;
            mpz_inits(n, n1, NULL);

            #pragma omp for schedule(dynamic, 8)
            for (long s = 0; s < num_samples; s++) {
                random_n(n, state, bits);
                if (mpz_even_p(n)) mpz_add_ui(n, n, 1); // 奇数にそろえる(必須ではないが安定させる)
                mpz_add_ui(n1, n, 1);

                long steps_n  = collatz_steps(n);
                long steps_n1 = collatz_steps(n1);
                if (steps_n == steps_n1) match++;
            }

            mpz_clears(n, n1, NULL);
            gmp_randclear(state);
        }

        double elapsed = omp_get_wtime() - t0;
        double rate = (double)match / (double)num_samples;
        printf("bits=%6lu (~10^%.0f) : match=%5ld/%5ld = %.2f%%   (%.1fs)\n",
               bits, bits * 0.30103, match, num_samples, 100.0 * rate, elapsed);
        fflush(stdout);
    }

    return 0;
}
