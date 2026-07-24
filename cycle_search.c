/*
 * コラッツ周期点 全数探索プログラム
 * =================================
 *
 * 長さ q のすべてのパリティベクトル（2^q通り）について、
 * 「そのビット列を強制的に24回...ではなくq回コラッツ操作した時、
 *   自分自身に戻ってくる有理数 z」を求め、
 * さらに z の"本当の"（動的に決まる）偶奇パターンが、
 * 最初に仮定したビット列と一致するか（自己無撞着かどうか）をチェックする。
 *
 * 一致していて、かつ z が正の整数なら、それは未知のコラッツサイクルの発見。
 * (現在の数学的知見では、q がよほど大きくない限り絶対に見つからないはずだが、
 *  それを自分の手で計算機的に確認する)
 *
 * ビルド方法:
 *   gcc -O3 -march=native -fopenmp cycle_search.c -o cycle_search -lgmp
 *
 * 実行方法:
 *   ./cycle_search <max_q>
 *   例: ./cycle_search 26
 *
 * 必要なライブラリ: GMP (多倍長演算)
 *   Ubuntu/Debian: sudo apt install libgmp-dev
 *   Mac (Homebrew): brew install gmp   (コンパイル時に -I/opt/homebrew/include -L/opt/homebrew/lib を追加)
 *   Windows: MSYS2 (pacman -S mingw-w64-x86_64-gmp) や WSL の利用を推奨
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <omp.h>
#include <time.h>

// bits[i] (i=0..q-1) に対して、周期点 z = D/(1-C) を計算し、
// 自己無撞着かどうか（本当の偶奇がbitsと一致するか）を判定する。
// 戻り値: 1=自己無撞着, 0=不一致
// out_z に結果の z を格納（自己無撞着だった場合のみ意味を持つ）
static int solve_and_check(const int *bits, int q, mpq_t out_z) {
    mpq_t C, D, tmp;
    mpq_inits(C, D, tmp, NULL);
    mpq_set_ui(C, 1, 1);
    mpq_set_ui(D, 0, 1);

    for (int i = 0; i < q; i++) {
        if (bits[i] == 1) {
            // C = 3*C ; D = 3*D + 1
            mpq_set_ui(tmp, 3, 1);
            mpq_mul(C, C, tmp);
            mpq_mul(D, D, tmp);
            mpq_set_ui(tmp, 1, 1);
            mpq_add(D, D, tmp);
        } else {
            // C = C/2 ; D = D/2
            mpq_set_ui(tmp, 2, 1);
            mpq_div(C, C, tmp);
            mpq_div(D, D, tmp);
        }
    }

    mpq_t one;
    mpq_init(one);
    mpq_set_ui(one, 1, 1);
    if (mpq_equal(C, one)) {
        mpq_clears(C, D, tmp, one, NULL);
        return 0;
    }

    // z = D / (1 - C)
    mpq_t oneMinusC;
    mpq_init(oneMinusC);
    mpq_sub(oneMinusC, one, C);
    mpq_div(out_z, D, oneMinusC);

    mpq_clears(C, D, tmp, one, oneMinusC, NULL);

    // 自己無撞着性チェック: z の動的な偶奇が bits と一致するか
    mpq_t cur;
    mpq_init(cur);
    mpq_set(cur, out_z);

    int consistent = 1;
    for (int i = 0; i < q; i++) {
        // 分母は常に奇数のはずなので、分子の偶奇だけ見ればよい
        int parity = mpz_odd_p(mpq_numref(cur)) ? 1 : 0;
        if (parity != bits[i]) {
            consistent = 0;
            break;
        }
        if (parity == 0) {
            mpq_t two;
            mpq_init(two);
            mpq_set_ui(two, 2, 1);
            mpq_div(cur, cur, two);
            mpq_clear(two);
        } else {
            mpq_t three, oneL;
            mpq_inits(three, oneL, NULL);
            mpq_set_ui(three, 3, 1);
            mpq_set_ui(oneL, 1, 1);
            mpq_mul(cur, cur, three);
            mpq_add(cur, cur, oneL);
            mpq_clears(three, oneL, NULL);
        }
    }
    if (consistent && !mpq_equal(cur, out_z)) {
        consistent = 0;
    }
    mpq_clear(cur);
    return consistent;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s <max_q> [min_q]\n", argv[0]);
        return 1;
    }
    int max_q = atoi(argv[1]);
    int min_q = (argc >= 3) ? atoi(argv[2]) : 1;

    printf("q からmax_q=%dまで全数探索します（GMP + OpenMP, %d threads）\n",
           max_q, omp_get_max_threads());
    fflush(stdout);

    for (int q = min_q; q <= max_q; q++) {
        long long total = 1LL << q;
        long long consistent_count = 0;
        long long pos_int_hits = 0;

        double t0 = omp_get_wtime();

        #pragma omp parallel
        {
            long long local_consistent = 0;
            long long local_hits = 0;

            #pragma omp for schedule(dynamic, 4096)
            for (long long pattern = 0; pattern < total; pattern++) {
                int bits[64]; // q は64未満を想定
                for (int i = 0; i < q; i++) {
                    bits[i] = (int)((pattern >> i) & 1LL);
                }

                mpq_t z;
                mpq_init(z);
                int ok = solve_and_check(bits, q, z);
                if (ok) {
                    local_consistent++;
                    if (mpz_cmp_ui(mpq_denref(z), 1) == 0 && mpq_sgn(z) > 0) {
                        local_hits++;
                        #pragma omp critical
                        {
                            gmp_printf("  *** HIT *** q=%d pattern=%lld z=%Qd\n", q, pattern, z);
                        }
                    }
                }
                mpq_clear(z);
            }

            #pragma omp atomic
            consistent_count += local_consistent;
            #pragma omp atomic
            pos_int_hits += local_hits;
        }

        double elapsed = omp_get_wtime() - t0;
        printf("q=%2d: self-consistent=%8lld / %10lld total   positive-int hits=%lld   (%.2fs)\n",
               q, consistent_count, total, pos_int_hits, elapsed);
        fflush(stdout);
    }

    return 0;
}
