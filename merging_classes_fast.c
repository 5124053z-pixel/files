/*
 * merging_classes_fast.c
 * =======================
 * n ≡ r (mod 2^k) のとき、nとn+1が「同じステップ数で、同じ値に合流する」
 * という現象が、rによらず"必ず"起きる剰余類(merging class)の個数を数える。
 *
 * k を大きくしていったとき、merging classの割合が何%に収束するかを
 * 調べることで、「隣接整数の約50%が同じ全停止ステップ数を持つ」という
 * 既知だが未証明の現象のうち、どれだけが「軌道の完全な合流」という
 * 単純なメカニズムで説明できるかを見る。
 *
 * ビルド:
 *   gcc -O3 -march=native -fopenmp merging_classes_fast.c -o merging_classes_fast
 *
 * 実行:
 *   ./merging_classes_fast <k_min> <k_max>
 *   例: ./merging_classes_fast 3 24
 *
 * 注意: modulus (2^k) が大きくなるほど、途中で使う整数も大きくなります。
 * unsigned __int128 を使っているので、k=30程度までは安全のはずです。
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

typedef unsigned __int128 u128;

// コラッツ軌道を辿り、traj[]に値を記録する。1に到達するか、
// max_steps に達したら終了。実際に格納した要素数を返す。
static int collatz_traj(u128 n, u128 *traj, int max_steps) {
    int len = 0;
    traj[len++] = n;
    while (n != 1 && len < max_steps) {
        if (n % 2 == 0) n = n / 2;
        else n = 3 * n + 1;
        traj[len++] = n;
    }
    return len;
}

// nの軌道の中に、n1の軌道と同じ値がいつ最初に現れるかを探す。
// (nから何ステップ目, n1から何ステップ目) を merge_a, merge_b に書き込む。
// 見つからなければ 0 を返す。
static int find_merge(u128 n, u128 n1, int max_steps, int *merge_a, int *merge_b) {
    static _Thread_local u128 traj_n[300];
    static _Thread_local u128 traj_n1[300];
    int len_n  = collatz_traj(n,  traj_n,  max_steps);
    int len_n1 = collatz_traj(n1, traj_n1, max_steps);

    for (int i = 0; i < len_n; i++) {
        for (int j = 0; j < len_n1; j++) {
            if (traj_n[i] == traj_n1[j]) {
                *merge_a = i;
                *merge_b = j;
                return 1;
            }
        }
    }
    return 0;
}

// r (mod modulus) が merging class かどうか判定。
// num_samples 個の代表元(j=2,3,...)について、常に同じ(merge_a, merge_b)で
// 合流するかを確認する。
static int is_merging_class(u128 r, u128 modulus, int num_samples, int max_steps) {
    int first_a = -1, first_b = -1;
    for (int s = 0; s < num_samples; s++) {
        u128 j = (u128)(s + 2);
        u128 n = r + j * modulus;
        if (n < 2) continue;
        int a, b;
        if (!find_merge(n, n + 1, max_steps, &a, &b)) {
            return 0;
        }
        if (first_a == -1) {
            first_a = a; first_b = b;
        } else if (a != first_a || b != first_b) {
            return 0;
        }
    }
    return 1;
}

int main(int argc, char **argv) {
    if (argc < 3) {
        fprintf(stderr, "usage: %s <k_min> <k_max>\n", argv[0]);
        return 1;
    }
    int k_min = atoi(argv[1]);
    int k_max = atoi(argv[2]);
    int num_samples = 4;
    int max_steps = 260;

    printf("k_min=%d, k_max=%d, threads=%d\n", k_min, k_max, omp_get_max_threads());
    fflush(stdout);

    for (int k = k_min; k <= k_max; k++) {
        u128 modulus = (u128)1 << k;
        long long total = (long long)1 << k;
        long long count = 0;

        double t0 = omp_get_wtime();

        #pragma omp parallel for schedule(dynamic, 1024) reduction(+:count)
        for (long long r = 0; r < total; r++) {
            if (is_merging_class((u128)r, modulus, num_samples, max_steps)) {
                count++;
            }
        }

        double elapsed = omp_get_wtime() - t0;
        double frac = (double)count / (double)total;
        printf("k=%2d modulus=%lld : merging=%lld/%lld = %.4f%%   (%.1fs)\n",
               k, (long long)modulus, count, total, 100.0 * frac, elapsed);
        fflush(stdout);
    }

    return 0;
}
