#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gmp.h>

int index_b(int i, int j){
    int k;

    k = i * (i + 1) / 2 + j;
    
    return k;

}

void bincoef_array(mpz_t b[], int K){
    for (int i = 0; i < K - 1; i++){
        for (int j = 0; j <= i; j++){
            mpz_init(b[index_b(i, j)]);
            mpz_bin_uiui(b[index_b(i, j)], i + 2, j + 1);  
        }
    }
}

void W(mpz_t w, int r, int l, mpz_t bincoef_a []){
    int i, j, k;
    mpz_t w_array [r], sum_w, product, bincoef, powil, final_result;

    mpz_init_set_ui(w_array[0], 1); 

    for (int i = 2; i <= r; ++i){
        mpz_init_set_ui(sum_w, 0); 

        for (int j=1; j < i; ++j){
            mpz_init(product);

            mpz_mul(product, w_array[j - 1], bincoef_a[index_b(i - 2, j - 1)]);

            mpz_add(sum_w, sum_w, product);

            mpz_clear(product);
        }

        mpz_init(powil);

        mpz_ui_pow_ui(powil, i, l);

        mpz_init(w_array[i - 1]);

        mpz_sub(w_array[i - 1], powil, sum_w);

        mpz_clear(powil);
        mpz_clear(sum_w);
        

    }
    mpz_init_set(w, w_array[r - 1]);


    for (int k=0; k < r; k++){
        mpz_clear(w_array[k]); 
    }

}

void pmf(mpf_t p, int n, int M, double pr, double pg, mpz_t b_array []){

    // returns P(X=n)

    int l;

    mpz_t w;

    mpf_t pgg, prg;
    mpf_t prob, sum_prob;

    mpf_init_set_d(pgg, pg);
    mpf_init_set_d(prg, pr);
    
    mpf_init(sum_prob);
    mpf_init(prob);

    mpf_t one_minus_pg, one_minus_pg_pow, pr_pow;
    mpf_t prod_0, prod_1, prod_2, prod_3, prod_4, prod_5;

    mpz_t bin_coef;
    mpf_t bin_coef_f;

    mpf_t w_f;
    
    mpf_init(one_minus_pg);
    mpf_init(one_minus_pg_pow);
    mpf_init(pr_pow);
    mpf_init(prod_0);
    mpf_init(prod_1);
    mpf_init(prod_2);
    mpf_init(prod_3);

    mpz_init(bin_coef);
    mpf_init(bin_coef_f);

    mpf_init(w_f);
    
    mpf_init(p);

    mpf_set_d(prob, 0.0);

    for (l = M - 1; l < n; l++){
        mpf_ui_sub(one_minus_pg, 1, pgg);
        mpf_pow_ui(one_minus_pg_pow,  one_minus_pg, n - l - 1);

        mpf_pow_ui(pr_pow, prg, l + 1);

        mpf_mul(prod_0, one_minus_pg_pow, pr_pow);

        mpf_mul_ui(prod_1, prod_0, M);

        mpz_bin_uiui(bin_coef, n - 1, l);

        mpf_set_z(bin_coef_f, bin_coef);

        mpf_mul(prod_2, prod_1, bin_coef_f);

        W(w, M - 1, l, b_array);

        mpf_set_z(w_f, w);

        mpf_mul(prod_3, prod_2, w_f);

        mpf_add(prob, prob, prod_3);


    }

    mpf_set(p, prob);
        
    mpf_clear(pgg);
    mpf_clear(prg);
    mpf_clear(sum_prob);
    mpf_clear(one_minus_pg);
    mpf_clear(one_minus_pg_pow);
    mpf_clear(pr_pow);
    mpf_clear(prod_0);
    mpf_clear(prod_1);
    mpz_clear(bin_coef);
    mpf_clear(prod_2);
    mpf_clear(w_f);
    mpf_clear(prod_3);

}

int main(int argc, char * argv[]){
    if (argc <= 2){
    printf ("Usage: %s <number M> <number n> \n", argv[0]);
    return 2;
    }

    const int save_frequency = 5;
    int M;
    int n;

    M = atoi(argv[1]); 


    n = atoi(argv[2]);
    assert( n >= 0);

    assert(n >= M);

    int l;

    double pg = 0.95;
    double pr = pg / M;

    char file_name[45];
    sprintf(file_name, "../results/pmf_cdf_%d_%d_%d.csv", M, n, (int) (pg * 1000));

    printf("nombre: %s", file_name);

    int size_b = (M - 2) * (M - 1) / 2;
    mpz_t b_array[size_b];

    bincoef_array(b_array, M - 1);

    mpf_t prob;
    mpf_t sum_prob;

    mpf_init(sum_prob);

    FILE *f;
    f = fopen(file_name, "w");

    fprintf(f, "t,p,F\n");

    for (l = M; l <= n; l++){
        pmf(prob, l, M, pr, pg, b_array); 
        
        mpf_add(sum_prob, sum_prob, prob);
        
        fprintf(f, "%d,", l);

        mpf_out_str(f, 10, 10, prob);

        fprintf(f, ",");

        mpf_out_str(f, 10, 10, sum_prob);

        fprintf(f, "\n");

        if ((l % save_frequency) == 0){
            printf("l: %d, out of [%d, %d]\n", l, M, n);
            fclose(f);
            f = fopen(file_name, "a");

        }

    }

    for (l = 0; l < size_b; l++){
            mpz_clear(b_array[l]);
    }

    mpf_clear(prob);
    mpf_clear(sum_prob);

    fclose(f);


  return 1;
}

