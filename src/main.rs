use std::env;

use gmp::mpz::{mpz_struct, Mpz};

const SAVE_FREQUENCY: i32 = 5;
const PG: f32 = 0.95;

fn validate_arguments(args: &[String]) -> (i64, i64) {
    let m: i64 = args[1].parse().expect("M to be a number");
    let n: i64 = args[2].parse().expect("N to be a number");

    assert!(n >= 0 && n >= m);
    (m, n)
}

fn bin_coef(b: &mut Vec<Mpz>, m: i64) {
    for i in 0..m {
        for j in 0..=i {
            // b.push()
        }
    }
}

fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() <= 2 {
        println!("Usage: {} <M> <N>", args[0]);
        return;
    }

    let (m, n) = validate_arguments(&args);
    let pr = PG / n as f32;

    let file_name = format!("../results/pmf_cdf_{}_{}_{}.csv", m, n, (PG * 1000.0) as i32);
    println!("File name: {}", file_name);

    let size_b = (m - 2) * (m - 1) / 2;
    let mut b: Vec<Mpz> = vec!();

    bin_coef(&mut b, m);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    #[should_panic]
    fn test_empty_args() {
        let args = vec!();
        validate_arguments(&args);
    }

    #[test]
    #[should_panic]
    fn test_m_n_non_numbers() {
        let args = vec!("rust".to_owned(), "rocks".to_owned());
        validate_arguments(&args);
    }

    #[test]
    #[should_panic]
    fn test_negative_n() {
        let n = "-1".to_owned();
        let m = "-2".to_owned();
        let args = vec!(m, n);

        validate_arguments(&args);
    }

    #[test]
    #[should_panic]
    fn test_bigger_m() {
        let n = "0".to_owned();
        let m = "1".to_owned();
        let args = vec!(m, n);

        validate_arguments(&args);
    }
}

