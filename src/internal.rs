use crate::{Error, fma::FusedMuladd};
use num_traits::{ConstOne, ConstZero, Float};

#[cfg(test)]
mod test;

/// compute N^(i,epsilon') in eq (51) in the B&M paper, where i is the index
/// of the corresponding pole and epsilon' is calculated from the given epsilon
/// as described in section 5.1, p 116/117
pub fn compute_truncation_indices<F, const NPOLES: usize>(
    poles: &[F; NPOLES],
    epsilon: F,
) -> Result<[usize; NPOLES], Error>
where
    F: Float + ConstOne + ConstZero + FusedMuladd,
{
    if epsilon <= F::ZERO || !epsilon.is_normal() {
        return Err(Error::Precision);
    }

    debug_assert!(
        poles
            .iter()
            .all(|pole| { *pole < F::ZERO && *pole > -F::ONE && pole.is_finite() })
    );

    let rhon = compute_rhon(poles);
    let mu = compute_mus(poles);
    // NOTE that in section 5.1 in the B&M paper, the authors say that they
    // use the truncation index from eqn (51) but with a modified \epsilon'
    // where \epsilon' = rho^(n) \epsilon/2.
    let epsilon_prime = epsilon * rhon / (F::ONE + F::ONE);
    // NOTE: we calculate ln(epsilon' * rho^(n)) here
    let ln_epsilon_prime_rhon = (epsilon_prime * rhon).ln();

    let mut n_trunc = [0; NPOLES];
    let mut mu_prod = F::ONE;

    for idx in (0..NPOLES).rev() {
        let pole = poles[idx];
        let ftrunc_idx = (ln_epsilon_prime_rhon
            + F::ln((F::ONE - pole) * (F::ONE - mu[idx]) * mu_prod))
            / F::ln(-pole);
        n_trunc[idx] = 1 + ftrunc_idx.floor().to_usize().ok_or(Error::Limits)?;
        mu_prod = mu_prod * mu[idx];
    }
    Ok(n_trunc)
}

/// compute `rho^n` from eq (52) from Briand and Monasse
// NOTE: I wish we had const fn in traits, because the compute_rhon and compute_mu
// functions could be comptime in theory...
#[inline(always)]
fn compute_rhon<F>(poles: &[F]) -> F
where
    F: Float + ConstOne + ConstZero,
{
    debug_assert!(
        poles
            .iter()
            .all(|pole| { *pole < F::ZERO && *pole > -F::ONE && pole.is_finite() }),
        "poles must be in (-1,0)"
    );

    let prod = poles.iter().copied().fold(F::ONE, |prod, zj| {
        let factor = (F::ONE + zj) / (F::ONE - zj);
        prod * factor
    });
    prod.powi(2)
}

/// compute the {mu_j}_{j \in [0,NPOLES)} according to eq (50) from the B&M paper
#[inline(always)]
fn compute_mus<F, const NPOLES: usize>(poles: &[F; NPOLES]) -> [F; NPOLES]
where
    F: Float + ConstOne + ConstZero + FusedMuladd,
{
    let mut mu = [F::ZERO; NPOLES];

    let mut sum = F::ZERO;
    let mut temp_pole = poles[0].mul(-F::ONE).ln();
    for idx in 1..poles.len() {
        sum = sum + temp_pole.powi(-1);
        temp_pole = poles[idx].mul(-F::ONE).ln();
        mu[idx] = temp_pole * sum / (temp_pole.fmuladd(sum, F::ONE));
    }

    mu
}
