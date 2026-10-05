#[cfg(test)]
mod test;

/// some helper code that contains builder logic to call the function because
/// generic parameters and const generics can make the function slightly awkward
/// to call.
mod builder;

pub(crate) mod coefficients;

pub use coefficients::coefficients2d;
pub use coefficients::coefficients2d_blocked;

pub use builder::Coefficients2d;
