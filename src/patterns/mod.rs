mod base;
pub use base::Pattern;

mod checker;
pub use checker::Checker;

mod gradient;
pub use gradient::Gradient;

mod ring;
pub use ring::Ring;

mod solid;
pub use solid::Solid;

mod stripe;
pub use stripe::Stripe;

#[cfg(test)]
mod test_pattern;
#[cfg(test)]
pub use test_pattern::TestPattern;
