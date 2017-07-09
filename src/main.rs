pub mod dual;

use dual::Dual;

fn main() {
    let x = Dual::new(std::f64::consts::FRAC_PI_2, 1f64);
    println!("{:?}", x.sin());
}
