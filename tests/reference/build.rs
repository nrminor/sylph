fn main() {
    println!(
        "cargo:rustc-env=REFERENCE_TARGET={}",
        std::env::var("TARGET").expect("Cargo target")
    );
}
