fn main() {
    let rdkit_root = std::env::var("RDROOT")
        .unwrap_or_else(|_| "/home/brent/omsf/clone/rdkit".to_owned());

    let include = "/home/brent/Projects/rdkit-sys/include";
    println!(
        "cargo:rustc-env=LD_LIBRARY_PATH={include}:{rdkit_root}/build/lib"
    );
}
