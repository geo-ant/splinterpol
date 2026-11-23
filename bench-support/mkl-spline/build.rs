#[cfg(not(feature = "mkl"))]
fn main() {}

#[cfg(feature = "mkl")]
fn main() {
    println!("cargo:rerun-if-env-changed=MKLROOT");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_PATH");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_LIBDIR");
    println!("cargo:rerun-if-env-changed=PKG_CONFIG_SYSROOT_DIR");
    println!("cargo:rerun-if-changed=c/spline.c");

    let config = "mkl-static-lp64-seq".parse().unwrap();
    let library = match intel_mkl_tool::Library::pkg_config(config) {
        Ok(Some(library)) => library,
        _ => {
            let mkl_root = std::env::var_os("MKLROOT")
                .expect("MKL not found via pkg-config; set MKLROOT to your MKL installation");
            let mkl_root = std::path::PathBuf::from(mkl_root)
                .canonicalize()
                .expect("MKLROOT must point to an existing MKL installation");
            intel_mkl_tool::Library::seek_directory(config, &mkl_root)
                .expect("failed to inspect MKLROOT")
                .expect("MKL headers and libraries not found under MKLROOT")
        }
    };

    cc::Build::new()
        .file("c/spline.c")
        .include(&library.include_dir)
        .compile("mkl_spline_shim");

    bindgen::Builder::default()
        .header("c/spline.c")
        .clang_arg(format!("-I{}", library.include_dir.display()))
        .allowlist_type("MklInterpolator_1d_f64")
        .allowlist_function("MklInterpolator_1d_f64_.*")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks::new()))
        .generate()
        .expect("failed to generate MKL spline bindings")
        .write_to_file(
            std::path::PathBuf::from(std::env::var_os("OUT_DIR").unwrap()).join("bindings.rs"),
        )
        .expect("failed to write MKL spline bindings");

    library
        .print_cargo_metadata()
        .expect("failed to configure MKL linking");
}
