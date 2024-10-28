#!/usr/bin/env bash
#!/usr/bin/env bash

cd workflow/scripts/pmd-mask
cargo clean
git checkout "v0.3.2"

export CARGO_INCREMENTAL=0
export CC="x86_64-conda-linux-gnu-gcc"
export CFLAGS="-I./ -I$CONDA_PREFIX/include -L./ -L $CONDA_PREFIX/lib"
export CPPFLAGS="${CFLAGS}"
export LDFLAGS="-L$CONDA_PREFIX/lib"
export PERL5LIB="$CONDA_PREFIX/lib/perl5"
CMAKE_C_FLAGS="-I./ -L./ -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib" \
RUSTFLAGS="-Ctarget-cpu=native" \
cargo install --all-features --path . --root $CONDA_PREFIX

