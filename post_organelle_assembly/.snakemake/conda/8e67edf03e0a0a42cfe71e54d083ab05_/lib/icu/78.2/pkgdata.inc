GENCCODE_ASSEMBLY_TYPE=-a gcc-darwin
SO=dylib
SOBJ=dylib
A=a
LIBPREFIX=lib
LIB_EXT_ORDER=.78.2.dylib
COMPILE=x86_64-apple-darwin13.4.0-clang -DU_ALL_IMPLEMENTATION -DU_ATTRIBUTE_DEPRECATED= -D_FORTIFY_SOURCE=2 -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_/include -mmacosx-version-min=10.13 -mmacosx-version-min=10.13    -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_/include -fdebug-prefix-map=/Users/runner/miniforge3/conda-bld/icu_1767969687066/work=/usr/local/src/conda/icu-78.2 -fdebug-prefix-map=/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -fno-common -c
LIBFLAGS=-I/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_/include  -dynamic
GENLIB=x86_64-apple-darwin13.4.0-clang -dynamiclib -dynamic -march=core2 -mtune=haswell -mssse3 -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_/include -fdebug-prefix-map=/Users/runner/miniforge3/conda-bld/icu_1767969687066/work=/usr/local/src/conda/icu-78.2 -fdebug-prefix-map=/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_/lib -L/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_/lib -Wl,-rpath,/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_/lib 
LDICUDTFLAGS=
LD_SONAME=-Wl,-compatibility_version -Wl,78 -Wl,-current_version -Wl,78.2 -install_name /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/8e67edf03e0a0a42cfe71e54d083ab05_/lib/
RPATH_FLAGS=
BIR_LDFLAGS=
AR=x86_64-apple-darwin13.4.0-ar
ARFLAGS=r -c
RANLIB=x86_64-apple-darwin13.4.0-ranlib
INSTALL_CMD=/usr/bin/install -c
