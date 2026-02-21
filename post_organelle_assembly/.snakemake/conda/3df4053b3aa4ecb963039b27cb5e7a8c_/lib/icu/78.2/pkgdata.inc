GENCCODE_ASSEMBLY_TYPE=-a gcc-darwin
SO=dylib
SOBJ=dylib
A=a
LIBPREFIX=lib
LIB_EXT_ORDER=.78.2.dylib
COMPILE=arm64-apple-darwin20.0.0-clang -DU_ALL_IMPLEMENTATION -DU_ATTRIBUTE_DEPRECATED= -D_FORTIFY_SOURCE=2 -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/include -mmacosx-version-min=11.0 -mmacosx-version-min=11.0    -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/include -fdebug-prefix-map=/Users/runner/miniforge3/conda-bld/icu_1767969576337/work=/usr/local/src/conda/icu-78.2 -fdebug-prefix-map=/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -fno-common -c
LIBFLAGS=-I/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/include  -dynamic
GENLIB=arm64-apple-darwin20.0.0-clang -dynamiclib -dynamic -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/include -fdebug-prefix-map=/Users/runner/miniforge3/conda-bld/icu_1767969576337/work=/usr/local/src/conda/icu-78.2 -fdebug-prefix-map=/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_=/usr/local/src/conda-prefix -std=c11 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/lib -L/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/lib -Wl,-rpath,/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/lib 
LDICUDTFLAGS=
LD_SONAME=-Wl,-compatibility_version -Wl,78 -Wl,-current_version -Wl,78.2 -install_name /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/lib/
RPATH_FLAGS=
BIR_LDFLAGS=
AR=arm64-apple-darwin20.0.0-ar
ARFLAGS=r -c
RANLIB=arm64-apple-darwin20.0.0-ranlib
INSTALL_CMD=/usr/bin/install -c
