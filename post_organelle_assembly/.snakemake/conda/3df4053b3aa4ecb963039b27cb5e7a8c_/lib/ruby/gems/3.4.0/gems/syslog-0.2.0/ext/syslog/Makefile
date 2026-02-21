
SHELL = /bin/sh

# V=0 quiet, V=1 verbose.  other values don't work.
V = 0
V0 = $(V:0=)
Q1 = $(V:1=)
Q = $(Q1:0=@)
ECHO1 = $(V:1=@ :)
ECHO = $(ECHO1:0=@ echo)
NULLCMD = :

#### Start of system configuration section. ####
top_srcdir = $(topdir)/.
gem_srcdir = $(top_srcdir)/.bundle/gems/syslog-0.2.0
srcdir = $(gem_srcdir)/ext/syslog
topdir = ../../../../..
hdrdir = $(top_srcdir)/include
arch_hdrdir = $(extout)/include/$(arch)
PATH_SEPARATOR = :
VPATH = $(srcdir):$(arch_hdrdir)/ruby:$(hdrdir)/ruby
RUBYLIB =
RUBYOPT = -
prefix = $(DESTDIR).
rubysitearchprefix = $(rubylibprefix)/$(sitearch)
rubyarchprefix = $(rubylibprefix)/$(arch)
rubylibprefix = $(libdir)/$(RUBY_BASE_NAME)
exec_prefix = $(prefix)
vendorarchhdrdir = $(vendorhdrdir)/$(sitearch)
sitearchhdrdir = $(sitehdrdir)/$(sitearch)
rubyarchhdrdir = $(extout)/include/$(arch)
vendorhdrdir = $(rubyhdrdir)/vendor_ruby
sitehdrdir = $(rubyhdrdir)/site_ruby
rubyhdrdir = $(top_srcdir)/include
vendorarchdir = $(vendorlibdir)/$(sitearch)
vendorlibdir = $(vendordir)/$(ruby_version)
vendordir = $(rubylibprefix)/vendor_ruby
sitearchdir = $(sitelibdir)/$(sitearch)
sitelibdir = $(sitedir)/$(ruby_version)
sitedir = $(rubylibprefix)/site_ruby
rubyarchdir = $(extout)/$(arch)
rubylibdir = $(extout)/common
sitearchincludedir = $(includedir)/$(sitearch)
archincludedir = $(includedir)/$(arch)
sitearchlibdir = $(libdir)/$(sitearch)
archlibdir = $(libdir)/$(arch)
ridir = $(datarootdir)/$(RI_BASE_NAME)
modular_gc_dir = $(DESTDIR)
mandir = $(datarootdir)/man
localedir = $(datarootdir)/locale
libdir = $(DESTDIR)/Users/runner/miniforge3/conda-bld/ruby_1765963609795/work
psdir = $(docdir)
pdfdir = $(docdir)
dvidir = $(docdir)
htmldir = $(docdir)
infodir = $(datarootdir)/info
docdir = $(datarootdir)/doc/$(PACKAGE)
oldincludedir = $(DESTDIR)/usr/include
includedir = $(SDKROOT)$(prefix)/include
runstatedir = $(localstatedir)/run
localstatedir = $(prefix)/var
sharedstatedir = $(prefix)/com
sysconfdir = $(prefix)/etc
datadir = $(datarootdir)
datarootdir = $(prefix)/share
libexecdir = $(exec_prefix)/libexec
sbindir = $(exec_prefix)/sbin
bindir = $(exec_prefix)/bin
archdir = $(rubyarchdir)
buildlibdir = $(DESTDIR)/Users/runner/miniforge3/conda-bld/ruby_1765963609795/work


CC_WRAPPER = 
CC = arm64-apple-darwin20.0.0-clang
CXX = arm64-apple-darwin20.0.0-clang++
LIBRUBY = $(LIBRUBY_SO)
LIBRUBY_A = lib$(RUBY_SO_NAME)-static.a
LIBRUBYARG_SHARED = -l$(RUBY_SO_NAME)
LIBRUBYARG_STATIC = -l$(RUBY_SO_NAME)-static -framework CoreFoundation $(MAINLIBS)
empty =
OUTFLAG = -o $(empty)
COUTFLAG = -o $(empty)
CSRCFLAG = $(empty)

RUBY_EXTCONF_H = extconf.h
cflags   = $(hardenflags) -fdeclspec  $(optflags) $(debugflags) $(warnflags)
cxxflags = 
optflags = -O3 -fno-fast-math
debugflags = -ggdb3
warnflags = -Wall -Wextra -Wextra-tokens -Wdeprecated-declarations -Wdivision-by-zero -Wdiv-by-zero -Wimplicit-function-declaration -Wimplicit-int -Wpointer-arith -Wshorten-64-to-32 -Wwrite-strings -Wold-style-definition -Wmissing-noreturn -Wno-cast-function-type -Wno-constant-logical-operand -Wno-long-long -Wno-missing-field-initializers -Wno-overlength-strings -Wno-parentheses-equality -Wno-self-assign -Wno-tautological-compare -Wno-unused-parameter -Wno-unused-value -Wunused-variable -Wmisleading-indentation -Wundef
cppflags = 
CCDLFLAGS = -fno-common
CFLAGS   = $(CCDLFLAGS) -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/include -fdebug-prefix-map=/Users/runner/miniforge3/conda-bld/ruby_1765963609795/work=/usr/local/src/conda/ruby-3.4.8 -fdebug-prefix-map=/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_=/usr/local/src/conda-prefix $(cflags) -fno-common $(ARCH_FLAG)
INCFLAGS = -I. -I$(arch_hdrdir) -I$(hdrdir) -I$(srcdir)
DEFS     = 
CPPFLAGS = -DRUBY_EXTCONF_H=\"$(RUBY_EXTCONF_H)\" -D_FORTIFY_SOURCE=2 -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/include -mmacosx-version-min=11.0 -mmacosx-version-min=11.0 -D_XOPEN_SOURCE -D_DARWIN_C_SOURCE -D_DARWIN_UNLIMITED_SELECT -D_REENTRANT $(DEFS) $(cppflags)
CXXFLAGS = $(CCDLFLAGS) -ftree-vectorize -fPIC -fstack-protector-strong -O2 -pipe -stdlib=libc++ -fvisibility-inlines-hidden -fmessage-length=0 -isystem /Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/include -fdebug-prefix-map=/Users/runner/miniforge3/conda-bld/ruby_1765963609795/work=/usr/local/src/conda/ruby-3.4.8 -fdebug-prefix-map=/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_=/usr/local/src/conda-prefix -fdeclspec $(ARCH_FLAG)
ldflags  = -L. -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/lib -L/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/lib -fstack-protector-strong
dldflags = -Wl,-headerpad_max_install_names -Wl,-dead_strip_dylibs -Wl,-rpath,/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/lib -L/Users/user/Downloads/04. Pipeline_NGS/pipeline_assembly/post_organelle_assembly/.snakemake/conda/3df4053b3aa4ecb963039b27cb5e7a8c_/lib -Wl,-multiply_defined,suppress -Wl,-undefined,dynamic_lookup 
ARCH_FLAG = -arch arm64
DLDFLAGS = $(ldflags) $(dldflags) $(ARCH_FLAG)
LDSHARED = $(CC) -dynamic -bundle
LDSHAREDXX = $(CXX) -dynamic -bundle
POSTLINK = dsymutil $@ 2>/dev/null; { test -z '$(RUBY_CODESIGN)' || codesign -s '$(RUBY_CODESIGN)' $@; }
AR = arm64-apple-darwin20.0.0-ar
LD = arm64-apple-darwin20.0.0-ld
EXEEXT = 

RUBY_INSTALL_NAME = $(RUBY_BASE_NAME)
RUBY_SO_NAME = ruby.3.4
RUBYW_INSTALL_NAME = 
RUBY_VERSION_NAME = $(RUBY_BASE_NAME)-$(ruby_version)
RUBYW_BASE_NAME = rubyw
RUBY_BASE_NAME = ruby

arch = arm64-darwin20
sitearch = $(arch)
ruby_version = 3.4.0
ruby = /Users/runner/miniforge3/conda-bld/ruby_1765963609795/_build_env/share/rubygems/bin/ruby --disable=gems -I. -rarm64-darwin20-fake  -I'$(topdir)'
RUBY = $(ruby)
BUILTRUBY = $(topdir)/miniruby
ruby_headers = $(hdrdir)/ruby.h $(hdrdir)/ruby/backward.h $(hdrdir)/ruby/ruby.h $(hdrdir)/ruby/defines.h $(hdrdir)/ruby/missing.h $(hdrdir)/ruby/intern.h $(hdrdir)/ruby/st.h $(hdrdir)/ruby/subst.h $(arch_hdrdir)/ruby/config.h $(RUBY_EXTCONF_H)

RM = rm -f
RM_RF = rm -fr
RMDIRS = rmdir -p
MAKEDIRS = mkdir -p
INSTALL = /usr/bin/install -c
INSTALL_PROG = $(INSTALL) -m 0755
INSTALL_DATA = $(INSTALL) -m 644
COPY = cp
TOUCH = exit >

#### End of system configuration section. ####

preload = 
EXTSO = 
libpath = . $(topdir)
LIBPATH = -L. -L$(topdir)
DEFFILE = 

CLEANFILES = mkmf.log
DISTCLEANFILES = 
DISTCLEANDIRS = 

extout = $(topdir)/.ext
extout_prefix = $(extout)$(target_prefix)/
target_prefix = 
LOCAL_LIBS = 
LIBS = $(LIBRUBYARG_SHARED)  -lpthread  
ORIG_SRCS = syslog.c
SRCS = $(ORIG_SRCS) 
OBJS = syslog.o
HDRS = 
LOCAL_HDRS = 
TARGET = syslog_ext
TARGET_NAME = syslog_ext
TARGET_ENTRY = Init_$(TARGET_NAME)
DLLIB = $(TARGET).bundle
EXTSTATIC = 
STATIC_LIB = $(TARGET).a

TARGET_TOPDIR = $(topdir)/.bundle
TIMESTAMP_DIR = $(TARGET_TOPDIR)/.timestamp
BINDIR        = $(extout)/bin
RUBYCOMMONDIR = $(extout)/common
RUBYLIBDIR    = $(RUBYCOMMONDIR)$(target_prefix)
RUBYARCHDIR   = $(extout)/$(arch)$(target_prefix)
HDRDIR        = $(extout)/include/ruby$(target_prefix)
ARCHHDRDIR    = $(extout)/include/$(arch)/ruby$(target_prefix)
TARGET_GEM_DIR = $(TARGET_TOPDIR)/extensions/$(gem_platform)/$(ruby_version)/syslog-0.2.0
TARGET_SO_DIR = $(TARGET_GEM_DIR)$(target_prefix)/
TARGET_SO     = $(TARGET_SO_DIR)$(DLLIB)
CLEANLIBS     = $(TARGET_SO) $(TARGET_SO:=.dSYM)
CLEANOBJS     = $(OBJS) $(OBJS:.o=.bc) $(OBJS:.o=.i) $(OBJS:.o=.s) *.bak
TARGET_SO_DIR_TIMESTAMP = $(TIMESTAMP_DIR)/$(gem_platform)/$(ruby_version)/.gems.-.syslog-0.2.0.time

gem_platform = arm64-darwin-20

# default target
all:

gem = syslog-0.2.0

build_complete = $(TARGET_GEM_DIR)/gem.build_complete
install-so: build_complete
clean-so:: clean-build_complete

build_complete: $(build_complete)
$(build_complete): $(TARGET_SO)
	$(Q) $(TOUCH) $@

clean-build_complete:
	-$(Q)$(RM) $(build_complete)

install: gemspec
clean: clean-gemspec

gemspec = $(TARGET_TOPDIR)/specifications/$(gem).gemspec
$(gemspec): $(gem_srcdir)/.bundled.$(gem).gemspec
	$(Q) $(MAKEDIRS) $(@D)
	$(Q) $(COPY) $(gem_srcdir)/.bundled.$(gem).gemspec $@

gemspec: $(gemspec)

clean-gemspec:
	-$(Q)$(RM) $(gemspec)

install-rb: gemlib
clean-rb:: clean-gemlib

LN_S = ln -s
CP_R = cp -r

gemlib = $(TARGET_TOPDIR)/gems/$(gem)/lib
gemlib:
	$(Q) $(NULLCMD) $(RUBY) $(top_srcdir)/tool/ln_sr.rb -q -f -T $(gem_srcdir)/lib $(gemlib)

clean-gemlib:
	$(Q) $(NULLCMD) $(gemlib)

all:    install
static: all
.PHONY: all install static install-so install-rb
.PHONY: clean clean-so clean-static clean-rb

clean-static::
clean-rb-default::
clean-rb::
clean-so::
clean: clean-so clean-static clean-rb-default clean-rb
		-$(Q)$(RM_RF) $(CLEANLIBS) $(CLEANOBJS) $(CLEANFILES) .*.time

distclean-rb-default::
distclean-rb::
distclean-so::
distclean-static::
distclean: clean distclean-so distclean-static distclean-rb-default distclean-rb
		-$(Q)$(RM) Makefile $(RUBY_EXTCONF_H) conftest.* mkmf.log exts.mk
		-$(Q)$(RM) core ruby$(EXEEXT) *~ $(DISTCLEANFILES)
		-$(Q)$(RMDIRS) $(DISTCLEANDIRS) 2> /dev/null || true

realclean: distclean
install: install-so install-rb

install-so: $(TARGET_SO)
clean-so::
	-$(Q)$(RM) $(TARGET_SO) $(TARGET_SO_DIR_TIMESTAMP)
	-$(Q)$(RM_RF) $(CLEANLIBS)
	-$(Q)$(RMDIRS) $(TARGET_SO_DIR) 2> /dev/null || true
clean-static::
	-$(Q)$(RM) $(STATIC_LIB)
install-rb: pre-install-rb do-install-rb install-rb-default
install-rb-default: pre-install-rb-default do-install-rb-default
pre-install-rb: Makefile
pre-install-rb-default: Makefile
do-install-rb:
do-install-rb-default:
pre-install-rb-default:
	@$(NULLCMD)
$(TARGET_SO_DIR_TIMESTAMP):
	$(Q) $(MAKEDIRS) $(@D) $(TARGET_SO_DIR)
	$(Q) $(TOUCH) $@

site-install: site-install-so site-install-rb
site-install-so: install-so
site-install-rb: install-rb

.SUFFIXES: .c .m .cc .mm .cxx .cpp .o .S

.cc.o:
	$(ECHO) compiling $(<)
	$(Q) $(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -c $(CSRCFLAG)$<

.cc.S:
	$(ECHO) translating $(<)
	$(Q) $(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -S $(CSRCFLAG)$<

.mm.o:
	$(ECHO) compiling $(<)
	$(Q) $(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -c $(CSRCFLAG)$<

.mm.S:
	$(ECHO) translating $(<)
	$(Q) $(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -S $(CSRCFLAG)$<

.cxx.o:
	$(ECHO) compiling $(<)
	$(Q) $(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -c $(CSRCFLAG)$<

.cxx.S:
	$(ECHO) translating $(<)
	$(Q) $(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -S $(CSRCFLAG)$<

.cpp.o:
	$(ECHO) compiling $(<)
	$(Q) $(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -c $(CSRCFLAG)$<

.cpp.S:
	$(ECHO) translating $(<)
	$(Q) $(CXX) $(INCFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(COUTFLAG)$@ -S $(CSRCFLAG)$<

.c.o:
	$(ECHO) compiling $(<)
	$(Q) $(CC) $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) $(COUTFLAG)$@ -c $(CSRCFLAG)$<

.c.S:
	$(ECHO) translating $(<)
	$(Q) $(CC) $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) $(COUTFLAG)$@ -S $(CSRCFLAG)$<

.m.o:
	$(ECHO) compiling $(<)
	$(Q) $(CC) $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) $(COUTFLAG)$@ -c $(CSRCFLAG)$<

.m.S:
	$(ECHO) translating $(<)
	$(Q) $(CC) $(INCFLAGS) $(CPPFLAGS) $(CFLAGS) $(COUTFLAG)$@ -S $(CSRCFLAG)$<

$(TARGET_SO): $(OBJS) Makefile $(TARGET_SO_DIR_TIMESTAMP)
	$(ECHO) linking shared-object $(DLLIB)
	-$(Q)$(RM) $(@)
	$(Q) $(LDSHARED) -o $@ $(OBJS) $(LIBPATH) $(DLDFLAGS) $(LOCAL_LIBS) $(LIBS)
	$(Q) $(POSTLINK)

$(STATIC_LIB): $(OBJS)
	-$(Q)$(RM) $(@)
	$(ECHO) linking static-library $(@)
	$(Q) $(AR) rcu $@ $(OBJS)
	-$(Q)arm64-apple-darwin20.0.0-ranlib $(@) 2> /dev/null || true

$(OBJS): $(HDRS) $(ruby_headers)
