Place the public header files in this directory. They will be
available to your code (and other modules) with

     #include <IMP/micos/myheader.h>

All headers should include `IMP/micos/micos_config.h` as their
first include and surround all code with `IMPMICOS_BEGIN_NAMESPACE`
and `IMPMICOS_END_NAMESPACE` to put it in the IMP::micos namespace
and manage compiler warnings.

Headers should also be exposed to SWIG in the `pyext/swig.i-in` file.
