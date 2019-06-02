This is a porting to DGSOL, in order to make it interact with modern C++ code. This required two steps:
- Replacing the old Makefile structure with CMake
- Instead of re-structuring the old Fortran code to make it type-safe (some elements are pointers), I preferred to make the Fortran invocation in C, and then to invoke the C code from C++

All the rights to the original authors (see the License file).
