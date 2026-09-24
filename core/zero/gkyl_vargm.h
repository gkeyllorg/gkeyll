#pragma once

/**
   These horrible looking set of macros allows choosing a macro based on
   number of arguments passed to it. It is ugly as hell and so I am
   putting it in its own file. Please do not use it unless you know what
   you are doing. See:

   https://stackoverflow.com/questions/11761703/overloading-macro-on-number-of-arguments/11763277#11763277
*/

// get number of arguments with __NARG__
#define __NARG__(...)  __NARG_I_(__VA_ARGS__,__RSEQ_N())
#define __NARG_I_(...) __ARG_N(__VA_ARGS__)
#define __ARG_N(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N,...) N
#define __RSEQ_N() 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0

#define _VFUNC_(name, n) name##n
#define _VFUNC(name, n) _VFUNC_(name, n)
// general definition for any function name
#define VFUNC(func, ...) _VFUNC(func, __NARG__(__VA_ARGS__)) (__VA_ARGS__)
// general definition for any function name (extra args)
#define VFUNC1(func, a, ...) _VFUNC(func, __NARG__(__VA_ARGS__)) (a, __VA_ARGS__)
