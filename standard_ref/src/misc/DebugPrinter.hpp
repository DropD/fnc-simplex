/*******************************************************************************
 * 
 * DebugPrinter.hpp - build 20130425
 * (C) 2011-2013 Donjan Rodic <donjan@dyx.ch>
 * Released under the WTFPL (www.wtfpl.net) with no warranty whatsoever.
 * 
 * 
 * Creates a global static object named  dout  and defines the macros  dout_HERE
 * and  dout_FUNC .
 * 
 * Link with -rdynamic in order to properly get stack() frame names.
 * 
 * Pass DEBUGPRINTER_NO_EXECINFO flag on Windows (disables stack method).
 * 
 * Pass DEBUGPRINTER_NO_CXXABI if you don't have a cxxabi demangle call in your
 * libc distribution (disables demangle method). stack() will then print raw
 * stack frame names and a c++filt-ready output.
 * 
 * Pass DEBUGPRINTER_OFF to turn off all functionality provided here. You can
 * leave the debug statements in your code, since all methods and macros become
 * inline and empty, and thus should be optimised away by your compiler.
 * 
 * Note: keep in mind that compiler optimisations may inline (-> shorter stack).
 * 
 * 
 * Usage:
 *   dout << "foo" << std::endl;
 *   dout << "foo ", 5, " bar ", 6 << " foobar ", 7;
 * 
 *   dout.stack();
 *   dout.stack(int count);       // print at most count frames
 *   dout.stack(int count, true); // compact format
 *   dout_FUNC                    // shortcut for  dout.stack(1, true);
 * 
 *   dout(anything);              // highlighted format, anything must have a
 *                                //   std::ostream::operator<< overload
 *   dout(label, anything);       // label must have std::ostream::operator<<
 *   dout_HERE                    // shortcut for  dout(__func__, __LINE__);
 * 
 *   dout << dout.demangle( typeid(object).name() );
 *   dout << dout.demangle( custom_string, cxx_demangle_status_int & );
 * 
 * Precision of float output (default == 5) can be set with
 *   dout.precision(12);
 * 
 * The output std::ostream (default == std::cout) can be set with
 *   dout = std::cerr;
 * 
 * The highlighted output bash color (default == "36" == cyan) can be set with
 *   dout.hcolor("31");
 * 
 ******************************************************************************/

#pragma once

#include <iostream>

#ifndef DEBUGPRINTER_OFF

#include <iomanip>
#include <cstdlib>

#ifndef DEBUGPRINTER_NO_EXECINFO
#include <execinfo.h>
#endif // DEBUGPRINTER_NO_EXECINFO

#ifndef DEBUGPRINTER_NO_CXXABI
#include <cxxabi.h>
#endif // DEBUGPRINTER_NO_CXXABI

#endif // DEBUGPRINTER_OFF


#ifndef DEBUGPRINTER_OFF

class DebugPrinter {

  public:

/*******************************************************************************
 * Ctor and friends
 */

  DebugPrinter() : outstream(&std::cout), prec(5), hcol("36") {}

  template <typename T>
  friend DebugPrinter & operator<<(DebugPrinter & d, const T& output);
  friend inline DebugPrinter & operator<<(DebugPrinter & d,
                                          std::ostream& (*pf)(std::ostream&));

/*******************************************************************************
 * Setters
 */

  void operator=(std::ostream & os) { outstream = &os; }
  void precision(int prec_) { prec = prec_; }
  void hcolor(std::string str) { hcol = str; }  // no chaining -> void

/*******************************************************************************
 * Parentheses operators
 */

  template <typename T, typename U>
  inline void operator()(const T& label, const U& line,
                         const std::string sc = ": ") const {
    *outstream << "\033[0;" << hcol << "m" << ">>> " << label << sc << line
               << "\033[0m" << std::endl;
  }
  template <typename T>
  inline void operator()(const T& line) const { operator()("", line, ""); }

/*******************************************************************************
 * Demangling methods
 */

  #ifndef DEBUGPRINTER_NO_CXXABI

  inline std::string demangle(const std::string & str, int & status) const {
    char * demangled = abi::__cxa_demangle(str.c_str(), 0, 0, &status);
    std::string out = (status==0) ? std::string(demangled) : str;
    free(demangled);
    return out;
  }

  #else // DEBUGPRINTER_NO_CXXABI

  inline std::string demangle(const std::string & str, int & status) const {
    return str;
  }

  #endif // DEBUGPRINTER_NO_CXXABI

  inline std::string demangle(const std::string & str) const {
    int dummy;
    return demangle(str, dummy);
  }

/*******************************************************************************
 * Stack methods
 */

  #ifndef DEBUGPRINTER_NO_EXECINFO

  void stack(
             const int backtrace_size = max_backtrace,
             const bool compact = false,
             const int begin = 1 /* should stay 1 except for special needs */
            ) const {

    std::ostream & out = *outstream;
    int _end = -1;               // ignore last

    void * stack[max_backtrace];
    int r = backtrace(stack, backtrace_size+begin-_end);
    char ** symbols = backtrace_symbols(stack, r);
    if(!symbols) return;

    int end = r + _end;
    if(compact == false)
      out << "\nObtained " << end-begin << " stack frames:" << std::endl;

    #ifndef DEBUGPRINTER_NO_CXXABI

    for(int i = begin; i < end; ++i) {
      std::string prog = "  " + prog_part(std::string(symbols[i])) + ":  ";
      std::string mangled = mangled_part(std::string(symbols[i]));
      std::string offset = "  +" + offset_part(std::string(symbols[i]));
      if(compact == true) { prog = ""; offset = ""; }
      int status;
      std::string demangled = demangle(mangled, status);
      switch (status) {
        case -1:
          out << " Error: Could not allocate memory!" << std::endl;
          break;
        case -2:  // invalid name under the C++ ABI mangling rules
          out << prog << mangled << offset << std::endl;
          break;
        case -3:
          out << " Error: Invalid argument to demangle()" << std::endl;
          break;
        default:
          out << prog << demangled << offset << std::endl;
      }
    }
    if(compact == false) out << std::endl;

    #else // DEBUGPRINTER_NO_CXXABI

    for(int i = begin; i < end; ++i) {
      if(compact == false)
        out << "  " << symbols[i] << std::endl;
      else
        out << mangled_part(std::string(symbols[i])) << std::endl;
    }

    if(compact == false) out << std::endl;
    out << "echo '' && c++filt";
    for(int i = begin; i < end; ++i)
      out << " " << mangled_part(std::string(symbols[i]));
    out << " && echo ''" << std::endl;
    if(compact == false) out << std::endl;

    #endif // DEBUGPRINTER_NO_CXXABI

    free(symbols);
  }

  // The dout_FUNC macro pollutes the namespace but is more convenient
  //~ #ifndef DEBUGPRINTER_NO_CXXABI
  //~ void func() { stack(1, true, 2); }
  //~ #endif // DEBUGPRINTER_NO_CXXABI

  #else // DEBUGPRINTER_NO_EXECINFO

  void stack(
             const int backtrace_size = max_backtrace,
             const bool compact = false,
             const int begin = 1 /* should stay 1 except for special needs */
            ) const {
    *outstream << "DebugPrinter::stack() not available" << std::endl;
  }

  #endif // DEBUGPRINTER_NO_EXECINFO

/*******************************************************************************
 * Private parts
 */

  private:

  std::ostream * outstream;
  int prec;
  std::string hcol;

  static const unsigned int max_backtrace = 50;
  static const unsigned int max_demangled = 4096;

  inline std::string prog_part(const std::string str) const {
    return str.substr(0, str.find("("));
  }
  inline std::string mangled_part(const std::string str) const {
    std::string::size_type pos = str.find("(") + 1;
    return str.substr(pos, str.find("+", pos) - pos);
  }
  inline std::string offset_part(const std::string str) const {
    std::string::size_type pos = str.find("+") + 1;
    return str.substr(pos, str.find(")", pos) - pos);
  }

};

/*******************************************************************************
 * std::ostream overloads
 */

template <typename T>
DebugPrinter & operator<<(DebugPrinter & d, const T& output) {
  std::ostream & out = *d.outstream;
  size_t savep = (size_t)out.precision();
  std::ios_base::fmtflags savef =
                 out.setf(std::ios_base::fixed, std::ios::floatfield);
  out << std::setprecision(d.prec) << std::fixed << output
            << std::setprecision(savep);
  out.setf(savef, std::ios::floatfield);
  out.flush();
  return d;
}

inline DebugPrinter & operator<<(DebugPrinter & d,      // manipulators overload
                                 std::ostream& (*pf)(std::ostream&)) {
  std::ostream & out = *d.outstream;
  out << pf;
  return d;
}

template <typename T>
inline DebugPrinter & operator,(DebugPrinter & d, const T& output) {
  return operator<<(d, output);
}

inline DebugPrinter & operator,(DebugPrinter & d,       // manipulators overload
                                std::ostream& (*pf)(std::ostream&)) {
  return operator<<(d, pf);
}

/*******************************************************************************
 * Globals
 */

static DebugPrinter dout;

#define dout_HERE dout(__func__, __LINE__);
#define dout_FUNC dout.stack(1, true);


/**
 * End
 ******************************************************************************/

#else // DEBUGPRINTER_OFF


class DebugPrinter {
  public:

  void operator=(std::ostream & os) {}
  void precision(int prec_) {}
  void hcolor(std::string str) {}
  template <typename T, typename U>
  inline void operator()(const T& label, const U& line,
                         const std::string sc = "") const {}
  template <typename T>
  inline void operator()(const T& line) const {}
  inline std::string demangle(const std::string & str,
                              int & status) const { return ""; }
  inline std::string demangle(const std::string & str) const { return ""; }
  inline void stack(const int backtrace_size = 0, const bool compact = false,
                    const int begin = 0) const {}
};

template <typename T>
inline DebugPrinter & operator<<(DebugPrinter & d, const T& output) {return d;}
inline DebugPrinter & operator<<(DebugPrinter & d,
                                 std::ostream& (*pf)(std::ostream&)) {return d;}
template <typename T>
inline DebugPrinter & operator,(DebugPrinter & d, const T& output) {return d;}
inline DebugPrinter & operator,(DebugPrinter & d,
                                std::ostream& (*out)(std::ostream&)) {return d;}

static DebugPrinter dout;

#define dout_HERE ;
#define dout_FUNC ;


#endif // DEBUGPRINTER_OFF

