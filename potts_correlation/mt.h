//---------------------------------------------------------------------------
//          Copyright H. Watanabe 2015.
// Distributed under the Boost Software License, Version 1.0.
//    (See accompanying file LICENSE_1_0.txt or copy at
//     http://www.boost.org/LICENSE_1_0.txt)
//----------------------------------------------------------------------
#ifndef mt_h
#define mt_h
//---------------------------------------------------------------------------
class MT {
public:
  static void SetSeed(int seed);
  static double GetDouble(void);
  static double GetGauss(void);
  static void Save(char * filename);
  static void Load(char * filename);
};
//---------------------------------------------------------------------------
#endif
