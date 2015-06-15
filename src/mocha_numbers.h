#ifndef _NUMBERS
#define	_NUMBERS

#include <iostream>
#include <math.h>
#include <stdexcept>
#include <string>

namespace mocha_numbers_nmspc
{
  
  using namespace std;
  
  
  class mocha_distance
  {
    
    public:
      double r;
      double pow(int n);
      double invpow(int n);
      
      inline mocha_distance operator+(const mocha_distance &d1)
      {
	mocha_distance out;
	out.r = r + d1.r;
	return out;
      }
      
      inline mocha_distance operator-(const mocha_distance &d1)
      {
	mocha_distance out;
	out.r = r - d1.r;
	return out;
      }
      
      inline mocha_distance operator*(const mocha_distance &d1)
      {
	mocha_distance out;
	out.r = r * d1.r;
	return out;
      }
      
      inline mocha_distance operator/(const mocha_distance &d1)
      {
	mocha_distance out;
	out.r = r / d1.r;
	return out;
      }
      
  };
  
  double mocha_distance::pow(int n)
  {
    
    int i;
    double temp_r = 1.0;
    
    for(i = 0; i < n; i++) temp_r = temp_r * r;
    
    return temp_r;
    
  }
  
  
  double mocha_distance::invpow(int n)
  {
    
    int i;
    double temp_r = 1.0;
    
    for(i = 0; i < n; i++) temp_r = temp_r * r;
    
    return 1.0 / temp_r;
    
  }
  
  
  
  
  
  
  
  
  
}






#endif




















