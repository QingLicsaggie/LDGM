#include<iostream>
#include"ldgm.hpp"

int main(){
  LDGM ldgm(3, 6, 6, 3, 40000);
  ldgm.runExperiment(100, 1.00);
  return 1;
}
