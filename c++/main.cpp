#include<iostream>
#include"ldgm.hpp"

int main(){
  LDGM ldgm(6, 3, 40);
  ldgm.runExperiment(100, 1.00, 0.01);
  return 1;
}
