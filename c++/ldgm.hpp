#ifndef LDGM_HPP
#define LDGM_HPP
#include<iostream>
#include<vector>
#include<algorithm>
#include<math.h>
#include<cfloat>
#include<unordered_set>
using namespace std;
class LDGM{
 private:
  int dv2;   //variable node degree for LDGM part
  int dc2;   //check    node degree for LDGM part
  
  int c;     //number of nodes
  int e;     //# of total edges
  int m;     //m = dv2*c, # of check nodes for LDGM.
  int n;     //n = dc2*c, # of variable nodes for LDPC


  vector<double>E; //message along each edge

  //Ccon2 is edge connection for check node of LDGM
  //Vcon  is edge connection for variable node
  vector<vector<int> >Ccon2, Vcon;
  
  vector<int> U; //U is the value in LDPC bit node, i.e., variable nodes.
 public:
  //constructor function.
  LDGM(int v2, int c2, int cc);
  
  //Variable-node update 
  void variableUpdate( );

 //Check-node update for LDGM 
  void checkUpdateLDGM(unordered_set<int>DeletedVariable, double beta, vector<int>x);
    
  //compute the maximal marginal probability 
  double Change(vector<double>Previous);

  //Decrimation
  void Decrimation(unordered_set<int>&DeletedVariable, double beta);
  

  //obtain check nodes of LDGM based on variable nodes of LDGM
  vector<int> obtainLDGMCheck(vector<int> variable);

  //running experiment now
  double runExperiment(int experimentNum, double beta, double e);
  
  //calculateRewritingCost
  double rewritingCost(vector<int> x, vector<int> y);

  //print out for debugging
  void printE();

  //compute the rate
  double rate();
	
  //compute the theoretical rewriting cost.
  double theoreticalRewritingCost(double rate);
};
#endif//ldgm.hpp

