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
  int dv1;   //variable node degree for LDPC part
  int dc1;   //check    node degree for LDPC part
  int dv2;   //variable node degree for LDGM part
  int dc2;   //check    node degree for LDGM part
  
  int c;     //number of nodes.
  int m1;    //n*(dv1/dc1); % # of check-nodes in LDPC
  int e1;    //dv1*n; % # of edges in LDPC partm 
  int e2;    //dv2*n # of edges for LDGM
  int e;     //# of total edges
  int m;     //m = dv2*c, # of check nodes for LDGM.
  int n;     //n = dc2*c, # of variable nodes for LDPC


  vector<double>E; //message along each edge
  
  //Ccon1 is edge connection for check node of LDPC
  //Ccon2 is edge connection for check node of LDGM
  //Vcon  is edge connection for variable node
  vector<vector<int> >Ccon1, Ccon2, Vcon;
  
  vector<int> U; //U is the value in LDPC bit node, i.e., variable nodes.
 public:
  //constructor function.
  LDGM(int v1, int c1, int v2, int c2, int cc);
  
  //Variable-node update 
  void variableUpdate( );

 //Check-node update for LDGM 
  void checkUpdateLDGM(unordered_set<int>DeletedVariable, double beta, vector<int>x);
  
  //Check-node update for LDPC
  void checkUpdateLDPC(unordered_set<int>DeletedVariable, double beta);
  
  //compute the maximal marginal probability 
  double computeMaxMarginal(double beta);

  //Decrimation
  void Decrimation(unordered_set<int>&DeletedVariable, double beta);
  
  //parity check error
  bool ParityError(unordered_set<int>&DeletedVariable);

  //check to see whether LDPC checks are satisfied or not
  bool LDPCCheck();

  //obtain check nodes of LDGM based on variable nodes of LDGM
  vector<int> obtainLDGMCheck(vector<int> variable);

  //running experiment now
  double runExperiment(int experimentNum, double beta);
  
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

