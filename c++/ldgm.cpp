#include<iostream>
#include<unordered_set>
#include<cmath>
#include<cstdlib>
#include<time.h>
#include"ldgm.hpp"
/*----------------------------------------------
Constructor for LDGM. No need to say more about it.
  ------------------------------------------------*/
LDGM:: LDGM(int v2, int c2, int cc){ 
    dv2 = v2;
    dc2 = c2;
    c   = cc;
    
    m  = dv2*c;           // # of check nodes for LDGM.
    n  = dc2*c;          // # of variable nodes for LDGM 
    e  = dv2*n;          //# of edges for LDGM

    E.resize(e + 1);
    U.resize(n);
   
    //Vcon=[reshape(1:e,dv2,[])]; % variable-node connections
    for(int i = 0; i <n; i++){
      vector<int> row(dv2, 0);
      Vcon.push_back(row);	
    }

    //Vcon(: , 1: dv2) = reshape(1: e, dv2, [])
    int index = 1;
    for(int i = 0; i < n; i++)
      for(int j = 0; j < dv2; j++)
	Vcon[i][j] = index++;
    
    
    for(int i = 0; i < m; i++){
      vector<int>row(dc2, 0);
      Ccon2.push_back(row);
    }
    
    //Ccon2=reshape(randperm(e1),[],dc2); % check-node connections in LDPC
    vector<int>randomVector(e);
    for(int i = 0; i < e; i++)
      randomVector[i] = i + 1;
    
    srand(time(0));
    random_shuffle(randomVector.begin(), randomVector.end());
    index = 0;
    for(int i = 0; i < m; i++)
      for(int j = 0; j < dc2; j++)
      	Ccon2[i][j] = randomVector[index++];
		
    /*Ccon2[0][0] = 5; Ccon2[0][1] = 11; Ccon2[0][2] = 17;
    Ccon2[1][0] = 4; Ccon2[1][1] = 10; Ccon2[1][2] = 1;
    Ccon2[2][0] = 12; Ccon2[2][1] = 7; Ccon2[2][2] = 2;
    Ccon2[3][0] = 3; Ccon2[3][1] = 18; Ccon2[3][2] = 14;
    Ccon2[4][0] = 9; Ccon2[4][1] = 6; Ccon2[4][2] = 16;
    Ccon2[5][0] = 8; Ccon2[5][1] = 13; Ccon2[5][2] = 15;
    //print out for debugging
    cout<<"Vcon:"<<endl;
    for(int i = 0; i < n; i++){
      for(int j = 0; j < dv2; j++)
    	cout<<Vcon[i][j]<<" ";
      cout<<endl;
    }
    cout<<endl;
       
    cout<<"Ccon2: "<<endl;
    for(int i = 0; i < m; i++){
      for(int k = 0; k < dc2; k++){
    	cout<<Ccon2[i][k]<<" ";
      }
      cout<<endl;
      }*/
  }
/*---------------------------------------------------------
FUNCTIONS:
           variableUpdate()
DESCRIPTION:
           E(Vcon)=repmat(sum(E(Vcon),2),1,dv2)-E(Vcon);
PARAMETERS:
      INPUT:
          --None
      OUTPUT:
          --Update E based on density evolution rule for variable,
	   the updated rule in matlab language is
            E(Vcon)=repmat(sum(E(Vcon),2),1,dv2)-E(Vcon);
RETURN VALUES:
        NONE
  -----------------------------------------------------------*/
void LDGM:: variableUpdate(){
  //updating rule for variable node v_i = u_0 + \sum^{dv - 1} u_j
  vector<double>copyE(E);
  for(int variable = 0; variable < Vcon.size(); variable++){
    for(int degree = 0; degree < Vcon[0].size(); degree++){
	for(int d = 0; d < Vcon[0].size(); d++)if(d != degree){
	  E[Vcon[variable][degree]] += copyE[Vcon[variable][d]];
      }//end of for(d)
    }//end of for(degree...)
  }//end of for(variable...)
}
/*-------------------------------------------------------------
FUNCTION:
         checkUpdateLDGM(unordered_set<int>DeletedVariable, double beta, vector<int>x)
DESCRIPTION:
            Aref's paper equation 20.
            E(Vcon(MDEL,dv2))=Inf; 
            E(Ccon2)=(repmat(((-1).^x)*tanh(beta),1,dc2)).*GenProd(E(Ccon2),beta);
            E(Ccon2)=(1/beta)*atanh(E(Ccon2)); 
PARAMETERS:
        INPUT:
	   DeletedVariable -- deleted variable nodes.
	   beta-- a real number and the optimal number is based on experiment
	   x   -- received vector with value either 0 or 1, the size of which is m.
        OUTPUT:
	   Update E. The update rule in matlab language is
     
	    E(Vcon(MDEL,dv2))=Inf; 
            E(Ccon2)=(repmat(((-1).^x)*tanh(beta),1,dc2)).*GenProd(E(Ccon2),beta);
            E(Ccon2)=atanh(E(Ccon2)),
	    
	    note that the function for GenProd(V, beta) in matlab language is
	    
	    function [Vout]=GenProd(Vin,beta)

	    n=size(Vin,2); m=size(Vin,1);
	    Vin=tanh(beta*Vin);

	    Vf=cumprod(Vin(:,1:(n-1)),2); Vf=[ones(m,1) Vf]; //forward
	    Vb=cumprod(Vin(:,n:-1:2),2); Vb=Vb(:,n-1:-1:1); Vb=[Vb ones(m,1)]; //backward
	    Vout=Vf.*Vb;
  -------------------------------------------------------------*/
void LDGM:: checkUpdateLDGM(unordered_set<int>DeletedVariable, double beta, vector<int>x){
  /*------------------------------
    first step is for deleted edgeds connecting to LDGM check node and set their values as infinity.
  i.e., this step is corresponding to E(Vcon(MDEL,dv2))=Inf;
  --------------------------------*/
  for(auto it = DeletedVariable.begin(); it != DeletedVariable.end(); it++){
    int index = *it;
    for(int i = 0; i < dv2; i++)
	E[Vcon[index][i]] =  DBL_MAX;
  }
  /*-------------------------------------------------
  the next step is to update it based on belief propagation for check node, i.e., this step is
  corresponding to E(Ccon2)=(repmat(((-1).^x)*tanh(beta),1,dc2)).*GenProd(E(Ccon2),beta);
  E(Ccon2)=atanh(E(Ccon2)); 
  -------------------------------------------------------*/
  vector<double>copyE(E);
  for(int a = 0; a < Ccon2.size(); a++){
    for(int i = 0; i < Ccon2[0].size(); i++){
	double temp = 1.0;
	//the following part is to implement GenProd(E(Ccon2),beta)
	for(int degree = 0; degree < Ccon2[0].size(); degree++)
	  if(degree != i)
	    temp *= tanh(beta*copyE[Ccon2[a][degree]]);
	  
	//the following part is to implement -1^x_i*tanh(beta)
	double tempPow = x[a] == 0?  tanh(beta): -1.0*tanh(beta);
	//the following part is to implement 1/beta*atanh((repmat(((-1).^x)*tanh(beta),1,dc2)).*GenProd(E(Ccon2),1))
	E[Ccon2[a][i]] = atanh(tempPow*temp)/beta;
    }//end of for(i..)
  }//end of for(a...)
}
/*----------------------------------------
FUNCTION:
       double change(vector<double>Previous)
DESCRIPTION:
            equation 23 of aref's paper 
       
PARAMETERS:
          INPUT:
	    E-- as in the class
	    Previous--previous e
          OUTPUT:
	     |E|^{-1} \sum{|Previous[i] - E[i]|}
RETURN VALUES
        B
  ------------------------------------------*/
double LDGM:: Change(vector<double> Previous){
  double ret = 0.0;

  for(int i = 1; i < E.size(); i++)
    ret += abs(Previous[i] - E[i]);
  return ret/e;
}
/*--------------------------------------
FUNCTION:
       int Decrimation(unordered_set<int>DeletedVariable, double beta)
DESCRIPTION:
       Decrimation process, and the details can be referred to algorithm 1 
       of Aref's paper.
PARAMETERS:
       INPUT:
          DeletedVariable--deleted variable nodes.
	  beta  ---
       OUTPUT:
          Determine a deleted variable node, add it to the DeletedVariable and storage the value in U
RETURN VALUES:
        return decrimated index
  ----------------------------------------*/
int LDGM::Decrimation(unordered_set<int>&DeletedVariable, double beta){
  vector<double>message(n, 0.0);
  //This part is to implement, bt=sum(E(Vcon),2);
  for(int i = 0; i < n; i++)
      for(int degree = 0; degree < Vcon[0].size(); degree++)
	  message[i] += E[Vcon[i][degree]];

  srand(time(NULL));
  
  //this is to implement  [B,i]=max(abs(bt(ND)));
  int index = 0;
  double tempMax = 0.0;
  for(int i = 0; i < n; i++)
    if(tempMax < abs(message[i]) && DeletedVariable.find(i) == DeletedVariable.end()){
      tempMax = abs(message[i]);
      index = i;
    }

  //cout<<"max is "<<tempMax<<" and index is "<<index<<endl;
  //from the underlying ensemble, we choose a node uniformly. From Aref's paper, i.e., underline eq. 25. 
  if(tempMax == 0.0){
     vector<int>temp;
     for(int i = 0; i < n; i++){
       if(DeletedVariable.find(i) == DeletedVariable.end())
	 temp.push_back(i);
     }
     index = rand()%temp.size();
     if(rand()%2 == 0)
       U[temp[index]] = 0;
     else
       U[temp[index]] = 1;
      //set V_{dec} = V_{dec} U{i}
     //cout<<"choose "<<temp[index]<<" and set to "<<U[temp[index]]<<endl;;
      DeletedVariable.insert(temp[index]);
      return temp[index];
  }
  else{
    /*------------
      set u_{i^*} to 0 or 1 with probability (1 + tanh(beta*m_{i^*}))/2 or (1 - tanh(beta*m_{i^*})/2, i.e., randomized decision from Aref's paper
      -------------*/
    vector<int>Max;
    for(int i = 0; i < n; i++){
      if(abs(message[i]) == tempMax)
	Max.push_back(i);
    }
    index = rand()%Max.size();
    if((double)rand()/(double)RAND_MAX <= (1 + tanh(beta*message[Max[index]]))/2)
      U[Max[index]] = 0;
    else
      U[Max[index]] = 1;
    //set V_{dec} = V_{dec} U{i}
    //cout<<"choose "<<Max[index]<<" and set to "<<U[Max[index]]<<endl;
    DeletedVariable.insert(Max[index]);
    return Max[index];
  }
}
/*------------------------------------------------
FUNCTION:
      void updatey(vector<int>&x, int index, int value)
DESCRIPTION:
      update y based on the newly obtain value of bit index
PARAMETERS:
      INPUT:
          x --
	  index--bit index
	  value -- bit value
      OUTPUT:
          updated LDGM check nodes.
RETURN VALUE:
      A vector of interger indicating the LDGM check nodes.
   --------------------------------------------------*/
void LDGM:: updatey(vector<int>&x, int index,  int value){
  if(value == 0) return;
  
  unordered_set<int>edges;
  for(int i = 0; i < dv2; i++)
    edges.insert(Vcon[index][i]);

  for(int i = 0; i < m; i++){
    for(int d = 0; d < dc2; d++){
      if(edges.find(Ccon2[i][d]) != edges.end())
	x[i] = (x[i] + 1)%2;
    }
  }
}
/*-------------------------------------------------
FUNCTION:
     double rewritingCost(vector<int> x)
DESCRIPTION:
     Obtain the hamming weight of x.
  ---------------------------------------------------*/
double LDGM:: rewritingCost(vector<int> x){
  int res = 0;
  for(int i = 0; i < m; i++)
    if(x[i])
      res++;
  return ((double) res)/m;
}
/*-------------------------------------------------
FUNCTION:
       void printE()
DESCRIPTION:
       Print out E for debugging.
  --------------------------------------------------*/
void LDGM:: printE(){
  //cout<<"-------------------"<<endl;
  cout<<"E is as follows"<<endl;
  for(int i = 1; i < E.size()-1; i++)
    cout<<i<<" ---- "<<E[i]<<" "<<endl;
  cout<<"-------------------"<<endl;
}
/*----------------------------------------------
FUNCTION:
	double rate()
DESCRIPTION:
	Compute the rate of the WEM, which is m/n
PARAMETERS:
	INPUT: 
	   m, n as the class variables
        OUTPUT:
	   rate
RETURN VALUE
	rate
-----------------------------------------------*/
 double LDGM:: rate(){
   return (double)n/m;
}
/*-----------------------------------------------
FUNCTION:
       double theoreticalRewritingCost(double rate)
DESCRIPTION:
       Given the WEM rate, compute the theoretical rewriting cost,  which is H^-(rate).
-------------------------------------------------*/
double LDGM::theoreticalRewritingCost(double rate){
   for(double res = 0.0001; res < 0.5; res += 0.0001){
	double entropy = -res*log2(res) - (1-res)*log2(1-res);
	//cout<<res<<"---"<<entropy<<"----"<<rate<<endl;
	if(abs(entropy - rate) < 0.001){
	  //  cout<<"cost is "<<res<<endl;
	return res;
	}
   }
}
/*---------------------------------------------------
FUNCTION:
       double runExperiment(int experimentNum, double beta, double e)
DESCRIPTION:
       Given the number of iterations we have to run, run 
       the experiment and the return value is the average rewriting cost.
PARAMETERS:
       INPUT:
           beta--inverse temperature.
	   ee   -- as the variable in equation 23
  ----------------------------------------------------*/
double LDGM:: runExperiment(int experimentNum, double beta, double ee){
  double res = 0.0, Rate = rate(), cost = theoreticalRewritingCost(1 - Rate);


  //for each iteration of the experiment
  for(int iteration = 0; iteration < experimentNum;   iteration++){
    cout<<"experiment "<<iteration<<endl;

    //generate a beroulli symmetric source word x
    vector<int>x(m);
    for(int i = 0; i < m; i++)
      x[i] = rand()%2;

    /*x[0] = 0; x[1] = 1; x[2] = 1; 
      x[3] = 0; x[4] = 0; x[5] = 1;*/

    fill_n(U.begin(), n, 0);
    fill_n(E.begin(), e + 1, 0.0);
    
    /*cout<<"x is "<<endl;
    for(int i = 0; i < m; i++)
      cout<<x[i]<<" ";
      cout<<endl;*/
    unordered_set<int>DeletedVariable;
    
    for(int count = 0; count < n; count ++){
      if(count % 100 == 0)
        cout<<"count now is "<<count<<" and n is "<<n<<endl;
	vector<double>previous(E);
      //for each bit, run a couple of times or equation 23 does not hold. 10 is suggested by Aref.
      for(int i = 0; i < 10; i++){
	//cout<<"iteration "<<i<<endl;
	variableUpdate();
	//cout<<"variable"<<endl;
	//printE();

	checkUpdateLDGM(DeletedVariable, beta, x);
	//cout<<"check"<<endl;
	//printE();
	
	/*E(Vcon(MDEL,:))=0;*/
        for(auto it = DeletedVariable.begin(); it != DeletedVariable.end();  it ++){
	  int index = *it;
	  for(int degree = 0; degree < Vcon[0].size(); degree ++)
	    E[Vcon[index][degree]] = 0.0;
	}
        //if(Change(previous) >= ee) break;
      }

      int bitIndex = Decrimation(DeletedVariable, beta);
      updatey(x, bitIndex, U[bitIndex]);
      
      /*cout<<"now x is "<<endl;
      for(int i = 0; i < m; i++)
	cout<<x[i]<<" ";
	cout<<endl;*/
    }//end for for(...count)
    res += rewritingCost(x);
    if(iteration%1==0) cout<<" So far average rewriting cost is "<<res/(1+iteration)<<",rate is "<<Rate<<" and theoretical value is "<<cost<<endl;
  }//end for for(... experimentNum ...)
  return res/experimentNum;
}
