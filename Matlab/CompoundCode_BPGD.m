clear; 
% clc;

% For the LDGM Part
dv2=6; % Variable-node degree
dc2=3; % Check-node degree
%%%%%%%%%%%

rate=(dc2/dv2); 
dist_sha=finverse(@(x)(BinaryEntropy(x)),1-rate,0,0.5);
fprintf('rate=%f,Shannon distortion is %f\n',rate,dist_sha);

c=40000;
m=dv2*(c); % Block Length / # of check-nodes in LDGM
n=dc2*(c); % # of variable-nodes
e=dv2*n; % # of edges in LDGM part


% Extra edge to handle irregularity in LDGM part
E=zeros(1,e); % This is the main object in the sim
E_check_code=zeros(1,e); % Object to check if the code-bits satisfy LDPC checks

Vcon= reshape(1:e, dv2, [])'; % variable-node connections

Ccon2=reshape(randperm(e),[],dc2);

M=10; % # of source sequences to average the distortion

beta=1; % inverse temperature
T=10; % Number of BP iterations before a forced decimation

Distortion=zeros(1,M);
tic;
for j=1:M
    UD=zeros(n,1);
    % Bernoulli symmetric source
    y=zeros(m,1); y(rand(m,1)<0.5)=1; % current instance of source seq.
    yr=y;
    E=zeros(1,e);
    count_n=n;
    MDEL=[]; % Variable bits deleted
    while(count_n > 0)
        for i=1:T
            % Variabel-node update
	    E(Vcon)=repmat(sum(E(Vcon),2),1,dv2)-E(Vcon);

            % Check-node LDGM update
            E(Vcon(MDEL, :))=Inf; 
            E(Ccon2)=(repmat(((-1).^yr)*tanh(beta),1,dc2)).*GenProd(E(Ccon2),beta);
            E(Ccon2)=(1/beta)*atanh(E(Ccon2));

            E(Vcon(MDEL,:))=0;
            
            % Possibility of decimation?
            bt=sum(E(Vcon),2);
            B=max(abs(bt));
        end
        if(B > 0)
            I=find(abs(bt)==B);
            nd=I(randi(length(I),1,1));
            ud=rand>(1+tanh(bt(nd)))/2;
	else
	    ND=1:n; ND(MDEL)=[];
            nd=ND(randi(length(ND),1,1));
            ud=(rand<0.5);
        end
        UD(nd)=ud;
        MDEL=[MDEL nd];
        count_n=count_n-1;
        yr=Updatey(yr,Vcon,Ccon2,nd,ud);
        if(mod(count_n,100)==0)
            fprintf('Finished %d/%d code-bits\n',n-count_n,n);
        end
    end
    Distortion(j)=sum(yr)/m;
    fprintf('Distortion for the curr. seq is %f, avg is %f \n',Distortion(j), mean(Distortion(1:j)) );
    fprintf('Shannon Distortion is %f\n',dist_sha);
end
toc;
