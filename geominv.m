% geominv.m - an attempt at a MATLAB-based inversion tool for geodetic data
%
% uses the fmincon function, with an appropriate algorithm to solve for the
% best-fitting inverse model for a bounded set of input parameters
%
% requires matlab optimization toolkit, i think
%
% early days
%
% gjf, 23-jan-2014
% genius in action (tm)
%
% version history
%
% 

% set-up information

clear;

nrestarts = 100;

faultapriorimodel  = [140; 50; -89; 1.0; 770*1e3; 4230*1e3; 10*1e3; 1*1e3; 10*1e3];
parameterbounds    = [ 10; 10;  20; 0.9;  20*1e3;   20*1e3;  5*1e3; 1*1e3;  3*1e3];

indatafile = 'dinar.okinv';

% select appropriate solver for fmincon and stop it from vomiting
% everything to screen

inversionoptions=optimoptions('fmincon','Algorithm','interior-point','Display','off');

% load data in okinv format

disp('geominv v1.00, written by gareth funning')
disp(sprintf(' reading data from input file %s',indatafile));
indatafid = fopen(indatafile);
tmpindata = fscanf(indatafid,'%f %f %f %f %f %f %*s');
ndata = length(tmpindata)/6;
indata = reshape(tmpindata,6,ndata)';
indata(:,1:2)=indata(:,1:2)*1e3;

% set up a model parameter matrix

modelparams = zeros(length(faultapriorimodel),3);
modelparams(:,1)=faultapriorimodel;

% and define which parameters will be inverted

parammapping=(1:9)';
inversionparams=zeros(length(parammapping),nrestarts);

% set up an anonymous function to pass values into the solver

anonyfaultmodel=@(x)faultmodel(x,indata,modelparams,parammapping);

% initialize output arrays

localminima=zeros(length(parammapping),nrestarts);
penalties=zeros(1,nrestarts);

% loop through different starting models

disp(sprintf('  starting fmincon solver, with %d monte carlo restarts',nrestarts));
tic;

parfor i=1:nrestarts


% select a random starting model within the bounds

    inversionparams(:,i) = (faultapriorimodel-parameterbounds) + (2*parameterbounds).*rand(length(faultapriorimodel),1);

% and optimise with fmincon

    [localminima(:,i),penalties(i)]=fmincon(anonyfaultmodel,inversionparams(:,i),[],[],[],[],faultapriorimodel-parameterbounds,faultapriorimodel+parameterbounds,[],inversionoptions);


end

t=toc;

disp(sprintf('   %d restarts completed in %f seconds',nrestarts,t));

% produce output -- find the minimum misfit model
% and calculate forward-computed displacements and fault coordinates

bestpen=min(penalties);
bestmodelno=find(penalties==bestpen);
bestmodel=localminima(:,bestmodelno);

disp(sprintf('    best penalty = %f (rms = %f)',bestpen,sqrt(bestpen/ndata)));

[displacements,surfacetrace,faultframe]=forwardmodel(bestmodel,indata,modelparams,parammapping);

save displacements.dat displacements -ascii;
save surfacetrace.dat surfacetrace -ascii;


%  (optionally) produce statistics on convergence

disp(' ')
disp('     convergence history:')
bestpen=1e30;
bestpens=zeros(1,nrestarts);

for i=1:nrestarts
    
    if (penalties(i)<bestpen)
        bestpen=penalties(i);
        disp(sprintf('      iteration #%d, best penalty now %f',i,bestpen));
    end
    
    bestpens(i)=bestpen;
    
end


plot(bestpens)

% (optionally) display model

disp(' ')
disp('     global miniumum model parameters:')
disp(sprintf('      strike: %f',bestmodel(1)));
disp(sprintf('      dip:    %f',bestmodel(2)));
disp(sprintf('      rake:   %f',bestmodel(3)));
disp(sprintf('      slip:   %f',bestmodel(4)));
disp(sprintf('      x:      %f',bestmodel(5)));
disp(sprintf('      y:      %f',bestmodel(6)));
disp(sprintf('      length: %f',bestmodel(7)));
disp(sprintf('      top:    %f',bestmodel(8)));
disp(sprintf('      bottom: %f',bestmodel(9)));
disp(' ')

