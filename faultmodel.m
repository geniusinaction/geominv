function penalty = faultmodel(inversionparams,indata,modelparams,parammapping)

% function penalty = faultmodel(inversionparams,indata,modelparams,parammapping)
%
% function that computes a penalty function based upon a set of input
% parameters 
%
% output: 
%   'penalty'
%       scalar total squared residual [sum((data-model_prediction)^2)]
%
% input: 
%   'inversionparams'
%       vector (1 row for each model or nuisance parameter being inverted)
%
% uses the following global variables 
%   'indata' 
%       matrix (6 columns - xpos, ypos, displ, losx, losy, losz)
%   'modelparams'
%       matrix (9 rows per fault, 4 columns)
%              (strike, dip, rake, slip, x, y, length, top, bottom)
%   'nuisanceparams'
%       matrix (3 rows per dataset, 4 columns)
%

% global variables

%global indata;              % the geodetic data
%global modelparams;         % the fault parameters and elastic parameters
%global parammapping;
%global nuisanceparams;      % gradients and DC shifts

% compute the forward slip model, then...

for i=1:length(inversionparams)

   modelparams(parammapping(i),1)=inversionparams(i);
    
end

% first, parse the input
 
[r,c]=size(modelparams);
nfaults=r/9;
[ndata,c]=size(indata);

obsdispl_los=indata(:,3);

outdisps=zeros(ndata,3);
modeldispl_los=zeros(ndata,1);

% loop through the fault segments

for i=1:nfaults
    
% grab the input parameters

    strike=modelparams(1+(i-1)*9,1);
    dip=modelparams(2+(i-1)*9,1);
    
    xupdip=modelparams(5+(i-1)*9,1);
    yupdip=modelparams(6+(i-1)*9,1);
    
    faultlength=modelparams(7+(i-1)*9,1);

    slip=modelparams(4+(i-1)*9,1);
    rake=modelparams(3+(i-1)*9,1);
    
    faultbottom=modelparams(9+(i-1)*9,1);
    faulttop=modelparams(8+(i-1)*9,1);

% do some annoying trig to account for the okada code's input geometry format
    
    width=(faultbottom-faulttop)/sind(dip);
    centroiddepth=(faultbottom+faulttop)/2;
    gamma=centroiddepth/tand(dip);
    
    xcentroid=xupdip-gamma*sind(strike-90);
    ycentroid=yupdip-gamma*cosd(strike-90);
    
    
    
    [outdisps(:,1),outdisps(:,2),outdisps(:,3)]=okada85(indata(:,1)-xcentroid,indata(:,2)-ycentroid,centroiddepth,strike,dip,faultlength,width,rake,slip,0,0.25);
    


% calculate los deformation

modeldispl_los=modeldispl_los+outdisps(:,1).*indata(:,4)+outdisps(:,2).*indata(:,5)+outdisps(:,3).*indata(:,6);

end

b=obsdispl_los-modeldispl_los;

A=[ones(ndata,1) indata(:,1)-xcentroid indata(:,2)-ycentroid];

x=A\b;

penalty = (b-A*x)'*(b-A*x);

%penalty = b'*b;

