%macro feedbackOut(L= ,K=, typemtxfut= , PredNamesTD=, DistributionTD=, mydatasorted=, outVar=, idVar=, timeVar= ,
optim=NLPCG, ContUpdate=TRUE);

data typemtxfut;
set &typemtxfut.;
run;

data predNamesTD;
set &PredNamesTD.;
run;

data DistributionTD;
set &DistributionTD.;
run;

data mydatasorted;
set &mydatasorted.;
run;

data L;
set &L.;
run;


*Obtain some summary stats for later use;
data _null_;
	set mydatasorted;
	timepts=max(&timeVar.);
	nsub=_N_ / timepts;
	call symput('Nsubs',nsub);
	call symput('Tpoints',timepts);
run;


%if &K.=1 %then %do;

proc iml;
use typemtxfut;
read all var _ALL_ into typeMTXFUT_M;
close typemtxfut;

use PredNamesTD;
read all var _ALL_ into PredNamesTD;
close PredNamesTD;

use L;
read all var _all_ into L;
close L;

USE mydatasorted;
read all var {&timeVar.} into time2;
read all var {&idVar.} into ID2;
CLOSE mydatasorted;

T=max(time2);
N=ncol(unique(ID2)); *Number of subjects;

%let m=1;

%let j=L[&m.];

/*use following line to create names of time dependent covariates*/
TDPredasOut=predNamesTD[&J.];/*j*/ /*output this dataset*/

/*macro to create names of time-dependent covariates as outcomes*/
start CreateMacroTD(values, macroName, delimiter=' ');
   if type(values)="N" then          
      y = rowvec( char(values) );   /* convert numeric to character */
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);       /* delimit and concatenate */
   s = substr(s, 1, nleng(s)-nleng(delimiter)); /* remove delimiter at end */
   call symputx(macroName, s);      /* create macro variable */
finish;

call CreateMacroTD(TDPredasOut, "PredVarTDasOut");


use DistributionTD;
read all var _ALL_ into DistributionTD;
close DistributionTD;

typematrix=shape(typeMtxFut_m[&J.+1,],T,T); /*j+1, T, T*/
TransposeTypeMTX=t(typematrix);
typeMtxN=shape(transposetypemtx, 1,T*T);

/* define helper functions ROW and COL */
start row(x);  /* return matrix m such that m[i,j] = i */
   return( repeat( T(1:x), 1, x ));
finish;
start col(x);  /* return matrix m such that m[i,j] = j */
   return( repeat(1:x, x) );
finish;

*Helper matrices;
rt = row(T);
ct = col(T);

TypeMtxComboTD = j(T,T*T,0); 

TypeMtxComboTD2 = j(T,T*T,0);
*Loop across the different types, starting with current;
DO i = 0 TO (T-1);
	*Identify the indices for the appropriate setting;
	Idx = loc(rt-ct = i);
	IdxDif = setdif(1:T*T, Idx);
		temp = typeMtxN;
	
	temp[,IdxDif] = 0;
	
	TypeMtxComboTD[(i*1+1):((i+1)*1),] = temp;

	*Shift the values to accomodate shifted Y-X relationship;
	if(i>0) then DO; 
		dummyZero = j(1, i, 0) || temp;
		TypeMtxComboTD2[(i*1+1):((i+1)*1),] = dummyZero[,1:(T*T)];
	END;
	
	ELSE TypeMtxComboTD2[(i*1+1):((i+1)*1),] = temp;

END; 


neqTD = TypeMtxComboTD2[,+]; /*ADDING UP NUMBER OF 1'S FOR EACH ROW OF TYPEMTXCOMBO2*/

*Print a note about omitted moment conditions;
if min(neqTD)=0 then do;
	zeroPredTD = loc(neqTD=0);
	Note = "There are no valid moment conditions for " +char(ncol(zeroPredTD))+" outcome on &PredVarTDasOut covariate relationship(s).";
	print Note[label="Moment Condition Notes"];
	print "These outcome on covariate relationships will be omitted in the analysis.";
end;
else print "All outcome on &PredVarTDasOut covariate relationships will be evaluated."[label="Moment Condition Notes"];


/*FOR THOSE ROWS RELATED TO LAGS OF THE COVARIATE WITH NO VALID MOMENTS WE DELETE THAT ROW*/
/*THIS BECAUSE WE WILL NOT BE USING THOSE LAGS IN OUR MODELS*/
keepPredTD = loc(neqTD>0);


TypeMtxComboTD3 = TypeMtxComboTD2[keepPredTD,]; /*WILL ONLY KEEP ROWS OF COVARIATES WITH VALID MOMENTS*/


*Append the intercept conditions to the matrix;
*Append the time independent conditions to the matrix, if available;
	TypeMtxComboTD3 = shape(I(T),1) // TypeMtxComboTD3; *Type I intercept;


*If there are time independent variables, they are represented after the intercept, before the time dependent;
CREATE TypeMtxComboTD3 from TypeMtxComboTD3;
APPEND from TypeMtxComboTD3;
CLOSE TypeMtxComboTD3;
print "Each row of TypeMtx is the shifted type vector for each of the predictors, by individual test";
print TypeMtxComboTD3;


USE mydatasorted;
read all var {&PredVarTDasOut.} into xnewTD; /*will delete this in loop inside macro*/
read all var {&outVar.} into Ynew;
read all var {&idVar. &timeVar.} into othersnew;
CLOSE mydatasorted;

/*creating lagged outcomes matrix*/

Y2New=j(N*T, T-1, 0);
do i=1 to T-1;
 Y2new[(i*N+1):N*T, (i-1)+1]=Ynew[1:(T-i)*N];
end;

Y2new=Ynew || Y2new;

*Remove the outcomes as predictors with no valid moments;
Y3new=Y2new[,keepPredTD];

*Add in the ID, outcome and time (and possible time independent variables;
TDcovar=othersnew||XnewTD||Y3new; /*dataset that will be used in second macro*/


*Create the adjusted variable names;
outNames = t(repeat(t({&outVar.}),T));
lagName = j(1,T,0);
DO i=1 to (T-1);
	lagName[1,(i*1+1):(i+1)*1]=j(1,1,i);/*nothing changes here because we only have one predictor outcome for each covariate X*/
END;
varNamesTD = catx("_", outNames, char(lagName));


*Retain only the variable names for those with valid moment conditions;
varNamesTD = varNamesTD[,keepPredTD];

/* use following lines to create names of distributions for time dependent covariates*/
distPredVarTD=DistributionTD[&J.];

/*macro to create names of distribution for time-dependent covariates as outcomes*/
start CreateMacroTD(values, macroName, delimiter=' ');
   if type(values)="N" then          
      y = rowvec( char(values) );   /* convert numeric to character */
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);       /* delimit and concatenate */
   s = substr(s, 1, nleng(s)-nleng(delimiter)); /* remove delimiter at end */
   call symputx(macroName, s);      /* create macro variable */
finish;

call CreateMacroTD(distPredVarTD, "distPredVarTD");

*Add in the names for the ID, outcome and time;
varnamesTD2={&idVar.} ||  {&timeVar.} || {&PredVarTDasOut.}|| varNamesTD; 

CREATE TDmodeldata from TDcovar[c=varNamesTD2];
APPEND from TDcovar;
close TDmodeldata;

/*sort the data*/
call sort(TDcovar, {1 2});

/*create data set to output*/
CREATE MydataTD3 from TDcovar[c=varNamesTD2];
APPEND from TDcovar;
close MydataTD3;

/*running GMM model for time dependent covariates*/
start CreateMacro(values, macroName, delimiter=' ');
   if type(values)="N" then          
      y = rowvec( char(values) );   /* convert numeric to character */
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);       /* delimit and concatenate */
   s = substr(s, 1, nleng(s)-nleng(delimiter)); /* remove delimiter at end */
   call symputx(macroName, s);      /* create macro variable */
finish;

call CreateMacro(varNamesTD, "adjPredVarTD");
quit;


*Obtain initial GEE estimates for optimization;
ods select none;
proc genmod data=mydataTD3 descend;
	class &idVar. &timeVar.;
	model &PredVarTDasOut. = &adjPredVarTD. / dist=&distPredVarTD.;
	repeated subject=&idVar. / within=&timeVar. corr=ind corrw;
	OUTPUT OUT=GEEout XBETA=xb RESRAW = rraw;
	ods output GEEEmpPEst=betaGEETD;
RUN;
ods output close;
ods select all;


*Obtain GMM Estimates;
proc iml;
USE mydataTD3;
READ all VAR {&adjPredVarTD.} INTO ZmatTD;
READ all VAR {&PredVarTDasOut.} INTO yvecTD;
READ all VAR {&idVar.} INTO ID;
READ all VAR {&timeVar.} INTO time;
CLOSE mydataTD3;

N=ncol(unique(ID)); *Number of subjects;
T = max(time); *Number of time points;
NpoutTD = ncol(ZmatTD); *Number of predictors, not including intercept;
PnTD=NpoutTD+1;                      * number of covariates/parameters TO estimate, including intercept;

*Create the X matrix;
int = j(N*T,1,1);
XmatTD =j(N*T,PnTD,.); XmatTD[,1]=int; XmatTD[,2:PnTD]=ZmatTD; *XmatTD is matrix of intercept in column 1 (vector of 1) and ZmatTD;

*Use GEE estimates as starting values;
USE betaGEETD;
READ all VAR {Estimate} INTO beta0TD;
CLOSE betaGEETD;
beta0TD = t(beta0TD);

/*Use the identified covariate types from V40 code (multiple test);*/
USE TypeMtxComboTD3;
READ all INTO TypeMtxTD;
CLOSE TypeMtxComboTD3;

*Number of valid moment conditions;
neqTD = TypeMtxTD[,+];

nlocTD = j(1,PnTD+1,0);          * nloc containing the starting/END positions of reg eqs brought by each covariate  ;
DO p =1 TO PnTD;
  nlocTD[p+1] = sum(neqTD[1:p]); 
END;

nregTD = sum(neqTD); * total number of reg_eq, 26;

WnTD = I(nregTD);                * I(26) initial weight matrix ;

*Module to calculate the objective function;
START TSGMM(beta) global(PnTD,T,N,XmatTD,yvecTD,nregTD,TypeMtxTD,nlocTD,WnTD);      
GnTD = j(nregTD,1,0);             * TO initialize Gn whenever objfun is called, b/c Gn depends on beta ;   
eqTD = j(nregTD,N,0);              * 26x1625 reg_eqs for a type I covariate, for all persons  ;
STD = j(nregTD,nregTD,0); 

DO i = 1 TO N;               * DO loops for each person;
  xTD = XmatTD[(i-1)*T+1:i*T,];  * T x Pn, READ data for each person across all times;
  yTD = yvecTD[(i-1)*T+1:i*T];   * T outcomes for person i;

if "&distPredVarTD." = "BIN" then 
  	mu = exp(xTD*t(beta)) / ( 1+exp(xTD*t(beta)) );
  else if "&distPredVarTD." = "NORMAL" then
 	mu = xTD*t(beta);

  Rsd = yTD - mu;                * T x 1;

  DO p = 1 TO PnTD; * Do loops for each predictor;
  	if "&distPredVarTD." = "BIN" then
    	D = xTD[,p]#mu#(1- mu);
	else if "&distPredVarTD." = "NORMAL" then
		D= xTD[,p]#(mu##(-1));
    Eqmtx = Rsd*t(D);	*Gives x#mu#(1-mu)#(y-mu); *Do the time points not matter?;
    eqTD[nlocTD[p]+1:nlocTD[p+1],i] = Eqmtx[loc(TypeMtxTD[p,]^=0)];    *Save calculations to eq;
  END; *End loop across predictors;
  STD = STD + eqTD[,i]*t(eqTD[,i]);
END;     * TO CLOSE "DO i = 1 TO N";
WnTD = ginv(STD/N);		* Update the weight matrix;
GnTD = eqTD[,:];               * row mean => col vector  neq x 1;
f = t(GnTD)*WnTD*GnTD;           * the objective fn TO be minimized; 
RETURN(f);
FINISH TSGMM;


tc = {2000 2000}; optn = {0 1};              * optn[1]=0 specifies a minimization problem;
*Optimization takes: response code, row vector of optimal pts, module, starting value, optimization details, threshold conditions;
  *Run with the GEE estimates as starting parameters;
if "&optim."="NLPNRA" then
	CALL NLPNRA(rc, xres,"TSGMM", beta0TD,optn, , tc);           * DO nonlinear optimization Newton Raphson;
else if "&optim."="NLPNRR" then
	CALL NLPNRR(rc, xres,"TSGMM", beta0TD,optn, , tc);           * DO nonlinear optimization Newton Raphson ridge;
else if "&optim."="NLPQN" then
	CALL NLPQN(rc, xres,"TSGMM", beta0TD,optn, , tc);           * DO nonlinear optimization dual quasi-newton;
else if "&optim."="NLPCG" then
	CALL NLPCG(rc, xres,"TSGMM", beta0TD,optn, , tc);           * DO nonlinear optimization conjugate gradient;
else if "&optim."="NLPQUA" then
	CALL NLPQUA(rc, xres,"TSGMM", beta0TD,optn, , tc);           * DO nonlinear optimization quadratic;
else if "&optim."="NLPTR" then
	CALL NLPTR(rc, xres,"TSGMM", beta0TD,optn, , tc);           * DO nonlinear optimization trust region;
else if "&optim."="NLPNMS" then
	CALL NLPNMS(rc, xres,"TSGMM", beta0TD,optn, , tc);           * DO nonlinear optimization nelder-mead simplex;
else do;
	print "Must use valid optimization, will run with NLPCG";
	CALL NLPCG(rc, xres,"TSGMM", beta0TD,optn, , tc);           * DO nonlinear optimization conjugate gradient;
end;
betaTD = xres;


* Module to calculate the asymptotic variance;
DG = j(nregTD,PnTD,.);    * nreg x Pn;
DO k = 1 TO PnTD; 
  DGi = j(nregTD,N,0);          
  DO i = 1 TO N;               * DO loops for each for each person;
    x = XmatTD[(i-1)*T+1:i*T,];  * T x Pn, READ data for each person across all times;
    y = yvecTD[(i-1)*T+1:i*T];  
	if "&distPredVarTD." = "BIN" then 
    	mu = exp(x*t(betaTD)) / ( 1+exp(x*t(betaTD)) );
	else if "&distPredVarTD." = "NORMAL" then
		mu = x*t(betaTD);
    Rsd = y - mu;                * T x 1;
	if "&distPredVarTD." = "BIN" then do;
    	Dk =  x[,k]#mu#(1- mu);
		Dkz =  x[,k]#(1- 2*mu);
		end;
	else if "&distPredVarTD." = "NORMAL" then do;
		Dk =  x[,k]#(mu##(-1));
		Dkz =  x[,k]#(-mu##(-2));
		end;
    DO p = 1 TO PnTD;
	  if "&distPredVarTD." = "BIN" then 
      	Dp = x[,p]#mu#(1- mu);
	  else if "&distPredVarTD." = "NORMAL" then
		Dp = x[,p]#(mu##(-1));
      Dkzp = Dkz#Dp;                  * T x 1;
	  DGmtx = Rsd*t(Dkzp)-Dk*t(Dp);
      DGi[nlocTD[p]+1:nlocTD[p+1],i] = DGmtx[loc(TypeMtxTD[p,]^=0)];   
    END;
  END;
  DG[,k]= DGi[,:];    * row mean => col vector  neq x 1;
END;   


AsymVar = (1/N)*ginv(t(DG)*WnTD*DG);   * Pn x Pn  note Wn is the inverse ;
AVvec = vecdiag(AsymVar); *Take only the diagonal elements for variance;
StdDev = sqrt(AVvec);

*Calculate the p-Values;
zvalue = t(betaTD)/StdDev;
pvalue = 2*(1-cdf('normal',abs(zvalue)));

*Create a matrix of results;
OutmtxTD = j(PnTD,4,.);
ColLabel={"Estimate" "StdDev" "Zvalue" "Pvalue"};
Varnames=t({"Intercept"} || {&adjPredVarTD.});
OutTitle = "Analysis Feedback from outcome to &PredVarTDasOut. covariate";

OutmtxTD[,1]=t(betaTD);
OutmtxTD[,2]=StdDev;
OutmtxTD[,3]=zvalue;
OutmtxTD[,4]=pvalue;
PRINT OutmtxTD[colname=colLabel rowname=Varnames label=OutTitle ];

*Calculate residuals with GMM estimates;
resvec = yvecTD - XmatTD*t(betaTD);

*Save the results as a dataset;
CREATE outdat1TD FROM OutmtxTD;
APPEND FROM OutmtxTD;
CLOSE outdat1TD;
QUIT;
%end;

%else %do;
proc iml;
use typemtxfut;
read all var _ALL_ into typeMTXFUT_M;
close typemtxfut;

use PredNamesTD;
read all var _ALL_ into PredNamesTD;
close PredNamesTD;

use L;
read all var _all_ into L;
close L;

USE mydatasorted;
read all var {&timeVar.} into time2;
read all var {&idVar.} into ID2;
CLOSE mydatasorted;

T=max(time2);
N=ncol(unique(ID2)); *Number of subjects;

%do m=1 %to &K.;

%let j=L[&m.];

/*use following line to create names of time dependent covariates*/
TDPredasOut=predNamesTD[&j.];/*j*/ /*output this dataset*/
print TDPredasOut;

/*macro to create names of time-dependent covariates as outcomes*/

start CreateMacroTD(values, macroName, delimiter=' ');
   if type(values)="N" then          
      y = rowvec( char(values) );   
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);      
   s = substr(s, 1, nleng(s)-nleng(delimiter)); 
   call symputx(macroName, s);      
finish;

call CreateMacroTD(TDPredasOut, "PredVarTDasOut");


use DistributionTD;
read all var _ALL_ into DistributionTD;
close DistributionTD;

typematrix=shape(typeMtxFut_m[&j.+1,],T,T); /*j+1, T, T*/

TransposeTypeMTX=t(typematrix);
typeMtxN=shape(transposetypemtx, 1,T*T);

/* define helper functions ROW and COL */

start row(x);  /* return matrix m such that m[i,j] = i */ 
   return( repeat( T(1:x), 1, x ));
finish;
start col(x);  /* return matrix m such that m[i,j] = j */ 
   return( repeat(1:x, x) );
finish;

*Helper matrices;
rt = row(T);
ct = col(T);

TypeMtxComboTD = j(T,T*T,0); 

TypeMtxComboTD2 = j(T,T*T,0);
*Loop across the different types, starting with current;
DO i = 0 TO (T-1);
	*Identify the indices for the appropriate setting;
	Idx = loc(rt-ct = i);
	IdxDif = setdif(1:T*T, Idx);
		temp = typeMtxN;
	
	temp[,IdxDif] = 0;
	
	TypeMtxComboTD[(i*1+1):((i+1)*1),] = temp;

	*Shift the values to accomodate shifted Y-X relationship;
	if(i>0) then DO; 
		dummyZero = j(1, i, 0) || temp;
		TypeMtxComboTD2[(i*1+1):((i+1)*1),] = dummyZero[,1:(T*T)];
	END;
	
	ELSE TypeMtxComboTD2[(i*1+1):((i+1)*1),] = temp;

END; 


neqTD = TypeMtxComboTD2[,+]; /*ADDING UP NUMBER OF 1'S FOR EACH ROW OF TYPEMTXCOMBO2*/

*Print a note about omitted moment conditions;

if min(neqTD)=0 then do;
	zeroPredTD = loc(neqTD=0);
	Note = "There are no valid moment conditions for " +char(ncol(zeroPredTD))+" outcome on &PredVarTDasOut covariate relationship(s).";
	print Note[label="Moment Condition Notes"];
	print "These outcome on covariate relationships will be omitted in the analysis.";
end;
else print "All outcome on &PredVarTDasOut covariate relationships will be evaluated."[label="Moment Condition Notes"];


/*FOR THOSE ROWS RELATED TO LAGS OF THE COVARIATE WITH NO VALID MOMENTS WE DELETE THAT ROW*/
/*THIS BECAUSE WE WILL NOT BE USING THOSE LAGS IN OUR MODELS*/

keepPredTD = loc(neqTD>0);


TypeMtxComboTD3 = TypeMtxComboTD2[keepPredTD,]; /*WILL ONLY KEEP ROWS OF COVARIATES WITH VALID MOMENTS*/


*Append the intercept conditions to the matrix;
*Append the time independent conditions to the matrix, if available;

TypeMtxComboTD3 = shape(I(T),1) // TypeMtxComboTD3; *Type I intercept;


*If there are time independent variables, they are represented after the intercept, before the time dependent;
CREATE TypeMtxComboTD3&m. from TypeMtxComboTD3;
APPEND from TypeMtxComboTD3;
CLOSE TypeMtxComboTD3&m.;
print "Each row of TypeMtx is the shifted type vector for each of the predictors, by individual test";


USE mydatasorted;
read all var {&PredVarTDasOut.} into xnewTD; /*will delete this in loop inside macro*/

read all var {&outVar.} into Ynew;
read all var {&idVar. &timeVar.} into othersnew;
CLOSE mydatasorted;

/*creating lagged outcomes matrix*/

Y2New=j(N*T, T-1, 0);
do i=1 to T-1;
 Y2new[(i*N+1):N*T, (i-1)+1]=Ynew[1:(T-i)*N];
end;

Y2new=Ynew || Y2new;

*Remove the outcomes as predictors with no valid moments;
Y3new=Y2new[,keepPredTD];

*Add in the ID, outcome and time (and possible time independent variables;
TDcovar=othersnew||XnewTD||Y3new; /*dataset that will be used in second macro*/


*Create the adjusted variable names;

outNames = t(repeat(t({&outVar.}),T));
lagName = j(1,T,0);
DO i=1 to (T-1);
	lagName[1,(i*1+1):(i+1)*1]=j(1,1,i);
END;
varNamesTD = catx("_", outNames, char(lagName));


*Retain only the variable names for those with valid moment conditions;
varNamesTD = varNamesTD[,keepPredTD];
/* use following lines to create names of distributions for time dependent covariates*/


/*distPredVarTD=DistributionTD[&j.];*/


/*macro to create names of distribution for time-dependent covariates as outcomes*/

/*
call CreateMacroTD(distPredVarTD, "distPredVarTD");
*/
*Add in the names for the ID, outcome and time;
varnamesTD2={&idVar.} ||  {&timeVar.} || {&PredVarTDasOut.}|| varNamesTD; 

CREATE TDmodeldata&m. from TDcovar[c=varNamesTD2];
APPEND from TDcovar;
close TDmodeldata&m.;

/*sort the data*/

call sort(TDcovar, {1 2});

/*create data set to output*/

CREATE MydataTD3&m. from TDcovar[c=varNamesTD2];
APPEND from TDcovar;
close MydataTD3&m.;

/*running GMM model for time dependent covariates*/

start CreateMacro(values, macroName, delimiter=' ');
   if type(values)="N" then          
      y = rowvec( char(values) );   
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);       
s = substr(s, 1, nleng(s)-nleng(delimiter));
   call symputx(macroName, s);      
finish;

/*
call CreateMacro(varNamesTD, "adjPredVarTD");
*/



num=putn(&m., "BEST.6");
name="adjPredVarTD"+strip(num);
call CreateMacro(varNamesTD,name);

%let d=2*&K.;
num2=putn(&K.+&m., "BEST.6");

nameout="PredVarTDasOut"+strip(num);
call CreateMacro(TDPredasOut,nameout);



distPredVarTD=DistributionTD[&j.];
print distPredVarTD;


namedistr="distPredVarTD"+strip(num);
call CreateMacro(distPredVarTD, namedistr);

namedistr2="distPredVarTD"+strip(num2);
call CreateMacro(distPredVarTD, namedistr);


CREATE distPredVarTDData&m. from distPredVarTD[c="distPredVarTD"];
APPEND from distPredVarTD;
close distPredVarTDData&m.;

%end
quit;


%do m=1 %to &K.;
*Obtain initial GEE estimates for optimization;
ods select none;
proc genmod data=mydataTD3&m. descend;
	class &idVar. &timeVar.;
	model &&PredVarTDasOut&m. = &&adjPredVarTD&m. / dist=&&distPredVarTD&m.;
	repeated subject=&idVar. / within=&timeVar. corr=ind corrw;
	OUTPUT OUT=GEEout XBETA=xb RESRAW = rraw;
	ods output GEEEmpPEst=betaGEETD&m.;
RUN;
ods output close;
ods select all;
%end;

*Obtain GMM Estimates;

PROC IML;
N=&Nsubs.; *Number of subjects;
T = &Tpoints.; *Number of time points;


*Module to obtain basic information for each model;
start baseCalcs(TypeMtx, Zmat, neq, nloc, nreg, Pn, Xmat ) global(N, T);
	neq = TypeMtx[,+];
	Pn = ncol(Zmat) + 1; *number of covariates/parameters to estimate, including intercept;
	
	nloc = j(1,Pn+1,0);          * nloc containing the starting/END positions of reg eqs brought by each covariate  ;
	DO p =1 TO Pn;
	  nloc[p+1] = sum(neq[1:p]); 
	END;
	
	nreg = sum(neq); * total number of reg_eq,;

		*Create the X matrix;
	int = j(N*T,1,1);
	Xmat =j(N*T,Pn,.); 
	Xmat[,1]=int; 
	Xmat[,2:Pn]=Zmat; *Xmat is matrix of intercept in column 1 (vector of 1) and Zmat;


finish baseCalcs;

*Initialize storage for the do loop;
beta0 = {0};
Pn = 0;
Pnvec = {0};
nreg = 0;
Xmat = j(N*T,1,0);
yvecs = j(N*T,&K.,0);
TypeMtx = j(1,T*T,0);
modelRep = {0};  *Vector indicating which indices pertain to each model;

*Iterate through the models for moment checks and data restructuring;
%do m=1 %to &K.;

*Assumes that the ID and times are the same for both, so will only read in one of each;
USE mydataTD3&m.;
READ all VAR {&&adjPredVarTD&m.} INTO Zmat&m.;
READ all VAR {&&predVarTDasOut&m.} INTO yvec&m.;
READ all VAR {&idVar.} INTO ID;
READ all VAR {&timeVar.} INTO time;
CLOSE mydataTD3&m.;

use distPredVarTDdata&m.;
read all var _all_ into distPredVarTDdata&m.;
close distPredVarTDdata&m.; 

*Use GEE estimates as starting values;

USE betaGEETD&M.;
READ all VAR {Estimate} INTO beta0&M.;
CLOSE betaGEETD&M.;
beta0&M. = t(beta0&M.);

*Column bind the two beta row vectors;
beta0 = beta0 || beta0&M.;
print beta0;

*Use the identified covariate types;
USE TypeMtxComboTD3&M.;
READ all INTO TypeMtx&M.;
CLOSE TypeMtxComboTD3&M.;

use DistributionTD;
read all var _ALL_ into DistributionTD;
close DistributionTD;


*Obtain basic informatnoi;
call baseCalcs(TypeMtx&M., Zmat&M., neq&M., nloc&M., nreg&M., Pn&M., Xmat&M.);


Pn = Pn + Pn&M.;

nreg = nreg + nreg&M.;

TypeMtx = TypeMtx // TypeMtx&M.;


Xmat = Xmat || Xmat&M.;
yvecs[,&M.] = yvec&M.;

Pnvec = Pnvec // Pn&M.;

modelRep = modelRep || j(1,Pn&M., &M.);

%end; *End loop through models;


*Remove the placeholders;
beta0 = remove(beta0,1);
TypeMtx = TypeMtx[2:Pn+1,];
Xmat = Xmat[,2:Pn+1];
Pnvec = remove(Pnvec,1);
modelRep = remove(modelRep,1);

Pnloc = j(1,&m.+1,0);          * Pnloc containing the starting/END positions of predictors for different models ;
nloc = {0};					* nloc containing the starting/END positions of reg eqs brought by each covariate  ;
shift = 0; 		*Amount to shift nloc;
%DO M =1 %TO &K.;
	Pnloc[&M.+1] = sum(Pnvec[1:&M.]); 
	nloc = nloc || nloc&M.[,2:(Pn&M.+1)]+shift;
	shift = shift+nloc&M.[,Pn&M.+1]; *Increment the location shift for nloc;
%END;

*Define initial values;
S12 = I(nreg);
Wn12 = I(nreg);

nModels = &K.;

*Apply either continuously updating or two-stage GMM;
%if &ContUpdate.=TRUE %then
 	%do;

*Define the module for optimization using continuously updating GMM;
START TSGMM(beta) global(Nmodels, Pnvec, Pnloc,Pn, modelRep, T, N, Xmat, yvecs, nreg, nloc, TypeMtx, Wn12, S12);
*Need an if statement to check whether it is first or second sample;
Gn12 = j(nreg,1,0);
eq = j(nreg,N,0);

DO i=1 to N;
	x = Xmat[(i-1)*T+1:i*T,];  * T x Pn, READ data for each predictor;
	mu = j(T,nModels,0); *Matrix to store mu;
	Rsd = j(T,nModels,0); *Matrix to store residuals;
	Rsdrep = j(T,Pn,0); *Matrix of repeated residuals based on predictors per model;
	D = j(T,Pn,0);
	*Calculate residuals;
	%DO m=1 %to &K.;
	    betaM = beta[,Pnloc[&m.]+1:Pnloc[&m.+1] ];
		xm = x[, Pnloc[&m.]+1:Pnloc[&m.+1] ];
        
		if "&&distPredVarTD&m." = "BIN" then
        mu[,&m.] = exp(xm*t(betaM)) / ( 1+exp(xm*t(betaM)) ); *The mean;
		else if "&&distPredVarTD&m." = "NORMAL" then
 	    mu[,&m.] = xm*t(betaM);


		Rsd[,&m.] = yvecs[(i-1)*T+1:i*T,&m.] - mu[,&m.]; *Residuals, Tx1;
		
		if "&&distPredVarTD&m." = "BIN" then
		D[,Pnloc[&m.]+1:Pnloc[&m.+1]] = xm#mu[,&m.]#(1-mu[,&m.]);
        else if "&&distPredVarTD&m." = "NORMAL" then
		D[,Pnloc[&m.]+1:Pnloc[&m.+1]]=xm#(mu[,&m.]##(-1));
      
	%END; *End loop across models;

	DO p=1 to Pn; *Loop across predictors;
		Eqmtx = Rsd[,modelRep[p]]*t(D[,p]);
		eq[nloc[p]+1:nloc[p+1],i] = Eqmtx[loc(TypeMtx[p,]^=0)];    *Save calculations to eq;
	END; *End loop across predictors;

END;

Gn12 = eq[,:];				* row mean => col vector  neq x 1;
S12 = eq*t(eq);		*Inner products;
Wn12 = ginv(S12/N); *Weight matrix is the inverse;

Q = t(Gn12)*Wn12*Gn12;		*Objective function;

RETURN(Q);
FINISH TSGMM;


tc = {3000 3000 }; optn = {0 2};              * optn[1]=0 specifies a minimization problem;
*NLRPNA takes: response code, row vector of optimal pts, module, starting value, optimization details, threshold conditions;
  *Run with the GEE estimates as starting parameters;
if "&optim."="NLPNRA" then
	CALL NLPNRA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson;
else if "&optim."="NLPNRR" then
	CALL NLPNRR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson ridge;
else if "&optim."="NLPQN" then
	CALL NLPQN(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization dual quasi-newton;
else if "&optim."="NLPCG" then
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
else if "&optim."="NLPQUA" then
	CALL NLPQUA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization quadratic;
else if "&optim."="NLPTR" then
	CALL NLPTR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization trust region;
else if "&optim."="NLPNMS" then
	CALL NLPNMS(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization nelder-mead simplex;
else do;
	print "Must use valid optimization, will run with NLPCG";
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
end;
beta = xres;
%end;

%else %do;

*Define the module for optimization using continuously updating GMM;
START TSGMM(beta) global(Nmodels, Pnvec, Pnloc,Pn, modelRep, T, N, Xmat, yvecs, nreg, nloc, TypeMtx, Wn12, S12);
*Need an if statement to check whether it is first or second sample;
Gn12 = j(nreg,1,0);
eq = j(nreg,N,0);

DO i=1 to N;
	x = Xmat[(i-1)*T+1:i*T,];  * T x Pn, READ data for each predictor;
	mu = j(T,nModels,0); *Matrix to store mu;
	Rsd = j(T,nModels,0); *Matrix to store residuals;
	Rsdrep = j(T,Pn,0); *Matrix of repeated residuals based on predictors per model;
	D = j(T,Pn,0);
	*Calculate residuals;
	*Calculate residuals;
	%DO m=1 %to &K.;
	    betaM = beta[,Pnloc[&m.]+1:Pnloc[&m.+1] ];
		xm = x[, Pnloc[&m.]+1:Pnloc[&m.+1] ];
               
		if "&&distPredVarTD&m." = "BIN" then
        mu[,&m.] = exp(xm*t(betaM)) / ( 1+exp(xm*t(betaM)) ); *The mean;
		else if "&&distPredVarTD&m." = "NORMAL" then
 	    mu[,&m.] = xm*t(betaM);


		Rsd[,&m.] = yvecs[(i-1)*T+1:i*T,&m.] - mu[,&m.]; *Residuals, Tx1;
		
		if "&&distPredVarTD&m." = "BIN" then
		D[,Pnloc[&m.]+1:Pnloc[&m.+1]] = xm#mu[,&m.]#(1-mu[,&m.]);
        else if "&&distPredVarTD&m." = "NORMAL" then
		D[,Pnloc[&m.]+1:Pnloc[&m.+1]]=xm#(mu[,&m.]##(-1));
      
	%END; *End loop across models;


	DO p=1 to Pn; *Loop across predictors;
		Eqmtx = Rsd[,modelRep[p]]*t(D[,p]);
		eq[nloc[p]+1:nloc[p+1],i] = Eqmtx[loc(TypeMtx[p,]^=0)];    *Save calculations to eq;
	END; *End loop across predictors;


END;

Gn12 = eq[,:];				* row mean => col vector  neq x 1;
S12 = eq*t(eq);		*Inner products;

Q = t(Gn12)*Wn12*Gn12;		*Objective function;

RETURN(Q);
FINISH TSGMM;

tc = {3000 3000 }; optn = {0 2};              * optn[1]=0 specifies a minimization problem;
*NLRPNA takes: response code, row vector of optimal pts, module, starting value, optimization details, threshold conditions;
  *Run with the GEE estimates as starting parameters;
if "&optim."="NLPCG" then DO;
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
	beta0 = xres; 
  	Wn12 = ginv(S12/N);
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
END;
else if "&optim."="NLPNRA" then DO;
	CALL NLPNRA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson;
	beta0 = xres; 
  	Wn12 = ginv(S12/N);  
	CALL NLPNRA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson;
END;
else if "&optim."="NLPNRR" then DO;
	CALL NLPNRR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson ridge;
	beta0 = xres; 
  	Wn12 = ginv(S12/N);
	CALL NLPNRR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson ridge;
END:
else if "&optim."="NLPQN" then DO;
	CALL NLPQN(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization dual quasi-newton;
	beta0 = xres; 
  	Wn12 = ginv(S12/N);
	CALL NLPQN(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization dual quasi-newton;	
END;
else if "&optim."="NLPQUA" then DO;
	CALL NLPQUA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization quadratic;
	beta0 = xres; 
  	Wn12 = ginv(S12/N);
	CALL NLPQUA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization quadratic;
END;
else if "&optim."="NLPTR" then DO;
	CALL NLPTR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization trust region;
	beta0 = xres; 
  	Wn12 = ginv(S12/N);
	CALL NLPTR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization trust region;
END;
else if "&optim."="NLPNMS" then DO;
	CALL NLPNMS(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization nelder-mead simplex;
	beta0 = xres; 
  	Wn12 = ginv(S12/N);
	CALL NLPNMS(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization nelder-mead simplex;
END;
else DO;
	print "Must use valid optimization, will run with NLPCG";
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
	beta0 = xres; 
  	Wn12 = ginv(S12/N);
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
END;

beta = xres;
Wn12 = ginv(S12/N); 

%end; *End check for continuously updating or two-stage GMM;



*Module for asymptotic variance;
/*Asymptotic variance should be done by covariate instead of all together*/
*Separate the betas based on model;
start=1;
end=0;
%DO M =1 %TO &K.;
	end = end + Pn&M.;
	beta&m. = beta[,start:end];
	start = start + Pn&M.;


DG&m. = j(nreg&m.,Pn&m.,.);    * nreg x Pn;
DO k = 1 TO Pn&m.; 
  DGi = j(nreg&m.,N,0);              * reg_eqs for a type I covariate, for all persons  ;
  DO i = 1 TO N;               * DO loops for each for each person;
    x = Xmat&m.[(i-1)*T+1:i*T,];  * T x Pn, READ data for each person across all times;
    y = yvec&m.[(i-1)*T+1:i*T];


if "&&distPredVarTD&m." = "BIN" then 
   	mu = exp(x*t(beta&m.)) / ( 1+exp(x*t(beta&m.)));
	else if "&&distPredVarTD&m." = "NORMAL" then
 	mu = x*t(beta&m.);

    Rsd = y - mu;                * T x 1;

if "&&distPredVarTD&m." = "BIN" then do;
    	Dk =  x[,k]#mu#(1- mu);
		Dkz =  x[,k]#(1- 2*mu);
		end;
	else if "&&distPredVarTD&m." = "NORMAL" then do;
		Dk =  x[,k]#(mu##(-1));
		Dkz =  x[,k]#(-mu##(-2));
		end;
   
    DO p = 1 TO Pn&m.;
    if "&&distPredVarTD&m." = "BIN" then 
      	Dp = x[,p]#mu#(1- mu);
	  else if "&&distPredVarTD&m." = "NORMAL" then
		Dp = x[,p]#(mu##(-1));

	  Dkzp = Dkz#Dp;                  * T x 1;
	  DGmtx = Rsd*t(Dkzp)-Dk*t(Dp);
      DGi[nloc&m.[p]+1:nloc&m.[p+1],i] = DGmtx[loc(TypeMtx&m.[p,]^=0)];   
    END;
  END;
  DG&m.[,k]= DGi[,:];    * row mean => col vector  neq x 1;
END;   

%END;

*Combine the derivatives into one matrix;
%if &K.=1 %then %do;
	DG = (DG1 || j(nreg1,Pn1,0)); 
	Varnames=t({"Intercept1"} || {&adjPredVarTD1.});
%end;

%else %if &K.=2 %then %do;
	DG = (DG1 || j(nreg1,Pn2,0)) // (j(nreg2,Pn1,0) || DG2); 
	Varnames=t({"Intercept1"} || {&adjPredVarTD1.} || {"Intercept2"} || {&adjPredVarTD2.});
%end;
%else %if &K.=3 %then %do;
	DG = (DG1 || j(nreg1,Pn2+Pn3,0)) // (j(nreg2,Pn1,0) || DG2 || j(nreg2,Pn3,0)) // (j(nreg3,Pn1+Pn2,0) || DG3 );
	Varnames=t({"Intercept1"} || {&adjPredVarTD1.} || {"Intercept2"} || {&adjPredVarTD2.} || {"Intercept3"} || {&adjPredVarTD3.});
%end;


AsymVar = (1/N)*ginv(t(DG)*Wn12*DG);   * Pn x Pn  note Wn12 is the inverse ;
AVvec = vecdiag(AsymVar); *Take only the diagonal elements for variance;
StdDev = sqrt(AVvec);


*Calculate the p-Values;
zvalue = t(beta)/StdDev;
pvalue = 2*(1-cdf('normal',abs(zvalue)));


*Create a matrix of results;
Outmtx = j(Pn,5,.);
ColLabel={"Estimate" "StdDev" "Zvalue" "Pvalue" "Model"};

OutTitle = "Analysis of Joint GMM Estimates";


Outmtx[,1]=t(beta);
Outmtx[,2]=StdDev;
Outmtx[,3]=zvalue;
Outmtx[,4]=pvalue;
Outmtx[,5]=t(modelRep);
PRINT Outmtx[colname=colLabel rowname=Varnames label=OutTitle ];

*Save the results as a dataset;
CREATE outdat1 FROM Outmtx;
APPEND FROM Outmtx;
CLOSE outdat1;

quit;
%end;
%mend feedbackOut;







/***************************************************************************************/
/**	SAS Code: %partitionedGMMPart1								                      **/
/** Originally Programmed by Kyle Irimata Modified by Elsa Vazquez for Feedback Model **/
/** Description: Performs partitioned GMM regression and provides Matrix with         **/
/** valid moment conditions for feedback of outcome to time dependent covariates      **/
/***************************************************************************************/
%macro partitionedGMMPart1(ds=., file=, timeVar=, outVar=, predVarTD=,distrPredVarTD=, idVar=, alpha=0.05, predVarTI=., distr=bin, optim=NLPCG, MC=LWY);

*Check if there is a library;
%if &ds. ^=. %then %do;
	LIBNAME DS &ds.;  
	DATA mydata; 
		SET DS.&file.; 
	RUN;
%end;
%else %do;
	DATA mydata;
		SET &file.;
	RUN;
%end;

TITLE "Partitioned GMM";
PROC SORT DATA=mydata OUT=mydatasorted; 
BY &timeVar.; RUN;

*Check if there is a time independent variable;
%if &predVarTI. ^=.  %then %do;
	%let predVar = &predVarTD &predVarTI;
%end;
%else %do;
	%let predVar = &predVarTD;
%end;

*Use either the LWY or LS moment check;
%if &MC=LWY %then %do; *LWY approach;

*Obtain residuals;
%if &distr.=normal %then %do;
	proc reg data=mydatasorted NOPRINT;
	BY &timeVar.;
	MODEL &outVar. = &predVar.;
	OUTPUT OUT=outpool3 PREDICTED=mu;
	RUN;

	DATA outpool3; SET outpool3;
	wt = 1/mu; 
	rsdraw = &outVar.-mu; 
%end;
%else %if &distr.=bin %then %do;
	PROC logistic DATA=mydatasorted NOPRINT; 
	BY &timeVar.;
	MODEL &outVar. (event='1') = &predVar. / aggregate scale=none;
	OUTPUT OUT=outpool3 P=mu;
	RUN;

	DATA outpool3; SET outpool3;
	wt = mu*(1-mu);
	rsdraw = &outVar.-mu; 
%end;
%else %do;
	%put ERROR must be normal or binomial distributions; 
	%return;
%end;

PROC SORT DATA=outpool3 OUT=outpool3 ;
  BY &idVar. &timeVar.; RUN;
quit;

PROC IML;
use outpool3;                                                           
%if &predVarTI. ^=. %then %do;
	read all VARIABLES {&predVarTD. &predVarTI.} into Zmat; 
	read all var {&predVarTI.} into timeInd;
	read all var {&predVarTD.} into timeDep;
%end;
%else %do;
	read all VARIABLES {&predVarTD.} into Zmat;
%end;
read all var {wt} into wt;
read all var {rsdraw} into rsd;
read all var {&idVar.} INTO ID;
read all var {&timeVar.} INTO time;
close outpool3;


use mydata;
read all var{&predVarTD.} into timeDep2;
close mydata;


N=ncol(unique(ID)); *Number of subjects;
T = max(time); *Number of time points;
Np = ncol(Zmat); *Number of predictors, not including intercept;

NpTD2 = ncol(timeDep2);

%if &predVarTI. ^=. %then %do;
	NpTI = ncol(timeInd);*The last NPTI columns of Zmat will be time independent;
	NpTD = ncol(timeDep); *The first NPTD columns of Zmat will be time dependent;
%end;

*Find the correlation between X (a) and Y (rsd) across the T time points;
start rho(a,rsd) global(N,T);
abm = j(N,2*T,.);
abm[,1:T] = shape(rsd,N);	* N x T - First T columns are the values of rsd sorted into T columns (by time);
abm[,T+1:2*T] = shape(a,N); * Remaining T columns are the values of a sorted into T columns;
corr = corr(abm);  
rho = corr[1:T,T+1:2*T];    * T x T - Only take the correlations between the X and Y in the T time points;
return(rho);
finish rho;


*Standard deviation for each correlation;
start stddev(a,rsd) global(N,T);
bm = shape(rsd,N);    		 * N x T;
bdev = bm-j(N,1,1)*bm[:,];   * bdev N x T,   bm[:,] Col Mean is a row vector   1 x T - each row of bm minus column means;
bdev2 = bdev#bdev;      
am = shape(a,N);   
adev = am-j(N,1,1)*am[:,];  *N x T;
adev2 = adev#adev;      
stddev = sqrt( (1/N)*t(bdev2)*adev2 );   * T x T;
return(stddev);
finish stddev;

* corrected standardization;
start stdzn(x) global(N,T);
xrows = shape(x,N);   *N x T - by shape default columns are T=nrow(x)/N;
y = xrows - xrows[:,];  *N x T - Each value minus the column mean (mean for that time point);
vcv = (1/(N-1))*t(y)*y; * T x T;
v = diag(vcv); * T x T diagonal elements of vcv;
sinv = sqrt(inv(v));
x2 = y*sinv;   *N x T;
x2  = shape(x2,N*T,1); *N*T x 1 vector of standardized values;
return(x2);
finish stdzn;

pvec = j(Np*T*T,1,.);sevec = j(Np*T*T,1,.);  * pvec   (T*T) x Np;
r4out = j(T,T*Np,.); se4out = j(T,T*Np,.); z4out = j(T,T*Np,.); p4out = j(T,T*Np,.); 

y = rsd;
y_std = stdzn(y);


DO i=1 TO Np;
x = wt#Zmat[,i]; 			 * (N*T) x 1;
x_std = stdzn(x);

*Find p-values for the correlation of X_i and Y;
r = rho(x_std,y_std);		 * T x T;
se = stddev(x_std,y_std);	 * T x T;
z = sqrt(N)*(r/se);
p = 2*(1-cdf("normal",abs(z)));  * T x T;

*Fill the corresponding T columns of each matrix with the calculated values;
r4out[,T*(i-1)+1:T*i] = r;
se4out[,T*(i-1)+1:T*i] = se;
z4out[,T*(i-1)+1:T*i] = z;
p4out[,T*(i-1)+1:T*i] = p;

DO j = 1 TO T;
p[j,j] = 1;  * Not going to test diagonal elements, set to 1;
END;
*Takes the values in p4out and se4out, across row, then down column and creates vectors;
pvec[T*T*(i-1)+1:T*T*i,1] = shape(p,T*T,1);    *(T*T) x 1;
sevec[T*T*(i-1)+1:T*T*i,1] = shape(se,T*T,1);   *(T*T) x 1;
END;


TypeVec2 = (pvec >= &alpha.*j(Np*T*T,1,1) );

Type2 = shape(TypeVec2, Np, T*T);


*Individual test approach;
x = wt; 			 * (N*T) x 1;
x_std = stdzn(x);
r = rho(x_std,y_std);		 * T x T;
se = stddev(x_std,y_std);	 * T x T;
z = sqrt(N)*(r/se);
p = 2*(1-cdf("normal",abs(z)));  * T x T;

T_wt = p>&alpha.;
T_wt = shape(T_wt,1);

TypeMtx3 = j(Np+1,T*T,.);
TypeMtx3[1,] = T_wt;
TypeMtx3[2:Np+1,] = Type2;
Type2[1,] = shape(I(T),1);


*Adjust the last NpTI rows of TypeMtx3 to account for time independent covariates;
%if &predVarTI. ^=. %then %do;
	TypeMtx3= TypeMtx3[1:(NPTD+1),]; *Subset TypeMtx3 to include only time dependent;
	TypeMtxIND = repeat(shape(I(T),1),NPTI); *Create rows of typeMtx for the time independent;
	Np = NpTD;
%end;



/* define helper functions ROW and COL */
start row(x);  /* return matrix m such that m[i,j] = i */
   return( repeat( T(1:x), 1, x ));
finish;
start col(x);  /* return matrix m such that m[i,j] = j */
   return( repeat(1:x, x) );
finish;

*Helper matrices;
rt = row(T);
ct = col(T);

*Indices of upper diagonal for removal;
upperIdx = loc(ct>rt);
lowerIdx= loc(ct<rt);

*Remove backwards conditions;
typeMtxNB = TypeMtx3;
typeMtxNB[,upperIdx] = 0;


/*Keep backwards conditions*/
typeMtxfut=typeMtx3;
typeMtxFut[,loweridx]=0;


/*CHECKING HOW TO CHANGE ORDER OF COVARIATES FOR OUTCOMES AS PREDICTORS*/

CREATE TypeMtxFut from TypeMtxFut;
APPEND from TypeMtxFut;
CLOSE TypeMtxFut;


*Create the combination TypeMtx, of maximum size to be pared down later, excluding the intercept;
TypeMtxCombo = j((Np*T),T*T,0); 

TypeMtxCombo2 = j((Np*T),T*T,0);
*Loop across the different types, starting with current;
DO i = 0 TO (T-1);
	*Identify the indices for the appropriate setting;
	Idx = loc(rt-ct = i);
	IdxDif = setdif(1:T*T, Idx);/*THIS RETURNS ALL THOSE INDICES OF THE MATRIX FOR WHICH RT-CT IS DIFFERENT FROM I*/
	temp = typeMtxNB[2:(Np+1),];
	temp[,IdxDif] = 0;
	TypeMtxCombo[(i*Np+1):((i+1)*Np),] = temp;

	*Shift the values to accomodate shifted Y-X relationship;
	if(i>0) then DO; 
		dummyZero = j(Np, i, 0) || temp;
		TypeMtxCombo2[(i*Np+1):((i+1)*Np),] = dummyZero[,1:(T*T)];
	END;
	
	ELSE TypeMtxCombo2[(i*Np+1):((i+1)*Np),] = temp;

END; 


*Identify predictors or relationships with no valid moment conditions;

neq = TypeMtxCombo2[,+]; /*ADDING UP NUMBER OF 1'S FOR EACH ROW OF TYPEMTXCOMBO2*/


*Print a note about omitted moment conditions;
if min(neq)=0 then do;
	zeroPred = loc(neq=0);
	Note = "There are no valid moment conditions for " +char(ncol(zeroPred))+" covariate relationship(s).";
	print Note[label="Moment Condition Notes"];
	print "These covariate relationships will be omitted in the analysis.";
end;
else print "All covariate relationships will be evaluated."[label="Moment Condition Notes"];


/*FOR THOSE ROWS RELATED TO LAGS OF THE COVARIATE WITH NO VALID MOMENTS WE DELETE THAT ROW*/
/*THIS BECAUSE WE WILL NOT BE USING THOSE LAGS IN OUR MODELS*/
keepPred = loc(neq>0);

TypeMtxCombo3 = TypeMtxCombo2[keepPred,]; /*WILL ONLY KEEP ROWS OF COVARIATES WITH VALID MOMENTS*/


*Append the intercept conditions to the matrix;
*Append the time independent conditions to the matrix, if available;
%if &predVarTI. ^=. %then %do;
/*	TypeMtxCombo3 = TypeMtxNB[1,] // TypeMtxIND // TypeMtxCombo3; *Tested intercept;*/
	TypeMtxCombo3 = shape(I(T),1) // TypeMtxIND // TypeMtxCombo3; *Type I intercept;
	
%end;
%else %do;
/*	TypeMtxCombo3 = TypeMtxNB[1,] // TypeMtxCombo3; *Tested intercept;*/
	TypeMtxCombo3 = shape(I(T),1) // TypeMtxCombo3; *Type I intercept;
%end;

*If there are time independent variables, they are represented after the intercept, before the time dependent;
CREATE TypeMtxCombo3 from TypeMtxCombo3;
APPEND from TypeMtxCombo3;
CLOSE TypeMtxCombo3;
print "Each row of TypeMtx is the shifted type vector for each of the predictors, by individual test";
print TypeMtxCombo3;

*Create a printable form of the type matrix;
Type2[,upperIdx] = 0;
Type4out = j(T,T*(Np),.);
DO i=1 to Np;
Type4out[,T*(i-1)+1:T*i] = shape(Type2[i,],T,T);
END;

CREATE Type4out from Type4out;
APPEND from Type4out;
CLOSE Type4out;


*Create the modified data set;
USE mydatasorted;
read all var {&predVarTD.} into xnew;
read all var {&outVar.} into Ynew;
read all var {&idVar. &timeVar.} into othersnew;
%if &predVarTI. ^=. %then %do;
	read all VAR {&idVar. &outVar. &timeVar. &predVarTI.} INTO otherVars;
%end;
%else %do;
	read all VAR {&idVar. &outVar. &timeVar.} INTO otherVars;
%end;
CLOSE mydatasorted;


*Creating the new  variables;
X2 = j(N*T,Np*(T-1),0);/*X2 has dimensions [total # obs, #of lagged predictors multiplied by number of lags, all values start at 0*/
DO i=1 TO T-1;
	X2[(i*N+1):N*T, ((i-1)*Np+1):i*Np] = Xnew[1:(T-i)*N,];
END;


X2 = Xnew || X2;


*Remove the predictors with no valid moments;
X3 = X2[,keepPred];


*Add in the ID, outcome and time (and possible time independent variables;
X3 = otherVars || X3;


*Create the adjusted variable names;
predNames = t(repeat(t({&predVarTD.}),T));
lagName = j(1,Np*T,0);
DO i=1 to (T-1);
	lagName[1,(i*Np+1):(i+1)*Np]=j(1,Np,i);
END;
varNames = catx("_", predNames, char(lagName));
varnamesTD=varnames;


/*creating vector for time dependent predictors names*/
predNamesTD=t(repeat(t({&predVarTD.}),1));

/*this part will be deleted in feedback macro*/
create PredNamesTD from predNamesTD;
append from PrednamesTD;
close PredNamesTD;


/*creating vector for time dependent predictors distributions either  normal or binary*/

DistributionTD=t(repeat(t({&distrPredVarTD.}),1));

/*this part will be deleted in feedback macro*/
create DistributionTD from DistributionTD;
appedn from DistributionTD;
close DistributionTD;


*Retain only the variable names for those with valid moment conditions;
varNames = varNames[,keepPred];
%if &predVarTI. ^=. %then %do;
	varNames = {&predVarTI.} ||  varNames;
%end;

*Add in the names for the ID, outcome and time;
varNames2 = {&idVar.} || {&outVar.} || {&timeVar.} ||  varNames;

*Sort the data;
call sort(X3 , {1 3});

CREATE Mydata3 from X3[c=varNames2];
APPEND from X3;
close Mydata3;

*Module as demonstrated at http://blogs.sas.com/content/iml/2016/01/18/create-macro-list-values.html;
start CreateMacro(values, macroName, delimiter=' ');
   if type(values)="N" then          
      y = rowvec( char(values) );   /* convert numeric to character */
   else 
      y = rowvec(values);
   s = rowcat(y + delimiter);       /* delimit and concatenate */
   s = substr(s, 1, nleng(s)-nleng(delimiter)); /* remove delimiter at end */
   call symputx(macroName, s);      /* create macro variable */
finish;

call CreateMacro(varNames, "adjPredVar");

QUIT;

%end;


*Obtain initial GEE estimates for optimization;
ods select none;
proc genmod data=mydata3 descend;
	class &idVar. &timeVar.;
	model &outVar. = &adjPredVar. / dist=&distr.;
	repeated subject=&idVar. / within=&timeVar. corr=ind corrw;
	OUTPUT OUT=GEEout XBETA=xb RESRAW = rraw;
	ods output GEEEmpPEst=betaGEE;
RUN;
ods output close;
ods select all; 


*Obtain GMM Estimates;
PROC IML;
*Do not import los_2 since there are no moment conditions; 
USE mydata3;
READ all VAR {&adjPredVar.} INTO Zmat;
READ all VAR {&outVar.} INTO yvec;
READ all VAR {&idVar.} INTO ID;
READ all VAR {&timeVar.} INTO time;
CLOSE mydata3;

N=ncol(unique(ID)); *Number of subjects;
T = max(time); *Number of time points;
Np = ncol(Zmat); *Number of predictors, not including intercept;
Pn=Np+1;                      * number of covariates/parameters TO estimate, including intercept;

*Create the X matrix;
int = j(N*T,1,1);
Xmat =j(N*T,Pn,.); Xmat[,1]=int; Xmat[,2:Pn]=Zmat; *Xmat is matrix of intercept in column 1 (vector of 1) and Zmat;

*Use GEE estimates as starting values;
USE betaGEE;
READ all VAR {Estimate} INTO beta0;
CLOSE betaGEE;
beta0 = t(beta0);

*Use the identified covariate types from V40 code (multiple test);
USE TypeMtxCombo3;
READ all INTO TypeMtx;
CLOSE TypeMtxCombo3;

*Number of valid moment conditions;
neq = TypeMtx[,+];

nloc = j(1,Pn+1,0);          * nloc containing the starting/END positions of reg eqs brought by each covariate  ;
DO p =1 TO Pn;
  nloc[p+1] = sum(neq[1:p]); 
END;

nreg = sum(neq); * total number of reg_eq, 26;

Wn = I(nreg);                * I(26) initial weight matrix ;

*Module to calculate the objective function;
START TSGMM(beta) global(Pn,T,N,Xmat,yvec,nreg,TypeMtx,nloc,Wn);      
Gn = j(nreg,1,0);             * TO initialize Gn whenever objfun is called, b/c Gn depends on beta ;   
eq = j(nreg,N,0);              * 26x1625 reg_eqs for a type I covariate, for all persons  ;
S = j(nreg,nreg,0); 

DO i = 1 TO N;               * DO loops for each person;
  x = Xmat[(i-1)*T+1:i*T,];  * T x Pn, READ data for each person across all times;
  y = yvec[(i-1)*T+1:i*T];   * T outcomes for person i;

  if "&distr." = "bin" then 
  	mu = exp(x*t(beta)) / ( 1+exp(x*t(beta)) );
  else if "&distr." = "normal" then
 	mu = x*t(beta);

  Rsd = y - mu;                * T x 1;
  DO p = 1 TO Pn; * Do loops for each predictor;
  	if "&distr." = "bin" then
    	D = x[,p]#mu#(1- mu);
	else if "&distr." = "normal" then
		D= x[,p]#(mu##(-1));
    Eqmtx = Rsd*t(D);	*Gives x#mu#(1-mu)#(y-mu); *Do the time points not matter?;
    eq[nloc[p]+1:nloc[p+1],i] = Eqmtx[loc(TypeMtx[p,]^=0)];    *Save calculations to eq;
  END; *End loop across predictors;
  S = S + eq[,i]*t(eq[,i]);
END;     * TO CLOSE "DO i = 1 TO N";
Wn = ginv(S/N);		* Update the weight matrix;
Gn = eq[,:];               * row mean => col vector  neq x 1;
f = t(Gn)*Wn*Gn;           * the objective fn TO be minimized; 
RETURN(f);
FINISH TSGMM;


tc = {2000 2000}; optn = {0 1};              * optn[1]=0 specifies a minimization problem;
*Optimization takes: response code, row vector of optimal pts, module, starting value, optimization details, threshold conditions;
  *Run with the GEE estimates as starting parameters;
if "&optim."="NLPNRA" then
	CALL NLPNRA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson;
else if "&optim."="NLPNRR" then
	CALL NLPNRR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization Newton Raphson ridge;
else if "&optim."="NLPQN" then
	CALL NLPQN(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization dual quasi-newton;
else if "&optim."="NLPCG" then
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
else if "&optim."="NLPQUA" then
	CALL NLPQUA(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization quadratic;
else if "&optim."="NLPTR" then
	CALL NLPTR(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization trust region;
else if "&optim."="NLPNMS" then
	CALL NLPNMS(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization nelder-mead simplex;
else do;
	print "Must use valid optimization, will run with NLPCG";
	CALL NLPCG(rc, xres,"TSGMM", beta0,optn, , tc);           * DO nonlinear optimization conjugate gradient;
end;
beta = xres;


* Module to calculate the asymptotic variance;
DG = j(nreg,Pn,.);    * nreg x Pn;
DO k = 1 TO Pn; 
  DGi = j(nreg,N,0);          
  DO i = 1 TO N;               * DO loops for each for each person;
    x = Xmat[(i-1)*T+1:i*T,];  * T x Pn, READ data for each person across all times;
    y = yvec[(i-1)*T+1:i*T];  
	if "&distr." = "bin" then 
    	mu = exp(x*t(beta)) / ( 1+exp(x*t(beta)) );
	else if "&distr." = "normal" then
		mu = x*t(beta);
    Rsd = y - mu;                * T x 1;
	if "&distr." = "bin" then do;
    	Dk =  x[,k]#mu#(1- mu);
		Dkz =  x[,k]#(1- 2*mu);
		end;
	else if "&distr." = "normal" then do;
		Dk =  x[,k]#(mu##(-1));
		Dkz =  x[,k]#(-mu##(-2));
		end;
    DO p = 1 TO Pn;
	  if "&distr." = "bin" then 
      	Dp = x[,p]#mu#(1- mu);
	  else if "&distr." = "normal" then
		Dp = x[,p]#(mu##(-1));
      Dkzp = Dkz#Dp;                  * T x 1;
	  DGmtx = Rsd*t(Dkzp)-Dk*t(Dp);
      DGi[nloc[p]+1:nloc[p+1],i] = DGmtx[loc(TypeMtx[p,]^=0)];   
    END;
  END;
  DG[,k]= DGi[,:];    * row mean => col vector  neq x 1;
END;   



AsymVar = (1/N)*ginv(t(DG)*Wn*DG);   * Pn x Pn  note Wn is the inverse ;
AVvec = vecdiag(AsymVar); *Take only the diagonal elements for variance;
StdDev = sqrt(AVvec);

*Calculate the p-Values;
zvalue = t(beta)/StdDev;
pvalue = 2*(1-cdf('normal',abs(zvalue)));

*Create a matrix of results;
Outmtx = j(Pn,4,.);
ColLabel={"Estimate" "StdDev" "Zvalue" "Pvalue"};
Varnames=t({"Intercept"} || {&adjPredVar.});
OutTitle = "Analysis of Partial GMM Estimates";

Outmtx[,1]=t(beta);
Outmtx[,2]=StdDev;
Outmtx[,3]=zvalue;
Outmtx[,4]=pvalue;
PRINT Outmtx[colname=colLabel rowname=Varnames label=OutTitle ];

*Calculate residuals with GMM estimates;
resvec = yvec - Xmat*t(beta);

*Save the results as a dataset;
CREATE outdat1 FROM Outmtx;
APPEND FROM Outmtx;
CLOSE outdat1;
QUIT;
%mend partitionedGMMPart1; *End macro code;


