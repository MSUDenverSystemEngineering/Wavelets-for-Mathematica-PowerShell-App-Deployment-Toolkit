(* Mathematica package *)

BeginPackage["WaveletWare`Filters`",{"WaveletWare`CommonUsage`","WaveletWare`MiscFunctions`"}]

(* Exported symbols added here with SymbolName::usage *)

(* Function Usage *)
Haar::usage="Haar[ ] returns the Haar lowpass filter."
Daub::usage="Daub[n] returns the length n Daubechies lowpass filter."
Coif::usage="Coif[k] returns the length 6k Coiflet lowpass filter."
LeGall::usage="LeGall[ ] returns the LeGall lowpass filter pair."
SplineFilters::usage="SplineFilters[m,n] returns biorthogonal spline lowpass filter pairs of length 2m+n-1 and n+1, respectively."
CDF97::usage="CDF97[ ] returns the 9/7 Cohen, Daubechies, Feauveau biorthogonal lowpass filter pair."

(* Options Usage *)

(* Options for filter functions *)
Options[Daub]={WorkingPrecision->$MachinePrecision};
Options[Coif]=Options[Daub];
Options[SplineFilters]={DisplayInfo->False};


Begin["`Private`"]
(* Implementation of the package *)

(* Functions *)
Haar[]=Sqrt[2]*{1,1}/2;

Daub[n_,OptionsPattern[]]:=Module[{prec,n2,p,y,k,alpha,r,z,f},
	If[!(IntegerQ[n]&&n>0&&EvenQ[n]),Message[Daub::"badinteger",n];Return[{}];];
	n2=n/2;
	
	prec=OptionValue[Daub,WorkingPrecision];
	If[!(NumberQ[prec]&&prec>0),Message[Daub::"badprecision"];prec=$MachinePrecision;];
	
	If[n2==1,Return[SetPrecision[Haar[],prec]]];
	y[x_]:=(1-x)*(1-1/x)/4;
	p[x_]:=Sum[Binomial[n-1,k]*y[x]^k*(1-y[x])^(n2-1-k),{k,0,n2-1}];
	r=Select[z/.NSolve[p[z]==0,z,WorkingPrecision->prec],(Norm[#]<1)&];
	alpha=Coefficient[p[z],z,Exponent[p[z],z]];
	f[x_]:=Sqrt[Abs[alpha]]*Apply[Times,Map[(x-#)&,r]]/Apply[Times,Sqrt[r]];
	Return[Reverse[Chop[Sqrt[2]*CoefficientList[((1+z)/2)^n2*f[z],z]]]];
];

Coif[k_,OptionsPattern[]]:=Module[{prec,idx,hh,h,H,p,orth,j,derivs0,derivsPi,s},
	If[!(IntegerQ[k] && k>0 && k<6),Message[Coif::"badinteger",k];Return[{}]];
	
	prec=OptionValue[Coif,WorkingPrecision];
	If[!(NumberQ[prec]&&prec>0),Message[Coif::"badprecision"];prec=$MachinePrecision;];
	
	idx = {2,1,5,6,12}[[k]];
	hh = Array[h, 6 k,-2 k];
	H[w_]:=hh.Table[E^(I*j*w),{j,-2 k,4 k - 1}];
	p=PadRight[hh, 2*Length[hh]];
	orth=Table[p.RotateRight[p, 2*j]==DiscreteDelta[0, j],{j,0,Length[hh]/2-1}];
	derivs0=Table[Derivative[j][H][0]==Sqrt[2]*DiscreteDelta[0,j],{j,0,2 k - 1}];
	derivsPi=Table[Derivative[j][H][Pi]==0,{j, 0, 2 k - 1}];
	s=Map[(hh/.#)&,NSolve[Join[orth,derivs0,derivsPi],hh,WorkingPrecision->prec]];
	Return[s[[idx]]];
];

LeGall[]:={{-1,2,6,2,-1}/8,{1,2,1}/2};

SplineFilters[n_,nt_,OptionsPattern[]]:=Module[{parity,l,lt,p,len,start,stop,t,h,ht,k},
	If[!AllTrue[{n, nt}, #>0&],Message[SplineFilters::"badinput","Input values cannot be negative."];Return[{}]];
	parity = Which[AllTrue[{n,nt},EvenQ],0,AllTrue[{n,nt},OddQ],1,True,Message[SplineFilters::"badinput","Input values must have the same parity."];Return[{}]];
	
	{l,lt}=({n, nt}-parity)/2;
	p[w_]:= Which[parity==0,1,True,E^(I*w/2)]*TrigToExp[Cos[w/2]^n*Sum[Binomial[l+lt+k+parity-1,k]*TrigToExp[Sin[w/2]^(2*k)],{k,0,l+lt+parity-1}]];
	len=2*n+nt-1;
	{start, stop}={(1+parity-len)/2,(len-1+parity)/2};
	h=Sqrt[2]*Table[Coefficient[p[t], E^(I*t),k], {k, start, stop}];
	ht=Sqrt[2]*Map[Binomial[nt,#]&,Range[0,nt]]/2^nt;
	
	(* Display Information if requested. *)
	If[TrueFalse[OptionValue[SplineFilters,DisplayInfo]],
		Print["Biorthogonal Spline Filter Information: \n"];
		Print["The length of the first filter is ",len," with starting and stopping indices ",start," and ",stop,", respectively."];
		{start,stop}={(parity-nt)/2,(nt+parity)/2};
		Print["The length of the second filter is ",nt+1," with starting and stopping indices ",start," and ",stop,", respectively."];
	];
	Return[{h,ht}];
];

CDF97[]:=Module[{rts, y,p, pw, H, Hw, a, w, h, ht},
  rts=Map[(y/.#)&,NSolve[1+4 y+10 y^2+20 y^3==0,y]];
  p[t_]:=a*(t-First[rts]);
  pw[t_]:=20*Apply[Times,Map[(t - #)&,Drop[rts, 1]]]/a;
  H[t_]:=Sqrt[2]*TrigToExp[Cos[t/2]^4]*p[TrigToExp[Sin[t/2]^2]]; 
  Hw[t_]:=Sqrt[2]*TrigToExp[Cos[t/2]^4]*pw[TrigToExp[Sin[t/2]^2]];
  Clear[a];
  a=First[a/.NSolve[H[0]==Sqrt[2],a]];
  h=Chop[Table[Coefficient[Hw[w],E^(I*w),k],{k,-4,4}]];
  ht=Chop[Table[Coefficient[H[w],E^(I*w),k],{k,-3,3}]];
  Return[{h,ht}];
  ];
  
(* Error Messages *)

Daub::badinteger="The value `1` is not an even, positive integer - Daub failed."
Daub::badprecision="Warning :: The value for working precision must be a positive number - resetting the working precision to the value of $MachinePrecision."
Coif::badinteger="The value `1` is not a member of {1,2,3,4,5} - Coif failed."
Coif::badprecision="Warning :: The value for working precision must be a positive number - resetting the working precision to the value of $MachinePrecision."
SplineFilters::badinput="`1` - SplineFilters failed."

  
  

End[]

EndPackage[]