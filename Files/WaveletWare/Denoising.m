(* Mathematica package *)

BeginPackage["WaveletWare`Denoising`",{"WaveletWare`CommonUsage`","WaveletWare`MiscFunctions`","WaveletWare`Filters`","WaveletWare`OrthogonalTransforms`","WaveletWare`BiorthogonalTransforms`"}]

(* Exported symbols added here with SymbolName::usage *)

(* Denoising functions *)

MAD::usage="MAD[data] computes the median absolute deviation of input an input vector or matrix."
ShrinkageFunction::usage="ShrinkageFunction[t,lambda] returns the soft threshold value given an input t and tolerance lambda."
DonohoSURE::usage="DonohoSURE[data] returns the SUREShrink tolerance lambda of a given input vector or matrix."
UniversalThreshold::usage="UniversalThreshold[data,sigma] takes a vector or matrix and an estimate of the noise in the input and returns the universal threshold lambda."
NoiseEstimate::usage="NoiseEstimate[data] takes a vector or matrix and a wavelet filter and returns an estimate of the noise."
WhiteNoise::usage="WhiteNoise[dim,sigma] takes a length of a vector or the dimensions dim of a matrix and a noise level sigma and returns a vector/matrix of white noise."
TestSparseness::usage="TestSparseness[data] takes a vector or matrix and determines if it is sparse, relative to the SUREShrink algorithm."
WaveletShrinkage::usage="WaveletShrinkage[data,filter,lambda,opts] takes a vector or matrix, an orthogonal filter or biorthogonal filter pair, and a tolerance(s) and returns a denoised approximation of the input."
VisuShrink::usage="VisuShrink[data,filter,opts] takes a noisy vector or matrix and an orthogonal filter or biorthogonal filter pair and returns a denoised approximation of the input."
SUREShrink::usage="SUREShrink[data,filter,opts] takes a vector or matrix and a wavelet filter, an orthogonal filter or biorthogonal filter pair, and returns a denoised approximation of the input."


(* Option Values for image processing functions *)

ThresholdByLevel::usage="ThresholdByLevel is a flag that tells VisuShrink to use threshold values for each individual iteration of the wavelet transform instead of the default single value."
(* Options for denoising functions *)

Options[WaveletShrinkage]={NumIterations->1,Computation->Numerical,Shift->0,Boundary->False,IntegerMap->False};
Options[VisuShrink]={NumIterations->1,ThresholdByLevel->False,Computation->Numerical,Shift->0,Boundary->False,IntegerMap->False};
Options[SUREShrink]={NumIterations->1,ThresholdByLevel->False,Computation->Numerical,Shift->0,Boundary->False,IntegerMap->False};

Begin["`Private`"]
(* Implementation of the package *)

(* Denoising functions *)

MAD[v_]:=Median[Abs[Flatten[v]-Median[Flatten[v]]]];

ShrinkageFunction[t_, l_]:=Module[{},
	(* Check that l is a nonnegative number *)
	If[!(NumericQ[l]&&NonNegative[l]),Message[ShrinkageFunction::"badlambda"];Return[0]];
	Return[Sign[t]*Which[Abs[t]>=l,Abs[t]-l,True,0]];
];
SetAttributes[ShrinkageFunction,Listable];

DonohoSURE[v_]:=Module[{n,a,b,c,s,r,idx},
	If[!MemberQ[{1,2},First[CheckData[v]]],Message[DonohoSURE::"badinput"];Return[0]];
	n=Length[Flatten[v]];
	a=Sort[Flatten[v]^2];
    b=Drop[FoldList[Plus,0,a],1];
    c=Range[n-1,0,-1];
    s=b+c*a;
    r=n-2*Range[n]+s;
    idx=First[Flatten[Position[r,Min[r]]]];
    Return[Sqrt[a[[idx]]]];
];

UniversalThreshold[v_,sigma_]:=Module[{},
	If[!VectorQ[Flatten[v],NumericQ],Message[UniversalThreshold::"badinput"];Return[0]];
	Return[sigma*Sqrt[2*Log[Length[Flatten[v]]]]];
];

NoiseEstimate[v_]:=Module[{},
	If[!MemberQ[{1,2,3},First[CheckData[v]]],Message[NoiseEstimate::"badinput"];Return[0]];
	Return[MAD[Flatten[v]]/.6745];
];

WhiteNoise[n_,sigma_]:=Module[{},
	If[!(IntegerQ[n]&&Positive[n])&&!(VectorQ[n,(IntegerQ[#]&&Positive[#])&]&& Length[n]==2),Message[WhiteNoise::"baddimensions"];Return[{}]];
	If[!(NumericQ[sigma]&& sigma>0),Message[WhiteNoise::"badsigma"];Return[{}]];
	Return[RandomVariate[NormalDistribution[0,sigma],n]];
];

TestSparseness[a_]:=Module[{n,s,u},
	(* First check if the data are input as either a vector or a matrix *)
	If[!MemberQ[{1,2},First[CheckData[a]]],Message[TestSparseness::"badinput"];Return[{}]];
	n=Length[Flatten[a]]; 
	s=Total[Flatten[a]^2]/n-1;
	u=3*Log[2,n]/2/Sqrt[n];
	(* If input is sparse, return True, else False *)
	Return[If[s<=u,True,False]];
];

WaveletShrinkage[a_,h_,lambda_,OptionsPattern[]]:=Module[{datatype,maxits,its,comp,offset,boundary,intmap,wt,wtl,lp,hp,iwt},
	(* First check if the data are input as either a vector or a matrix *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{1,2},First[CheckData[a]]],Message[WaveletShrinkage::"badinput"];Return[{}]];
	
	(* Check validity of the filter or filter pair h *)
	If[!MemberQ[{1,2},CheckFilter[h]],Message[WaveletShrinkage::"badfilter"];Return[{}]];
	
	(* Read options and determine if they are valid*)
	{its,comp,offset,boundary,intmap}=OptionValue[WaveletShrinkage,{NumIterations,Computation,Shift,Boundary,IntegerMap}];
	{boundary,intmap}=Map[TrueFalse,{boundary,intmap}];
	
	(* Check the validity of the NumIterations value and the size of the data with regard to filter length and number of iterations of the wavelet transform. *)
	If[!CheckIterations[its,maxits],Message[WaveletShrinkage,"baditerations"];its=0];
	If[!CheckDataFilterSize[datatype,a,its,h],Message[WaveletShrinkage::"badsizes"];Return[{}]];
	
	(* Check validity of tolerance lambda *)
	If[!((NumericQ[lambda] && lambda>=0)||(VectorQ[lambda,NumericQ] && Length[lambda]==its && AllTrue[lambda,(#>=0)&])),Message[WaveletShrinkage::"badlambda"];Return[{}]];
	
	wt=Which[
		N[h]==={.5,.5},HWT[a,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->False],
		h==={1,1}/Sqrt[2],HWT[a,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->True],
		h===LeGall[],LWT[a,NumIterations->its,IntegerMap->intmap,Computation->comp],
		VectorQ[h],WT[a,h,NumIterations->its,Shift->offset,Computation->comp],
		True,BWT[a,h,NumIterations->its,Boundary->boundary,Computation->comp]
	];
	
	wtl=Which[VectorQ[a],WaveletVectorToList[wt,NumIterations->its],True,WaveletMatrixToList[wt,NumIterations->its]];
	{lp,hp}={First[wtl],Rest[wtl]};	
	
	hp=Which[VectorQ[lambda],MapThread[ShrinkageFunction[#1,#2]&,{hp,lambda}],True,Map[ShrinkageFunction[#,lambda]&,hp]];
	
	(* If the transform is LWT and IntegerMap is set to True, round the highpass data *)
	hp=If[intmap&&h===LeGall[],Round[hp],hp];
	
	wt=Which[VectorQ[a],Flatten[{lp,hp}],True,WaveletListToMatrix[Prepend[hp,lp],NumIterations->its]];
	
	iwt=Which[
		N[h]==={.5,.5},InverseHWT[wt,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->False],
		h==={1,1}/Sqrt[2],InverseHWT[wt,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->True],
		h===LeGall[],InverseLWT[wt,NumIterations->its,IntegerMap->intmap,Computation->comp],
		VectorQ[h],InverseWT[wt,h,NumIterations->its,Shift->offset,Computation->comp],
		True,InverseBWT[wt,h,NumIterations->its,Boundary->boundary,Computation->comp]
	];
	
	Return[iwt];
];

VisuShrink[a_,h_,OptionsPattern[]]:=Module[{datatype,maxits,its,lvl,comp,offset,boundary,intmap,wt,wtl,lp,hp,sigma,lambda,iwt},
	(* First check if the data are input as either a vector or a matrix *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{1,2},First[CheckData[a]]],Message[VisuShrink::"badinput"];Return[{}]];
	
	(* Check validity of the filter or filter pair h *)
	If[!MemberQ[{1,2},CheckFilter[h]],Message[VisuShrink::"badfilter"];Return[{}]];
	
	(* Read options and determine if they are valid*)
	{its,lvl,comp,offset,boundary,intmap}=OptionValue[VisuShrink,{NumIterations,ThresholdByLevel,Computation,Shift,Boundary,IntegerMap}];
	{lvl,boundary,intmap}=Map[TrueFalse,{lvl,boundary,intmap}];
	
	(* Check the validity of the NumIterations value and the size of the data with regard to filter length and number of iterations of the wavelet transform. *)
	If[!CheckIterations[its,maxits],Message[VisuShrink,"baditerations"];its=0];
	If[!CheckDataFilterSize[datatype,a,its,h],Message[VisuShrink::"badsizes"];Return[{}]];
	
	wt=Which[
		N[h]==={.5,.5},HWT[a,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->False],
		h==={1,1}/Sqrt[2],HWT[a,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->True],
		h===LeGall[],LWT[a,NumIterations->its,IntegerMap->intmap,Computation->comp],
		VectorQ[h],WT[a,h,NumIterations->its,Shift->offset,Computation->comp],
		True,BWT[a,h,NumIterations->its,Boundary->boundary,Computation->comp]
	];
	
	wtl=If[datatype==1,WaveletVectorToList[wt,NumIterations->its],WaveletMatrixToList[wt,NumIterations->its]];
	{lp,hp}={First[wtl],Rest[wtl]};	

		
	sigma=NoiseEstimate[Last[hp]];
	lambda=sigma*Sqrt[2*Log[Which[VectorQ[a],If[lvl,Map[Length,hp],Length[Last[hp]]],True,If[lvl,Map[Apply[Times,Dimensions[First[#]]]&,hp],Apply[Times, Dimensions[First[Last[hp]]]]]]]];
	
	hp=If[lvl,MapThread[ShrinkageFunction[#1,#2]&,{hp,lambda}],Map[ShrinkageFunction[#1,lambda]&,hp]];

	(* If the transform is LWT and IntegerMap is set to True, round the highpass data *)
	hp=If[intmap&&h===LeGall[],Round[hp],hp];

	wt=Which[VectorQ[a],Flatten[{lp,hp}],True,WaveletListToMatrix[Prepend[hp,lp],NumIterations->its]];
	
	iwt=Which[
		N[h]==={.5,.5},InverseHWT[wt,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->False],
		h==={1,1}/Sqrt[2],InverseHWT[wt,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->True],
		h===LeGall[],InverseLWT[wt,NumIterations->its,IntegerMap->intmap,Computation->comp],
		VectorQ[h],InverseWT[wt,h,NumIterations->its,Shift->offset,Computation->comp],
		True,InverseBWT[wt,h,NumIterations->its,Boundary->boundary,Computation->comp]
	];
	
	Return[iwt];
];

SUREShrink[a_,h_,OptionsPattern[]]:=Module[{datatype,maxits,its,comp,offset,boundary,intmap,wt,wtl,lp,hp,sigma,sparse,lambda,p,iwt},
	(* First check if the data are input as either a vector or a matrix *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{1,2},First[CheckData[a]]],Message[SUREShrink::"badinput"];Return[{}]];
	
	(* Check validity of the filter or filter pair h *)
	If[!MemberQ[{1,2},CheckFilter[h]],Message[SUREShrink::"badfilter"];Return[{}]];
	
	(* Read options and determine if they are valid*)
	{its,comp,offset,boundary,intmap}=OptionValue[SUREShrink,{NumIterations,Computation,Shift,Boundary,IntegerMap}];
	{boundary,intmap}=Map[TrueFalse,{boundary,intmap}];
	
	(* Check the validity of the NumIterations value and the size of the data with regard to filter length and number of iterations of the wavelet transform. *)
	If[!CheckIterations[its,maxits],Message[SUREShrink,"baditerations"];its=0];
	If[!CheckDataFilterSize[datatype,a,its,h],Message[VisuShrink::"badsizes"];Return[{}]];
	
	wt=Which[
		N[h]==={.5,.5},HWT[a,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->False],
		h==={1,1}/Sqrt[2],HWT[a,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->True],
		h===LeGall[],LWT[Round[a],NumIterations->its,IntegerMap->intmap,Computation->comp],
		VectorQ[h],WT[a,h,NumIterations->its,Shift->offset,Computation->comp],
		True,BWT[a,h,NumIterations->its,Boundary->boundary,Computation->comp]
	];
	
	wtl=Which[VectorQ[a],WaveletVectorToList[wt,NumIterations->its],True,WaveletMatrixToList[wt,NumIterations->its]];
	{lp,hp}={First[wtl],Rest[wtl]};	

	sigma=NoiseEstimate[Last[hp]];
	
	hp=hp/sigma;
	sparse=Map[TestSparseness,hp,Which[VectorQ[a],{1},True,{2}]];
	
	p=If[VectorQ[a],1,2];
	lambda=MapThread[If[#1,DonohoSURE[#2],UniversalThreshold[#2,sigma]]&,{sparse,hp},p];
	hp=sigma*MapThread[ShrinkageFunction[#1,#2]&,{hp,lambda},p];
	
	(* If the transform is LWT and IntegerMap is set to True, round the highpass data *)
	hp=If[intmap&&h===LeGall[],Round[hp],hp];

	wt=If[p==1,Flatten[{lp,hp}],WaveletListToMatrix[Prepend[hp,lp],NumIterations->its]];
	
	iwt=Which[
		N[h]==={.5,.5},InverseHWT[wt,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->False],
		h==={1,1}/Sqrt[2],InverseHWT[wt,NumIterations->its,Shift->offset,Computation->comp,Orthogonal->True],
		h===LeGall[],InverseLWT[wt,NumIterations->its,IntegerMap->intmap,Computation->comp],
		VectorQ[h],InverseWT[wt,h,NumIterations->its,Shift->offset,Computation->comp],
		True,InverseBWT[wt,h,NumIterations->its,Boundary->boundary,Computation->comp]
	];
	
	Return[iwt];
];


(* Error messages *)

DonohoSURE::badinput="The input must either be a vector or a matrix of numerical values - returning 0."
ShrinkageFunction::badlambda="Warning :: The value(s) of lambda must be a real nonnegative number - returning 0."
UniversalThreshold::badinput="The input must either be a vector or a matrix of numerical values - returning 0."
NoiseEstimate::badinput="The input must either be a vector or a matrix of numerical values - returning 0."
WhiteNoise::baddimensions="The input must either be a vector or a matrix of numerical values - WhiteNoise failed."
WhiteNoise::badsigma="The value of sigma must be a positive real number - WhiteNoise failed."
TestSparseness::badinput="The input must either be a vector or a matrix of numerical values - TestSparseness failed."
WaveletShrinkage::badinput="The input must either be a vector or a matrix - WaveletShrinkage failed."
WaveletShrinkage::badfilter="The input filter must be an orthogonal lowpass filter or a biorthogonal lowpass filter pair whose lengths are both even or both odd - WaveletShrinkage failed."
WaveletShrinkage::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
WaveletShrinkage::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - WaveletShrinkage failed."
WaveletShrinkage::badlambda="The value(s) of lambda must be a real nonnegative number."
VisuShrink::badinput="The input must either be a vector or a matrix - VisuShrink failed."
VisuShrink::badfilter="The input filter must be an orthogonal lowpass filter or a biorthogonal lowpass filter pair whose lengths are both even or both odd - VisuShrink failed."
VisuShrink::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
VisuShrink::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - VisuShrink failed."
SUREShrink::badinput="The input must either be a vector or a matrix - SUREShrink failed."
SUREShrink::badfilter="The input filter must be an orthogonal lowpass filter or a biorthogonal lowpass filter pair whose lengths are both even or both odd - SUREShrink failed."
SUREShrink::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
SUREShrink::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - SUREShrink failed."


End[]
EndPackage[]
