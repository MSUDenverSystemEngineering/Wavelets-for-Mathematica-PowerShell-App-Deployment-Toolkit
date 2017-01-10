(* Mathematica package *)

BeginPackage["WaveletWare`OrthogonalTransforms`",{"WaveletWare`CommonUsage`","WaveletWare`MiscFunctions`","WaveletWare`Filters`","WaveletWare`BiorthogonalTransforms`"}]

(* Exported symbols added here with SymbolName::usage *)

(* Orthogonal transform functions *)

WT::usage="WT[data,h,opts] computes the orthogonal wavelet transform of a vector, a matrix or a list of three matrices."
HWT::usage="HWT[data,opts] computes the Haar wavelet transform of a vector, a matrix or a list of three matrices."
LeftWT::usage="LeftWT[a,h,opts] computes Wa, where W is a discrete wavelet transform matrix and a is a matrix with even row dimension."
LeftHWT::usage="LeftHWT[a,opts] computes Wa, where W is a discrete wavelet transform matrix associated with the Haar filter and a is a matrix or list of three matrices each of the same dimensions.  The row dimension must be even."
RightWT::usage="RightWT[a,opts] computes aW^T where W is a discrete wavelet transform matrix and a is a matrix or list of three matrices each of the same dimensions.  The column dimension must be even."
RightHWT::usage="RightHWT[a,opts] computes aW^T, where W is a discrete wavelet transform matrix associated with the Haar filter and a is a matrix or list of three matrices each of the same dimensions.  The column dimension must be even."
WTFull::usage="WTFull[data,h,opts] takes a vector, a matrix or a list of three matrices and a (biorthogonal) filter and returns the (biorthogonal) wavelet transform for iterations 0,1,2,..."

InverseWT::usage="InverseWT[a,h,opts] computes the inverse orthogonal wavelet transform of a vector, a matrix or a list of three matrices of equal dimensions."
InverseHWT::usage="InverseHWT[a,opts] computes the inverse Haar wavelet transform of a vector, a matrix or a list of three matrices of equal dimensions."
LeftInverseWT::usage="LeftInverseWT[a,h,opts] computes W^Ta, where W is a discrete wavelet transform matrix and a is a wavelet-transformed matrix or list of three matrices each of the same dimension.  The row dimension must be even."
LeftInverseHWT::usage="LeftInverseHWT[a,opts] computes W^Ta, where W is a discrete Haar wavelet transform matrix and a is a wavelet-transformed matrix or list of three matrices each of the same dimension.  The row dimension must be even."
RightInverseWT::usage="RightInverseWT[a,h,opts] computes aW, where W is a discrete wavelet transform matrix and a is a wavelet-transformed matrix or list of three matrices each of the same dimension.  The column dimension must be even."
RightInverseHWT::usage="RightInverseHWT[a,opts] computes aW, where W is a discrete Haar wavelet transform matrix and a is a wavelet-transformed matrix or list of three matrices each of the same dimension.  The column dimension must be even."

(* Options for orthogonal transform functions *)

Options[WT]={NumIterations->1,Shift->0,Computation->Numerical};
Options[HWT]=Append[Options[WT],Orthogonal->False];
Options[LeftWT]={Shift->0,Computation->Numerical};
Options[LeftHWT]=Append[Options[LeftWT],Orthogonal->False];
Options[RightWT]=Options[LeftWT];
Options[RightHWT]=Options[LeftHWT];
Options[WTFull]={NumIterations->1,Shift->0,Boundary->False,IntegerMap->False,Computation->Numerical};

Options[InverseWT]=Options[WT];
Options[InverseHWT]=Options[HWT];
Options[LeftInverseWT]=Options[LeftWT];
Options[LeftInverseHWT]=Append[Options[LeftInverseWT],Orthogonal->False];
Options[RightInverseWT]=Options[RightWT];
Options[RightInverseHWT]=Append[Options[RightInverseWT],Orthogonal->False];

Begin["`Private`"]
(* Implementation of the package *)

(* Orthogonal transform functions *)

WT[a_,h_,OptionsPattern[]]:=Module[{datatype,maxits,its,offset,comp,p,tmp,g,f},
	
	(* Check the validity of the data *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{1,2,3},datatype],Message[WT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  1 = orthogonal filter *)
	If[CheckFilter[h]!=1,Message[WT::"badfilter"];Return[{}];];
	
	(* Read options and determine if they are valid*)
	{its,offset,comp}=OptionValue[WT,{NumIterations,Shift,Computation}];
	If[!CheckIterations[its,maxits],Message[WT::"baditerations",maxits];its=0;];
	If[!IntegerQ[offset],Message[WT::"badvalue","Shift",offset];offset=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,its,h],Message[WT::"badsizes"];Return[{}]];
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];
	
	(* Take care of the trivial case *)
	If[its==0,Return[p*a]];
	
	(* Do the vector case *)
	
	If[datatype==1,
	
		(* Create the highpass filter *)
		tmp=ConstantArray[1,Length[h]/2];
		g=Riffle[-tmp,tmp]*h;
		
		f[x_]:=Module[{y},
			y=Partition[RotateLeft[First[x],offset],Length[h],2,{1,2}];
			Return[{y.Reverse[h],y.g,Drop[x,1]}];
		];
		
		Return[Flatten[Nest[f,{p*a},its]]];
	];
	
	(* Here is a function for performing iterations of the wavelet transform on a matrix *)
	
	f[x_]:=Module[{r},
		r[t_]:=Join[WaveletMatrixToList[RightWT[LeftWT[First[t],h,Shift->offset,Computation->comp],h,Shift->offset,Computation->comp]],Drop[t,1]];
		Return[WaveletListToMatrix[Nest[r,x,its],NumIterations->its]];
	];
	
	If[datatype==2,
		Return[f[{p*a}]];
	];
	
	If[datatype==3,
		Return[Map[f[{#}]&,p*a]];
	];
];

HWT[a_,OptionsPattern[]]:=Module[{filter,its,offset,comp},
	filter=TrueFalse[OptionValue[HWT,Orthogonal]];
	{its,offset,comp}=OptionValue[HWT,{NumIterations,Shift,Computation}];
	filter=If[filter,Haar[],{1,1}/2];
	Return[WT[a,filter,NumIterations->its,Shift->offset,Computation->comp]];
];

LeftWT[a_,h_,OptionsPattern[]]:=Module[{datatype,maxits,offset,comp,p,g,tmp,f},
	
	(* Check the validity of the data *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[LeftWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  1 = orthogonal filter *)
	If[CheckFilter[h]!=1,Message[LeftWT::"badfilter"];Return[{}];];
	
	(* Read options and determine if they are valid*)
	{offset,comp}=OptionValue[LeftWT,{Shift,Computation}];
	If[!IntegerQ[offset],Message[LeftWT::"badvalue","Shift",offset];offset=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,h],Message[LeftWT::"badsizes"];Return[{}]];
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	tmp=ConstantArray[1,Length[h]/2];
	g=Riffle[-tmp,tmp]*h;
	
	f[x_]:=Module[{r},
		r[t_]:=Flatten[Join[Map[(Partition[RotateLeft[t,offset],Length[h],2,{1, 2}].#)&,{Reverse[h],g}]]];
		Return[Transpose[Map[r,Transpose[p*x]]]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

LeftHWT[a_,OptionsPattern[]]:=Module[{filter,offset,comp},
	filter=TrueFalse[OptionValue[LeftHWT,Orthogonal]];
	{offset,comp}=OptionValue[LeftHWT,{Shift,Computation}];
	filter=If[filter,Haar[],{1,1}/2];
	Return[LeftWT[a,filter,Shift->offset,Computation->comp]];
];

RightWT[a_,h_,OptionsPattern[]]:=Module[{datatype,maxits,offset,comp,p,g,tmp,f},
	
	(* Check the validity of the data *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[RightWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  1 = orthogonal filter *)
	If[CheckFilter[h]!=1,Message[RightWT::"badfilter"];Return[{}];];
	
	(* Read options and determine if they are valid*)
	{offset,comp}=OptionValue[RightWT,{Shift,Computation}];
	If[!IntegerQ[offset],Message[RightWT::"badvalue","Shift",offset];offset=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,h],Message[LeftWT::"badsizes"];Return[{}]];
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	tmp=ConstantArray[1,Length[h]/2];
	g=Riffle[-tmp,tmp]*h;
	
	f[x_]:=Module[{r},
		r[t_]:=Flatten[Join[Map[(Partition[RotateLeft[t,offset],Length[h],2,{1,2}].#)&,{Reverse[h],g}]]];
		Return[Map[r,p*x]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

RightHWT[a_,OptionsPattern[]]:=Module[{filter,offset,comp},
	filter=TrueFalse[OptionValue[RightHWT,Orthogonal]];
	{offset,comp}=OptionValue[RightHWT,{Shift,Computation}];
	filter=If[filter,Haar[],{1,1}/2];
	Return[RightWT[a,filter,Shift->offset,Computation->comp]];
];

WTFull[v_,h_,OptionsPattern[]]:=Module[{datatype,its,offset,boundary,intmap,comp,maxits,f,wt},
	
	(* Check the validity of the data *)
	{datatype,maxits}=CheckData[v];
	If[!MemberQ[{1,2,3},datatype],Message[WTFull::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  1 = orthogonal filter, 2 = biorthogonal filter *)
	If[!MemberQ[{1,2},CheckFilter[h]],Message[WTFull::"badfilter"];Return[{}];];
	
	(* Get options and check integrity *)
	{its,offset,boundary,intmap,comp}=OptionValue[WTFull,{NumIterations,Shift,Boundary,IntegerMap,Computation}];
	{boundary,intmap}=Map[TrueFalse,{boundary,intmap}];	
	If[!CheckIterations[its,maxits],Message[WTFull::"baditerations",maxits];its=0;];
	If[!IntegerQ[offset],Message[WTFull::"badvalue","Shift",offset];offset=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,v,its,h],Message[WTFull::"badsizes"];Return[{}]];
	
	(* Compute the transformations *)
	f[t_,i_]:=Which[
		N[h]==={.5,.5},HWT[t,NumIterations->i,Shift->offset,Computation->comp,Orthogonal->False],
		h==={1,1}/Sqrt[2],HWT[t,NumIterations->i,Shift->offset,Computation->comp,Orthogonal->True],
		h===LeGall[],LWT[t,NumIterations->i,IntegerMap->intmap,Computation->comp],
		VectorQ[h],WT[t,h,NumIterations->i,Shift->offset,Computation->comp],
		True,BWT[t,h,NumIterations->i,Boundary->boundary,Computation->comp]
	];
	wt=Map[f[v,#]&,Range[0,its]];
	(*Return[MapThread[WaveletVectorToList[#1,NumIterations->#2]&,{wt,Range[0,its]}]];*)
	Return[wt];
];

InverseWT[a_,h_,OptionsPattern[]]:=Module[{datatype,its,offset,comp,maxits,p,tmp,g,data,f},
	
	(* Check the validity of the data *)
	{datatype,maxits}=CheckData[a];
	
	If[!MemberQ[{1,2,3},datatype],Message[InverseWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  1 = orthogonal filter *)
	If[CheckFilter[h]!=1,Message[InverseWT::"badfilter"];Return[{}];];
	
	(* Read options and determine if they are valid*)
	{its,offset,comp}=OptionValue[InverseWT,{NumIterations,Shift,Computation}];
	If[!CheckIterations[its,maxits],Message[InverseWT::"baditerations",maxits];its=0;];
	If[!IntegerQ[offset],Message[InverseWT::"badvalue","Shift",offset];offset=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,its,h],Message[InverseWT::"badsizes"];Return[{}]];
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];
	
	(* Take care of the trivial case *)
	If[its==0,Return[p*a]];
	
	(* Do the vector transformation *)
	If[datatype==1,
		
		(* Create the highpass filter *)
		tmp=ConstantArray[1,Length[h]/2];
		g=Riffle[-tmp,tmp]*h;
		
		data=WaveletVectorToList[p*a,NumIterations->its];
		f[t_] := Module[{s, q, r},
   			s=Map[Partition[#,Length[h]/2,1,{Length[h]/2,Length[h]/2}]&,Take[t,2]];
   			q=Map[Reverse[Transpose[Partition[#,2]]]&,{h,Reverse[g]}];
   			r=MapThread[Apply[Riffle,Transpose[#1.Transpose[#2]]]&,{s, q}];
   			Return[Prepend[Drop[t,2],RotateRight[First[r]+Last[r],offset]]];
   		];
		Return[First[Nest[f,data,its]]];
	];
	
	f[x_]:=Module[{d,r},
		d=WaveletMatrixToList[x,NumIterations->its];
		
		r[t_]:=Module[{y},
			y=RightInverseWT[LeftInverseWT[WaveletListToMatrix[Take[t,2]],h,Shift->offset,Computation->comp],h,Shift->offset,Computation->comp];
			Return[Prepend[Drop[t,2],y]];
		];
		
		Return[First[Nest[r,d,its]]];
	];
		
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];


InverseHWT[a_,OptionsPattern[]]:=Module[{filter,its,offset,comp},
	filter=TrueFalse[OptionValue[InverseHWT,Orthogonal]];
	{its,offset,comp}=OptionValue[InverseHWT,{NumIterations,Shift,Computation}];
	
	filter=If[filter,Haar[],{1,1}];
	Return[InverseWT[a,filter,NumIterations->its,Shift->offset,Computation->comp]];
];

LeftInverseWT[a_,h_,OptionsPattern[]]:=Module[{datatype,maxits,offset,comp,p,tmp,g,f},
	
	(* Check the validity of the data *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[LeftInverseWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  1 = orthogonal filter *)
	If[CheckFilter[h]!=1,Message[LeftInverseWT::"badfilter"];Return[{}];];
	
	(* Read options and determine if they are valid*)
	{offset,comp}=OptionValue[LeftInverseWT,{Shift,Computation}];
	If[!IntegerQ[offset],Message[LeftInverseWT::"badvalue","Shift",offset];offset=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,h],Message[LeftInverseWT::"badsizes"];Return[{}]];
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	tmp=ConstantArray[1,Length[h]/2];
	g=Riffle[-tmp,tmp]*h;
	
	(* A function to apply the left transform to a matrix *)
	f[x_]:=Module[{r,data},
		r[t_]:=Module[{s,q,v},
   			s=Map[Partition[#,Length[h]/2,1,{Length[h]/2,Length[h]/2}]&,Take[t,2]];
   			q=Map[Reverse[Transpose[Partition[#,2]]]&,{h,Reverse[g]}];
   			v=MapThread[Apply[Riffle,Transpose[#1.Transpose[#2]]]&,{s, q}];
   			Return[RotateRight[First[v]+Last[v],offset]];
   		];
   		data=Map[Partition[#,Length[#]/2]&,Transpose[p*x]];
   		Return[Transpose[Map[r,data]]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

LeftInverseHWT[a_,OptionsPattern[]]:=Module[{filter,offset,comp},
	filter=TrueFalse[OptionValue[LeftInverseHWT,Orthogonal]];
	{offset,comp}=OptionValue[LeftInverseHWT,{Shift,Computation}];
	filter=If[filter,Haar[],{1,1}];
	Return[LeftInverseWT[a,filter,Shift->offset,Computation->comp]];
];

RightInverseWT[a_,h_,OptionsPattern[]]:=Module[{datatype,maxits,offset,comp,p,tmp,g,f},
	
	(* Check the validity of the data *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[RightInverseWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  1 = orthogonal filter *)
	If[CheckFilter[h]!=1,Message[RightInverseWT::"badfilter"];Return[{}];];
	
	(* Read options and determine if they are valid*)
	{offset,comp}=OptionValue[RightInverseWT,{Shift,Computation}];
	If[!IntegerQ[offset],Message[RightInverseWT::"badvalue","Shift",offset];offset=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,h],Message[RightInverseWT::"badsizes"];Return[{}]];
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	tmp=ConstantArray[1,Length[h]/2];
	g=Riffle[-tmp,tmp]*h;
	
	
	(* A function to apply the left transform to a matrix *)
	f[x_]:=Module[{r,data},
		r[t_]:=Module[{s,q,v},
   			s=Map[Partition[#,Length[h]/2,1,{Length[h]/2,Length[h]/2}]&,Take[t,2]];
   			q=Map[Reverse[Transpose[Partition[#,2]]]&,{h,Reverse[g]}];
   			v=MapThread[Apply[Riffle,Transpose[#1.Transpose[#2]]]&,{s, q}];
   			Return[RotateRight[First[v]+Last[v],offset]];
   		];
   		data=Map[Partition[#,Length[#]/2]&,p*x];
   		Return[Map[r,data]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

RightInverseHWT[a_,OptionsPattern[]]:=Module[{filter,offset,comp},
	filter=TrueFalse[OptionValue[RightInverseHWT,Orthogonal]];
	{offset,comp}=OptionValue[RightInverseHWT,{Shift,Computation}];
	filter=If[filter,Haar[],{1,1}];
	Return[RightInverseWT[a,filter,Shift->offset,Computation->comp]];
];

(* Error messages *)

WT::badinput="The input for WT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - WT failed."
WT::badfilter="The input filter must be an orthogonal lowpass filter of even length - WT failed."
WT::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
WT::badvalue="Warning :: The value `2` for `1` is not an integer - resetting `1` = 0." 
WT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - WT failed."
LeftWT::badinput="The input for LeftWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - LeftWT failed."
LeftWT::badfilter="The input filter must be an orthogonal lowpass filter of even length - LeftWT failed."
LeftWT::badvalue="Warning :: The value for `1`, `2`, is not an integer - resetting `1` to its default value." 
LeftWT::baddimensions="The `1` is not longer than the filter length - LeftWT failed."
LeftWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - LeftWT failed."
RightWT::badinput="The input for RightWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - RightWT failed."
RightWT::badfilter="The input filter must be an orthogonal lowpass filter of even length - RightWT failed."
RightWT::badvalue="Warning :: The value for `1`, `2`, is not an integer - resetting `1` to its default value." 
RightWT::baddimensions="The `1` is not longer than the filter length - RightWT failed."
RightWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - RightWT failed."
WTFull::badinput="The input for WTFull is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - WTFull failed."
WTFull::badfilter="The input filter must be an orthogonal lowpass filter of even length - WTFull failed."
WTFull::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
WTFull::badvalue="Warning :: The value `2` for `1` is not an integer - resetting `1` = 0." 
WTFull::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - WTFull failed."

InverseWT::badinput="The input for InverseWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - InverseWT failed."
InverseWT::badfilter="The input filter must be an orthogonal lowpass filter of even length - InverseWT failed."
InverseWT::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
InverseWT::badvalue="Warning :: The value `2` for `1` is not an integer - resetting `1` = 0." 
InverseWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - InverseWT failed."
LeftInverseWT::badinput="The input for LeftInverseWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - LeftInverseWT failed."
LeftInverseWT::badfilter="The input filter must be an orthogonal lowpass filter of even length - LeftInverseWT failed."
LeftInverseWT::badvalue="Warning :: The value `2` for `1` is not an integer - resetting `1` = 0." 
LeftInverseWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - LeftInverseWT failed."
RightInverseWT::badinput="The input for RightInverseWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - RightInverseWT failed."
RightInverseWT::badfilter="The input filter must be an orthogonal lowpass filter of even length - RightInverseWT failed."
RightInverseWT::badvalue="Warning :: The value `2` for `1` is not an integer - resetting `1` = 0." 
RightInverseWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - RightInverseWT failed."

End[]
EndPackage[]
