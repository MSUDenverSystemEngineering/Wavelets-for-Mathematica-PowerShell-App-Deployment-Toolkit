(* Mathematica package *)

BeginPackage["WaveletWare`WaveletFunctions`",{"WaveletWare`CommonUsage`"}]

(* Exported symbols added here with SymbolName::usage *)

(* Continuous wavelet functions *)

ScalingFunctionData::usage="ScalingFunctionData[h,opts] takes a lowpass filter and applies the cascade algorithm to create a discrete vector of samples of the scaling function associated with the filter."
WaveletFunctionData::usage="WaveletFunctionData[h,opts] or WaveletFunctionData[{h,ht},opts] takes a lowpass filter (or filter pair), creates the associated highpass filter and then applies the cascade algorithm to create a discrete vector of samples of the wavelet function associated with the filter."
PacketFunctionData::usage="PacketFunctionData[h,n,opts] or PacketFunctionData[{h,ht},n,opts] takes a lowpass filter (pair) and a nonnegative integer n and creates the nth wavelet packet function."

(* Options *)

Resolution::usage="Resolution is an option for the ScalingFunctionData, WaveletFunctionData and PacketFunctionData to indicate the number of sample points to use per integer length of the support."
FilterStart::usage="FilterStart is an option for ScalingFunctionData, WaveletFunctionData and PacketFunctionData to identify the starting index of the input filter."

(* Options for continuous wavelet functions *)

Options[ScalingFunctionData]={NumIterations->10,Resolution->50,FilterStart->0};
Options[WaveletFunctionData]=Options[ScalingFunctionData];
Options[PacketFunctionData]=Options[ScalingFunctionData];

Begin["`Private`"]
(* Implementation of the package *)

ScalingFunctionData[h_,OptionsPattern[]]:=Module[{its,ppi,start,supp,v0,v1,pre,post,unitsize,k},
	
	(* Check the integrity of the filter *)
	If[!VectorQ[h,NumericQ],Message[ScalingFunctionData::"badfilter"];Return[{}];];
	
	(* Get option values and check integrity *)
	{its,ppi,start}=OptionValue[ScalingFunctionData,{NumIterations,Resolution,FilterStart}];
	
	If[!IntegerQ[its] || Negative[its],Message[ScalingFunctionData::"iterations"];its=NumIterations/.Options[ScalingFunctionData]];
	If[!IntegerQ[ppi] || ppi<1,Message[ScalingFunctionData::"resolution"];ppi=Resolution/.Options[ScalingFunctionData];];
	If[!IntegerQ[start],Message[ScalingFunctionData::"start"];start=FilterStart/.Options[ScalingFunctionData];];
		
	supp=Length[h]-1;
	
	(* Take care of the trivial case *)
	If[its==0,
		Return[Transpose[{Most[Range[start, supp + start, 1/ppi]],ConstantArray[1,ppi*supp]}]];
	];
	
	(* For iterations bigger than 1: *)
	unitsize=supp*ppi;
	pre=ConstantArray[0,unitsize];
	post=ConstantArray[0,2*unitsize-ppi-1];
	v0=Flatten[{pre,ConstantArray[1,ppi],post}];
	post=Drop[post,unitsize-ppi];
	For[k=1,k<=its,k++,
		v1=Sqrt[2.]*h.Transpose[Table[Reverse[Take[RotateLeft[v0,2*k],{1,unitsize+supp,ppi}]],{k,0,unitsize-1}]];
		v0=Flatten[{pre,v1,post}];
	];
	v0=start+Range[0,unitsize-1]/ppi; (* Table[N[k/ppi],{k,0,unitsize-1}]; *)
	Return[Transpose[{v0,v1}]];
];

WaveletFunctionData[h_,OptionsPattern[]]:=Module[{its,res,start,data,stop,low,a,b,g,n,nt,k,x,y},
	
	(* Check validity of the filter or filter pair h *)
	If[!((VectorQ[h] && EvenQ[Length[h]]) || (AllTrue[h,VectorQ] && Length[h]==2 && (AllTrue[Map[Length,h],EvenQ] || AllTrue[Map[Length,h],OddQ]))),Message[WaveletFunctionData::"badfilter"];Return[{}]];
	
	(* Read in the options and check integrity *)
	{its,res,start}=OptionValue[WaveletFunctionData,{NumIterations,Resolution,FilterStart}];
	If[!IntegerQ[its] || Negative[its],Message[WaveletFunctionData::"iterations"];its=NumIterations/.Options[WaveletFunctionData]];
	If[!IntegerQ[res]||!Positive[res],Message[WaveletFunctionData::"resolution"];res=Resolution/.Options[WaveletFunctionData];];
	If[!IntegerQ[start],Message[WaveletFunctionData:"start"];Return[{}]];
	
	If[VectorQ[h],
		nt=Length[h];
		stop=nt+start-1;
		low=h;
		g=Table[(-1)^k,{k,start,stop}]*Reverse[h];
		{a,b}={1-nt/2,nt/2},
		{n,nt}=Map[Length,h];
		{start,stop}=If[OddQ[nt],{(3-nt)/2,(nt+1)/2},{1-nt/2,nt/2}];
		low=First[h];
		{a,b}=Transpose[Map[If[OddQ[#],{(1-#)/2,(#-1)/2},{1-#/2,#/2}]&,{n,nt}]];
		{a,b}={(First[a]-Last[b]+1)/2,(First[b]-Last[a]+1)/2};
		g=(-1)^nt*Table[(-1)^k,{k,start,stop}]*Reverse[Last[h]];
	];
	
	x = Most[Range[a,b,1/res]];
	data=ScalingFunctionData[low,NumIterations->its,Resolution->res,FilterStart->start];
	
	(* Create the ordinates and downsample to get the resolution correct *) 
	y = Sqrt[2]*Total[g*Table[RotateRight[Join[Map[Last, data],ConstantArray[0,res*(nt-1)]],k*res],{k,0,nt-1}]];
	y=Take[y,{1,Length[y],2}];
	
	Return[Transpose[{x,y}]];
];

PacketFunctionData[h_,m_,OptionsPattern[]]:=Module[{its,res,start,k,low,g,n,nt,n1,n2,a,b,adj,x,y,data,filter},
	
	(* Check validity of the filter or filter pair h *)
	If[!((VectorQ[h] && EvenQ[Length[h]]) || (AllTrue[h,VectorQ] && Length[h]==2 && (AllTrue[Map[Length,h],EvenQ] || AllTrue[Map[Length,h],OddQ]))),Message[PacketFunctionData::"badfilter"];Return[{}]];
	
	(* Check that n is a nonnegative integer *)
	If[!IntegerQ[m] || Negative[m],Message[PacketFunctionData:"badinteger"];m=0];
	
	(* Read in the options and check integrity *)
	{its,res,start}=OptionValue[PacketFunctionData,{NumIterations,Resolution,FilterStart}];
	If[!IntegerQ[its] || Negative[its],Message[PacketFunctionData::"iterations"];its=NumIterations/.Options[PacketFunctionData]];
	If[!IntegerQ[res]||!Positive[res],Message[PacketFunctionData::"resolution"];res=Resolution/.Options[PacketFunctionData];];
	If[!IntegerQ[start],Message[PacketFunctionData:"start"];Return[{}]];
	
	If[VectorQ[h],
		nt=Length[h];
		{n1,n2}={start,nt+start-1};
		low=h;
		g=Table[(-1)^k,{k,n1,n2}]*Reverse[h];
		{a,b}={1-nt/2,nt/2},
		{n,nt}=Map[Length,h];
		{n1,n2}=If[OddQ[nt],{(3-nt)/2,(nt+1)/2},{1-nt/2,nt/2}];
		g=(-1)^nt*Table[(-1)^k,{k,n1,n2}]*Reverse[Last[h]];
		low=First[h];
		{a,b}=Transpose[Map[If[OddQ[#],{(1-#)/2,(#-1)/2},{1-#/2,#/2}]&,{n,nt}]];
		{a,b}={(First[a]-Last[b]+1)/2,(First[b]-Last[a]+1)/2};
	];
	
	
	If[m==0,
		{x,y}=Transpose[ScalingFunctionData[low,NumIterations->its,Resolution->res,FilterStart->start]];
		data=Transpose[{x-start,y}];
	];
	If[m==1,
		{x,y}=Transpose[WaveletFunctionData[h,NumIterations->its,Resolution->res,FilterStart->start]];
		data=Transpose[{x-a,y}];
	];
		
	If[m>1,
		{adj,filter}=Which[EvenQ[m],{m/2,low},True,{(m-1)/2,g}];
		data=Map[Last,PacketFunctionData[h,adj,NumIterations->its,Resolution->res,FilterStart->start]];
		y=Join[data,ConstantArray[0,res*(Length[filter]-1)]];
		y=Sqrt[2]*Total[filter*Table[RotateRight[y,k*res],{k,0,Length[filter]-1}]];
		y=Take[y,{1,Length[y],2}];
		x=Most[Range[a,b,(b-a)/Length[y]]]-a;
		data=Transpose[{x,y}];
	];
	Return[data];
];

(* Error Messages *)

ScalingFunctionData::badfilter="A valid filter is a vector of numeric values - ScalingFunctionData failed."
ScalingFunctionData::iterations="Warning :: The number of iterations must be a nonnegative integer - resetting NumIterations to its default value."
ScalingFunctionData::resolution="Warning :: The value for the resolution should be a positive integer - resetting Resolution to its default value."
ScalingFunctionData::start="Warning :: The value for FilterStart must be an integer - resetting FilterStart to its default value."
WaveletFunctionData::badfilter="A valid filter is a vector or a pair of vectors of numeric values - WaveletFunctionData failed."
WaveletFunctionData::iterations="Warning :: The number of iterations must be a nonnegative integer - resetting NumIterations to its default value."
WaveletFunctionData::resolution="Warning :: The value for the resolution should be a positive integer - resetting Resolution to its default value."
WaveletFunctionData::start="Warning :: The value for FilterStart must be an integer - resetting FilterStart to its default value."
PacketFunctionData::badinteger="The value `1` is not valid for a packet function number.  This value must be a nonnegative integer - PacketFunctionData failed."
PacketFunctionData::badfilter="Either a vector of numeric values or a pair of vectors with numeric values is expected for the filter input - PacketFunctionData failed."
PacketFunctionData::start="Warning :: The value for FilterStart must be an integer - resetting FilterStart to its default value."
PacketFunctionData::iterations="Warning :: The number of iterations must be a nonnegative integer - resetting NumIterations to its default value."
PacketFunctionData::resolution="Warning :: The value for the resolution should be a positive integer - resetting Resolution to its default value."

End[]
EndPackage[]