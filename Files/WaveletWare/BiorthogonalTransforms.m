(* Mathematica package *)

BeginPackage["WaveletWare`BiorthogonalTransforms`",{"WaveletWare`CommonUsage`","WaveletWare`MiscFunctions`","WaveletWare`Filters`"}]

(* Exported symbols added here with SymbolName::usage *)

(* Biorthogonal transform functions *)

LeftBWT::usage="LeftBWT[a,{h,ht},opts] computes ~Wa, where ~W is a discrete biorthogonal wavelet transform matrix and a is a matrix or a list of three matrices each of the same dimensions.  The row dimension must be even."
RightBWT::usage="RightBWT[a,{h,ht},opts] computes aW^T, where W is a discrete biorthogonal wavelet transform matrix and a is a matrix  or a list of three matrices each of the same dimensions.  The column dimension must be even."
BWT::usage="BWT[data,{h,ht},opts] returns the biorthogonal wavelet transform of a vector, a matrix or a list of three matrices of equal dimension."
InverseBWTAux::usage="InverseBWTAux[v,h,type] is an auxiliary function used by InverseBWT, LeftInverseBWT and RightInverseBWT for computing the inverse biorthogonal wavelet transform."
LeftInverseBWT::usage="LeftInverseBWT[a,{h,ht},opts] computes Wa, where W is a discrete biorthogonal wavelet transform matrix and a is a matrix  or a list of three matrices each of the same dimensions.  The row dimension must be even."
RightInverseBWT::usage="RightInverseBWT[a,{h,ht},opts] computes a(~W^T), where ~W is a discrete biorthogonal wavelet transform matrix and a is a matrix  or a list of three matrices each of the same dimensions.  The column dimension must be even."
InverseBWT::usage="InverseBWT[data,{h,ht},opts] returns the inverse biorthogonal wavelet transform of a vector, matrix or a list of three matrices of equal dimensions."
LeftLWT::usage="LeftLWT[a,opts] computes ~Wa, where ~W is a discrete wavelet transform matrix constructed from the LeGall filter pair and a is a matrix  or a list of three matrices each of the same dimensions.  The row dimension must be even."
RightLWT::usage="RightLWT[a,opts] computes aW^T, where ~W is a discrete wavelet transform matrix constructed from the LeGall filter pair and a is a matrix  or a list of three matrices each of the same dimensions.  The column dimension must be even."
LWT::usage="LWT[data,opts] takes a vector, a matrix or a list of three matrices of equal dimensions and employs lifting to compute the biorthogonal wavelet transform using the LeGall filter pair."
LeftInverseLWT::usage="LeftInverseLWT[a,opts] computes W^Ta, where W is a discrete wavelet transform matrix constructed from the LeGall filter pair and a is a matrix  or a list of three matrices each of the same dimensions.  The row dimension must be even."
RightInverseLWT::usage="RightInverseLWT[a,opts] computes a(~W), where ~W is a discrete wavelet transform matrix constructed from the LeGall filter pair and a is a matrix  or a list of three matrices each of the same dimensions.  The column dimension must be even."
InverseLWT::usage="InverseLWT[data,opts] takes a vector, a matrix or a list of three matrices with equal dimensions and employs lifting to compute the inverse biorthogonal wavelet transform using the LeGall filter pair."

(* Option Values for packet transforms *)

(* Options for orthogonal transform functions *)

Options[LeftBWT]={Computation->Numerical,Boundary->False};
Options[RightBWT]=Options[LeftBWT];
Options[BWT]={NumIterations->1,Computation->Numerical,Boundary->False};
Options[LeftInverseBWT]={Computation->Numerical,Boundary->False};
Options[RightInverseBWT]=Options[LeftInverseBWT];
Options[InverseBWT]=Options[BWT];
Options[LeftLWT]={IntegerMap->False,Computation->Numerical};
Options[RightLWT]=Options[LeftLWT];
Options[LWT]={NumIterations->1,IntegerMap->False,Computation->Numerical};
Options[LeftInverseLWT]={IntegerMap->False,Computation->Numerical};
Options[RightInverseLWT]=Options[LeftInverseLWT];
Options[InverseLWT]=Options[LWT];

Begin["`Private`"]
(* Implementation of the package *)

(* Biorthogonal transform functions *)

LeftBWT[a_,{h_,ht_},OptionsPattern[]]:=Module[{datatype,maxits,n,nt,parity,l,lt,comp,boundary,p,gt,k,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[LeftBWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  2 = biorthogonal filter *)
	If[CheckFilter[{h,ht}]!=2,Message[LeftBWT::"badfilter"];Return[{}];];
	
	(* Check the parity of the filter pair *)
	{n,nt}=Map[Length,{h,ht}];
	parity=CheckParity[{n,nt}];
	If[CheckParity[{n,nt}]==0,Message[LeftBWT::"badparity"];Return[{}]];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,{h,ht}],Message[LeftBWT::"badsizes"];Return[{}]];
	
	(* Get offsets for computation *)
	{l,lt}={n, nt}/2 - parity/2;
	
	(* Read options and determine if they are valid*)
	{comp,boundary}=OptionValue[LeftBWT,{Computation,Boundary}];
	boundary=TrueFalse[boundary];
		
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	gt=h*Table[(-1)^k,{k,2-parity-l,l+1}];
	
	(* Create a function to perform the left BWT on a matrix *)
	f[x_]:=Module[{r},
		
		r[t_]:=Module[{lp,hp,s,m=Length[t]},
			If[boundary===False,
				s=t,
				s=Join[t,Reverse[Take[t,{3-parity,m+parity-2}]]];
			];
		
			{lp, hp} = MapThread[Partition[s,#1,2,#2]&,{{nt,n},{{lt+1,lt+2},{l+parity-1,l+2*(parity-1)}}}];
		
			If[boundary===False,
				Return[Join[lp.ht,hp.gt]],
				Return[Flatten[Map[Take[#,m/2]&,{lp.ht,hp.gt}]]];
			];
		];
	
		Return[Transpose[Map[r,Transpose[p*x]]]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

RightBWT[a_,{h_,ht_},OptionsPattern[]]:=Module[{datatype,maxits,n,nt,parity,l,lt,comp,boundary,p,gt,k,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[RightBWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  2 = biorthogonal filter *)
	If[CheckFilter[{h,ht}]!=2,Message[RightBWT::"badfilter"];Return[{}];];
	
	(* Check the parity of the filter pair *)
	{n,nt}=Map[Length,{h,ht}];
	parity=CheckParity[{n,nt}];
	If[CheckParity[{n,nt}]==0,Message[RightBWT::"badparity"];Return[{}]];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,{h,ht}],Message[RightBWT::"badsizes"];Return[{}]];
	
	(* Get offsets for computation *)
	{l,lt}={n, nt}/2 - parity/2;
	
	(* Read options and determine if they are valid*)
	{comp,boundary}=OptionValue[RightBWT,{Computation,Boundary}];
	boundary=TrueFalse[boundary];
			
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	gt=h*Table[(-1)^k,{k,2-parity-l,l+1}];
	
	(* Create a function to perform the left BWT on a matrix *)
	f[x_]:=Module[{r},
		
		r[t_]:=Module[{lp,hp,s,m=Length[t]},
			If[boundary===False,
				s=t,
				s=Join[t,Reverse[Take[t,{3-parity,m+parity-2}]]];
			];
		
			{lp, hp} = MapThread[Partition[s,#1,2,#2]&,{{nt,n},{{lt+1,lt+2},{l+parity-1,l+2*(parity-1)}}}];
		
			If[boundary===False,
				Return[Join[lp.ht,hp.gt]],
				Return[Flatten[Map[Take[#,m/2]&,{lp.ht,hp.gt}]]];
			];
		];
	
		Return[Map[r,p*x]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

BWT[a_,{h_,ht_},OptionsPattern[]]:=Module[{datatype,maxits,n,nt,parity,its,comp,boundary,l,lt,p,gt,k,data,f},
		
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{1,2,3},datatype],Message[BWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  2 = biorthogonal filter *)
	If[CheckFilter[{h,ht}]!=2,Message[BWT::"badfilter"];Return[{}];];
	
	(* Check the parity of the filter pair *)
	{n,nt}=Map[Length,{h,ht}];
	parity=CheckParity[{n,nt}];
	If[CheckParity[{n,nt}]==0,Message[BWT::"badparity"];Return[{}]];
	
	(* Read options and determine if they are valid*)
	{its,comp,boundary}=OptionValue[BWT,{NumIterations,Computation,Boundary}];
	boundary=TrueFalse[boundary];
	If[!CheckIterations[its,maxits],Message[BWT::"baditerations",maxits];its=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,its,{h,ht}],Message[BWT::"badsizes"];Return[{}]];
	
	(* Get offsets for computation *)
	{l,lt}={n, nt}/2 - parity/2;
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	gt=h*Table[(-1)^k,{k,2-parity-l,l+1}];
	
	(* Take care of the trivial case *)
	If[its==0,Return[p*a]];
	
	(* Do the vector case first *)
	If[datatype==1,
		data={p*a};
			
		f[x_]:=Module[{lp,hp,t,m=Length[First[x]]},
			If[boundary===False,
				t=First[x],
				t=Join[First[x],Reverse[Take[First[x],{3-parity,m+parity-2}]]];
			];
			
			{lp, hp} = MapThread[Partition[t,#1,2,#2]&,{{nt,n},{{lt+1,lt+2},{l+parity-1,l+2*(parity-1)}}}];
			
			If[boundary==False,
				Return[{lp.ht,hp.gt,Drop[x,1]}],
				Return[{Take[lp.ht,m/2],Take[hp.gt,m/2],Drop[x,1]}];
			];
		];
		
		Return[Flatten[Nest[f,data,its]]];
	];
	
	(* Here is a function for performing iterations of the wavelet transform on a matrix *)
	f[x_]:=Module[{r},
		r[t_]:=Join[WaveletMatrixToList[RightBWT[LeftBWT[First[t],{h,ht},Computation->comp,Boundary->boundary],{h,ht},Computation->comp,Boundary->boundary]],Drop[t,1]];
		Return[WaveletListToMatrix[Nest[r,x,its],NumIterations->its]];
	];
	
	(* Now do the matrix and list of matrices cases *)
	If[datatype==2,
		Return[f[{p*a}]];
	];
	
	If[datatype==3,
		Return[Map[f[{#}]&,p*a]];
	];
];

InverseBWTAux[v_,f_,type_]:=Module[{n=Length[f],l,idx,oddf,evenf,evenl,oddl,evenv,oddv},
	(* v is an input vector, f is a filter, and type is 0 if the filter is lowpass and 1 if the filter is highpass. *)
   	l=n/2-If[OddQ[n],1/2,1];
   	
   	idx=If[(type==1&&OddQ[n]),If[OddQ[l],{1,2},{2,1}],If[OddQ[l],{2,1},{1,2}]];
   	{evenf,oddf}=Map[Reverse,Map[Take[f,{#,n,2}]&,idx]];

   	{evenl,oddl}=Map[Length,{evenf,oddf}];   	
   	
   	If[OddQ[n],
   		{evenv,oddv}=MapThread[Partition[v,#1,1,{#1,#1}/2+#2]&,{{evenl,oddl},{1/2,0}+type/2}],
   		{evenv,oddv}=MapThread[Partition[v,#1,1,{#1,#1}/2+#2]&,{{evenl,oddl},{If[EvenQ[evenl],1,1/2],If[EvenQ[oddl],0,1/2]}}];
   	];
   		
   	Return[Riffle[evenv.evenf,oddv.oddf]];
];

LeftInverseBWT[a_,{h_,ht_},OptionsPattern[]]:=Module[{datatype,maxits,n,nt,parity,l,lt,comp,boundary,p,g,k,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[LeftInverseBWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  2 = biorthogonal filter *)
	If[CheckFilter[{h,ht}]!=2,Message[LeftInverseBWT::"badfilter"];Return[{}];];
	
	(* Check the parity of the filter pair *)
	{n,nt}=Map[Length,{h,ht}];
	parity=CheckParity[{n,nt}];
	If[CheckParity[{n,nt}]==0,Message[LeftInverseBWT::"badparity"];Return[{}]];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,{h,ht}],Message[LeftInverseBWT::"badsizes"];Return[{}]];
	
	(* Get offsets for computation *)
	{l,lt}={n, nt}/2 - parity/2;
	
	(* Read options and determine if they are valid*)
	{comp,boundary}=OptionValue[LeftInverseBWT,{Computation,Boundary}];
	boundary=TrueFalse[boundary];
		
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	g=ht*Table[(-1)^k,{k,2-parity-lt,lt+1}];
	
	(* Create a function to perform the left BWT on a matrix *)
	f[x_]:=Module[{r},
		r[t_]:=Module[{y,u,m=Length[t]},
			If[boundary===False,
				u=Partition[t,m/2],
				If[parity==1,
					u=MapThread[Join[#1,Reverse[Take[#1,#2]]]&,{Partition[t,m/2],{{2,m/2},{1,m/2-1}}}], 
  					u=MapThread[Join[#1,Reverse[#1*#2]]&,{Partition[t,m/2],{1, -1}}];
  				];
			];
			y = Total[MapThread[InverseBWTAux[#1, #2, #3] &, {u, {h, g}, {0, 1}}]];
			If[boundary===False,
				Return[y],
				Return[Take[y,m]];
			];
		];
		
		Return[Transpose[Map[r,Transpose[p*x]]]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

RightInverseBWT[a_,{h_,ht_},OptionsPattern[]]:=Module[{datatype,maxits,n,nt,parity,l,lt,comp,boundary,p,g,k,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[RightInverseBWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  2 = biorthogonal filter *)
	If[CheckFilter[{h,ht}]!=2,Message[RightInverseBWT::"badfilter"];Return[{}];];
	
	(* Check the parity of the filter pair *)
	{n,nt}=Map[Length,{h,ht}];
	parity=CheckParity[{n,nt}];
	If[CheckParity[{n,nt}]==0,Message[RightInverseBWT::"badparity"];Return[{}]];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,{h,ht}],Message[RightInverseBWT::"badsizes"];Return[{}]];
	
	(* Get offsets for computation *)
	{l,lt}={n, nt}/2 - parity/2;
	
	(* Read options and determine if they are valid*)
	{comp,boundary}=OptionValue[RightInverseBWT,{Computation,Boundary}];
	boundary=TrueFalse[boundary];
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	g=ht*Table[(-1)^k,{k,2-parity-lt,lt+1}];
	
	(* Create a function to perform the right BWT on a matrix *)
	f[x_]:=Module[{r},
		r[t_]:=Module[{y,u,m=Length[t]},
			If[boundary===False,
				u=Partition[t,m/2],
				If[parity==1,
					u=MapThread[Join[#1,Reverse[Take[#1,#2]]]&,{Partition[t,m/2],{{2,m/2},{1,m/2-1}}}], 
  					u=MapThread[Join[#1,Reverse[#1*#2]]&,{Partition[t,m/2],{1, -1}}];
  				];
			];
			y = Total[MapThread[InverseBWTAux[#1, #2, #3] &, {u, {h, g}, {0, 1}}]];
			If[boundary===False,
				Return[y],
				Return[Take[y,m]];
			];
		];
		
		Return[Map[r,p*x]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

InverseBWT[a_,{h_,ht_},OptionsPattern[]]:=Module[{datatype,maxits,n,nt,parity,its,comp,boundary,l,lt,p,g,k,data,f},
		
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{1,2,3},datatype],Message[InverseBWT::"badinput"];Return[{}]];
	
	(* Check the validity of the filter.  2 = biorthogonal filter *)
	If[CheckFilter[{h,ht}]!=2,Message[InverseBWT::"badfilter"];Return[{}];];
	
	(* Check the parity of the filter pair *)
	{n,nt}=Map[Length,{h,ht}];
	parity=CheckParity[{n,nt}];
	If[CheckParity[{n,nt}]==0,Message[InverseBWT::"badparity"];Return[{}]];
	
	(* Read options and determine if they are valid*)
	{its,comp,boundary}=OptionValue[InverseBWT,{NumIterations,Computation,Boundary}];
	boundary=TrueFalse[boundary];
	If[!CheckIterations[its,maxits],Message[InverseBWT::"baditerations",maxits];its=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,its,{h,ht}],Message[InverseBWT::"badsizes"];Return[{}]];
	
	(* Get offsets for computation *)
	{l,lt}={n, nt}/2 - parity/2;
	
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];	
	
	(* Create the highpass filter *)
	g=ht*Table[(-1)^k,{k,2-parity-lt,lt+1}];
	
	(* Take care of the trivial case *)
	If[its==0,Return[p*a]];
	
	(* Do the vector case first *)
	If[datatype==1,
		data=WaveletVectorToList[p*a,NumIterations->its];
		
		f[x_]:=Module[{y,u,m=Length[First[x]]},
			If[boundary===False,
				u=Take[x,2],
				If[parity==1,
					u=MapThread[Join[#1,Reverse[Take[#1,#2]]]&,{Take[x,2],{{2,m},{1,m-1}}}],
					u=MapThread[Join[#1,Reverse[Take[#1*#2,{1,m}]]]&,{Take[x,2],{1,-1}}];
				];
			];
			
			y=Total[MapThread[InverseBWTAux[#1,#2,#3]&,{u,{h,g},{0,1}}]];
			
			If[boundary===False,
				Return[Prepend[Drop[x,2],y]],
				Return[Prepend[Drop[x,2],Take[y,2*m]]];
			];
		];
		
		Return[First[Nest[f,data,its]]];
	];
	
	(* Create a function to perform the InverseBWT on a matrix *)
	
	f[x_]:=Module[{d,r},
		d=WaveletMatrixToList[x,NumIterations->its];
		
		r[t_]:=Module[{y},
			y=RightInverseBWT[LeftInverseBWT[WaveletListToMatrix[Take[t,2]],{h,ht},Boundary->boundary,Computation->comp],{h,ht},Boundary->boundary,Computation->comp];
			Return[Prepend[Drop[t,2],y]];
		];
		
		Return[First[Nest[r,d,its]]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

LeftLWT[a_,OptionsPattern[]]:=Module[{datatype,maxits,intmap,comp,p,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[Left=LeftLWT::"badinput"];Return[{}]];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,LeGall[]],Message[LeftLWT::"badsizes"];Return[{}]];

	(* Read options and determine if they are valid*)
	{intmap,comp}=OptionValue[LeftLWT,{IntegerMap,Computation}];
	intmap=TrueFalse[intmap];
				
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];
	
	(* Override the type of computation if IntegerMap is set to True. *)
	If[intmap===True,p=1];
	
	(* If IntegerMap is true, then make sure the input elements are integers *)
	If[intmap===True && !AllTrue[Flatten[a],IntegerQ],Message[LeftLWT::"notintegerdata"];intmap=False;];
	
	(* A function to apply the left LWT to a matrix *)
	f[x_]:=Module[{r},
		r[t_]:=Module[{e,o,u,d,s},
			{e,o}=Map[Take[p*t,{#,Length[t],2}]&,{2,1}];
   			o=Append[o,Last[o]];
   			u=Map[Total,Partition[o,2,1]]/2;
   			d=e-If[intmap===True,Floor[u],u];
   			d=Prepend[d,First[d]];
   			u=Map[Total,Partition[d,2,1]]/4;
   			s=Drop[o,-1]+If[intmap===True,Floor[u+1/2],u];
   			Return[Join[s,Drop[d,1]]];
		];
	
		Return[Transpose[Map[r,Transpose[x]]]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

RightLWT[a_,OptionsPattern[]]:=Module[{datatype,maxits,intmap,comp,p,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[Left=RightLWT::"badinput"];Return[{}]];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,LeGall[]],Message[RightLWT::"badsizes"];Return[{}]];

	(* Read options and determine if they are valid*)
	{intmap,comp}=OptionValue[RightLWT,{IntegerMap,Computation}];
	intmap=TrueFalse[intmap];
				
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];
	
	(* Override the type of computation if IntegerMap is set to True. *)
	If[intmap===True,p=1];
	
	(* If IntegerMap is true, then make sure the input elements are integers *)
	If[intmap===True && !AllTrue[Flatten[a],IntegerQ],Message[RightLWT::"notintegerdata"];intmap=False;];
	
	f[x_]:=Module[{r},
		r[t_]:=Module[{e,o,u,d,s},
			{e,o}=Map[Take[p*t,{#,Length[t],2}]&,{2,1}];
   			o=Append[o,Last[o]];
   			u=Map[Total,Partition[o,2,1]]/2;
   			d=e-If[intmap===True,Floor[u],u];
   			d=Prepend[d,First[d]];
   			u=Map[Total,Partition[d,2,1]]/4;
   			s=Drop[o,-1]+If[intmap===True,Floor[u+1/2],u];
   			Return[Join[s,Drop[d,1]]];
		];
	
		Return[Map[r,x]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

LWT[a_,OptionsPattern[]]:=Module[{datatype,maxits,its,intmap,comp,p,data,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{1,2,3},datatype],Message[LWT::"badinput"];Return[{}]];
	
	(* Read options and determine if they are valid*)
	{its,intmap,comp}=OptionValue[LWT,{NumIterations,IntegerMap,Computation}];
	intmap=TrueFalse[intmap];
		
	If[!CheckIterations[its,maxits],Message[LWT::"baditerations",maxits];its=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,its,LeGall[]],Message[LWT::"badsizes"];Return[{}]];
			
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];
	
	(* Override the type of computation if IntegerMap is set to True. *)
	If[intmap==True,p=1];
	
	(* If IntegerMap is true, then make sure the input elements are integers *)
	If[(intmap && !AllTrue[Flatten[a],IntegerQ]),Message[LWT::"notintegerdata"];intmap=False;];
	
	If[its==0,Return[p*a]];
	
	(* Do the vector case first *)
	If[datatype==1,
		data={p*a};
		
		f[x_]:=Module[{e,o,m=Length[First[x]],t,d,s},
			{e,o}=Map[Take[p*First[x],{#,m,2}]&,{2,1}];
   			o=Append[o,Last[o]];
   			t=Map[Total,Partition[o,2,1]]/2;
   			d=e-If[intmap===True,Floor[t],t];
   			d=Prepend[d,First[d]];
   			t=Map[Total,Partition[d,2,1]]/4;
   			s=Drop[o,-1]+If[intmap===True,Floor[t+1/2],t];
   			Return[{s,Drop[d,1],Drop[x,1]}];
		];
		
		Return[Flatten[Nest[f,data,its]]];
		If[its==1,
			Return[f[data]],
			Return[Flatten[Nest[f,data,its]]];
		];
	];
	
	(* Here is a function for performing iterations of the LeGall wavelet transform on a matrix *)
	f[x_]:=Module[{d,r},
		d={p*x};
		r[t_]:=Join[WaveletMatrixToList[RightLWT[LeftLWT[First[t],IntegerMap->intmap,Computation->comp],IntegerMap->intmap,Computation->comp]],Drop[t,1]];
		Return[WaveletListToMatrix[Nest[r,d,its],NumIterations->its]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

LeftInverseLWT[a_,OptionsPattern[]]:=Module[{datatype,maxits,intmap,comp,p,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[Left=LeftInverseLWT::"badinput"];Return[{}]];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,LeGall[]],Message[LeftInverseLWT::"badsizes"];Return[{}]];

	(* Read options and determine if they are valid*)
	{intmap,comp}=OptionValue[LeftInverseLWT,{IntegerMap,Computation}];
	intmap=TrueFalse[intmap];
			
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];
	
	(* Override the type of computation if IntegerMap is set to True. *)
	If[intmap===True,p=1];
	
	(* If IntegerMap is true, then make sure the input elements are integers *)
	If[intmap===True && !AllTrue[Flatten[a],IntegerQ],Message[LeftInverseLWT::"notintegerdata"];intmap=False;];
	
	(* A function to apply the left LWT to a matrix *)
	f[x_]:=Module[{r},
		r[t_]:=Module[{e,o,u,d,s},
			{s,d}=Partition[p*t,Length[t]/2];
    		d=Prepend[d,First[d]];
    		u=Map[Total,Partition[d,2,1]]/4;
    		o=s-If[intmap===True,Floor[u+1/2],u];
    		o=Append[o,Last[o]];
    		u=Map[Total,Partition[o,2,1]]/2;
    		e=Drop[d,1]+If[intmap===True,Floor[u],u];
    		Return[Flatten[Transpose[{Drop[o,-1],e}]]];
		];
		Return[Transpose[Map[r,Transpose[x]]]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

RightInverseLWT[a_,OptionsPattern[]]:=Module[{datatype,maxits,intmap,comp,p,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{2,3},datatype],Message[Left=RightInverseLWT::"badinput"];Return[{}]];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,1,LeGall[]],Message[RightInverseLWT::"badsizes"];Return[{}]];

	(* Read options and determine if they are valid*)
	{intmap,comp}=OptionValue[RightInverseLWT,{IntegerMap,Computation}];
	intmap=TrueFalse[intmap];
			
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];
	
	(* Override the type of computation if IntegerMap is set to True. *)
	If[intmap===True,p=1];
	
	(* If IntegerMap is true, then make sure the input elements are integers *)
	If[intmap===True && !AllTrue[Flatten[a],IntegerQ],Message[RightInverseLWT::"notintegerdata"];intmap=False;];
	
	(* A function to apply the right LWT to a matrix *)
	f[x_]:=Module[{r},
		r[t_]:=Module[{e,o,u,d,s},
			{s,d}=Partition[p*t,Length[t]/2];
    		d=Prepend[d,First[d]];
    		u=Map[Total,Partition[d,2,1]]/4;
    		o=s-If[intmap===True,Floor[u+1/2],u];
    		o=Append[o,Last[o]];
    		u=Map[Total,Partition[o,2,1]]/2;
    		e=Drop[d,1]+If[intmap===True,Floor[u],u];
    		Return[Flatten[Transpose[{Drop[o,-1],e}]]];
		];
		Return[Map[r,x]];
	];
	
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

InverseLWT[a_,OptionsPattern[]]:=Module[{datatype,maxits,its,intmap,comp,p,data,f},
	
	(* First check data integrity *)
	{datatype,maxits}=CheckData[a];
	If[!MemberQ[{1,2,3},datatype],Message[InverseLWT::"badinput"];Return[{}]];
	
	(* Read options and determine if they are valid*)
	{its,intmap,comp}=OptionValue[InverseLWT,{NumIterations,IntegerMap,Computation}];
	intmap=TrueFalse[intmap];
	If[!CheckIterations[its,maxits],Message[InverseLWT::"baditerations",maxits];its=0;];
	
	(* Check that the filter size, iterations and the length of the data are compatible for doing the transform *)
	If[!CheckDataFilterSize[datatype,a,its,LeGall[]],Message[InverseLWT::"badsizes"];Return[{}]];
			
	(* Determine if the computation is to be numerical or symbolic *)
	p=Which[comp===Symbolic,1,True,1.];
	
	(* Override the type of computation if IntegerMap is set to True. *)
	If[intmap==True,p=1];
	
	(* If IntegerMap is true, then make sure the input elements are integers *)
	If[intmap===True && !AllTrue[Flatten[a],IntegerQ],Message[InverseLWT::"notintegerdata"];intmap=False;];
	
	If[its==0,Return[p*a]];
	
	(* Do the vector case first *)

	If[datatype==1,
		data=WaveletVectorToList[p*a,NumIterations->its];
		
		f[x_]:=Module[{e,o,t,d,s},
			{s,d}=Take[p*x,2];
    		d=Prepend[d,First[d]];
    		t=Map[Total,Partition[d,2,1]]/4;
    		o=s-If[intmap===True,Floor[t+1/2],t];
    		o=Append[o,Last[o]];
    		t=Map[Total,Partition[o,2,1]]/2;
    		e=Drop[d,1]+If[intmap===True,Floor[t],t];
    		Return[Prepend[Drop[x,2],Flatten[Transpose[{Drop[o,-1],e}]]]];
		];
		
		Return[Flatten[Nest[f,data,its]]];
	];
	
	(* Here is a function for performing iterations of the inverse LeGall wavelet transform on a matrix *)
	f[x_]:=Module[{d,r},
		d=WaveletMatrixToList[x,NumIterations->its];
		r[t_]:=Module[{y},
			y=LeftInverseLWT[RightInverseLWT[WaveletListToMatrix[Take[t,2]],IntegerMap->intmap,Computation->comp],IntegerMap->intmap,Computation->comp];
			Return[Prepend[Drop[t,2],y]];
		];
		Return[First[Nest[r,d,its]]];
	];
		
	If[datatype==2,Return[f[a]],Return[Map[f,a]]];
];

(* Error messages *)

LeftBWT::badinput="The input for LeftBWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - LeftBWT failed."
LeftBWT::badfilter="The input filter must be a biorthogonal lowpass filter pair whose lengths are both even or both odd - LeftBWT failed."
LeftBWT::badparity="The lengths of the filters must either be both even or both odd - LeftBWT failed."
LeftBWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - LeftBWT failed."
RightBWT::badinput="The input for LeftBWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - RightBWT failed."
RightBWT::badfilter="The input filter must be a biorthogonal lowpass filter pair whose lengths are both even or both odd - RightBWT failed."
RightBWT::badparity="The lengths of the filters must either be both even or both odd - RightBWT failed."
RightBWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - RightBWT failed."
BWT::badinput="The input for BWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - BWT failed."
BWT::badfilter="The input filter must be a biorthogonal lowpass filter pair whose lengths are both even or both odd - BWT failed."
BWT::badparity="The lengths of the filters must either be both even or both odd - BWT failed."
BWT::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
BWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - BWT failed."
LeftInverseBWT::badinput="The input for LeftInverseBWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - LeftInverseBWT failed."
LeftInverseBWT::badfilter="The input filter must be a biorthogonal lowpass filter pair whose lengths are both even or both odd - LeftInverseBWT failed."
LeftInverseBWT::badparity="The lengths of the filters must either be both even or both odd - LeftInverseBWT failed."
LeftInverseBWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - LeftBWT failed."
RightInverseBWT::badinput="The input for LeftInverseBWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - RightInverseBWT failed."
RightInverseBWT::badfilter="The input filter must be a biorthogonal lowpass filter pair whose lengths are both even or both odd - RightInverseBWT failed."
RightInverseBWT::badparity="The lengths of the filters must either be both even or both odd - RightInverseBWT failed."
RightInverseBWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - RightBWT failed."
InverseBWT::badinput="The input for InverseBWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - InverseBWT failed."
InverseBWT::badfilter="The input filter must be a biorthogonal lowpass filter pair whose lengths are both even or both odd - InverseBWT failed."
InverseBWT::badparity="The lengths of the filters must either be both even or both odd - InverseBWT failed."
InverseBWT::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
InverseBWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - InverseBWT failed."
LeftLWT::badinput="The input for LeftLWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - LeftLWT failed."
LeftLWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - LeftLWT failed."
LeftLWT::notintegerdata="Warning :: The option IntegerMap has been set to True but the input data are not integer-valued.  Resetting IntegerMap to False."
RightLWT::badinput="The input for RightLWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - RightLWT failed."
RightLWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - RightLWT failed."
RightLWT::notintegerdata="Warning :: The option IntegerMap has been set to True but the input data are not integer-valued.  Resetting IntegerMap to False."
LWT::badinput="The input for LWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - LWT failed."
LWT::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
LWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - LWT failed."
LWT::notintegerdata="Warning :: The option IntegerMap has been set to True but the input data are not integer-valued.  Resetting IntegerMap to False."
LeftInverseLWT::badinput="The input for LeftInverseLWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - LeftInverseLWT failed."
LeftInverseLWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - LeftInverseLWT failed."
LeftInverseLWT::notintegerdata="Warning :: The option IntegerMap has been set to True but the input data are not integer-valued.  Resetting IntegerMap to False."
RightInverseLWT::badinput="The input for RightInverseLWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - RightInverseLWT failed."
RightInverseLWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - RightInverseLWT failed."
RightInverseLWT::notintegerdata="Warning :: The option IntegerMap has been set to True but the input data are not integer-valued.  Resetting IntegerMap to False."
InverseLWT::badinput="The input for InverseLWT is numeric and must be a vector, matrix or a list of three matrices with equal dimensions - InverseLWT failed."
InverseLWT::baditerations="Warning :: The value for NumIterations, due to the length/dimensions of the input, must be a nonnegative integer less than or equal to `1` - resetting the value to 0."
InverseLWT::badsizes="Either the length/dimensions of the input data are not of the appropriate size or the length of the filter(s) are too large for the number of iterations requested  - InverseLWT failed."
InverseLWT::notintegerdata="Warning :: The option IntegerMap has been set to True but the input data are not integer-valued.  Resetting IntegerMap to False."

End[]
EndPackage[]