(* Mathematica package *)

BeginPackage["WaveletWare`AudioFunctions`",{"WaveletWare`CommonUsage`","WaveletWare`MiscFunctions`"}]

(* Exported symbols added here with SymbolName::usage *)

(* Audio functions *)

AudioFormats::usage="AudioFormats is a WaveletWare package variable that returns a list of acceptable Mathematica audio file formats."
AudioInfo::usage="AudioInfo[opts] returns the absolute path name for each audio file in the WaveletWare package.  It can also supply information about the files."
AudioRead::usage="AudioRead[file,opts] is a function that reads a Mathematica audio file returning the raw data and the sample rate."
WaveletVectorPlay::usage="WaveletVectorPlay[v,opts] or WaveletVectorPlay[v,tree,opts] plays all or parts of the wavelet (packet) transform of a vector."

(* Options for audio functions *)

Options[AudioInfo]={DisplayInfo->False};
Options[AudioRead]={DisplayInfo->False,PowersOfTwo->0};
Options[WaveletVectorPlay]=Join[Options[ListPlay],{NumIterations->1,Iteration->All}];
SetOptions[WaveletVectorPlay,SampleRate->11025];

Begin["`Private`"]
(* Implementation of the package *)

(* Audio Formats *)

AudioFormats={"AIFF","AU","FLAC","MP3","Ogg","SND","WAV","Wave64","MIDI"};

(* Audio Functions *)

AudioInfo[OptionsPattern[]]:=Module[{HomeDirectory=Directory[],AudioDirectory,files,filenames,channels,rate,data,size,type,its},
	
	AudioDirectory=PackageDirectory[DataType->Sounds];
	SetDirectory[AudioDirectory];
	files=Map[StringJoin[AudioDirectory,#]&,Union[FileNames[],FileNames["*",{"*"},Infinity]]];
	files=Select[files,MemberQ[AudioFormats,FileFormat[#]]&];
	SetDirectory[HomeDirectory];
	
	(* Display audio info if requested. *)
	If[TrueFalse[OptionValue[AudioInfo,DisplayInfo]],
		{channels,data,rate}=Transpose[Map[Import[#,{{"AudioChannels","Data","SampleRate"}}]&,files]];
		size=Map[Length[Which[VectorQ[#],#,True,First[#]]]&,data];
		its=Map[MaxIts,data];
		type=Map[FileFormat,files];
		filenames=Map[FileNameDrop[#,FileNameDepth[AudioDirectory]]&,files];
		Print["Information for Audio Files in the WaveletWare Package\n",TableForm[Transpose[{Range[Length[filenames]],filenames,size,type,channels,rate,its}],TableHeadings->{None,{"File Number","Audio File","Size","Type","Channels","Sample Rate","Max Iterations"}},TableAlignments->{Center,Center},TableSpacing->{5,5}]];
	];
	Return[files];
];

AudioRead[f_,OptionsPattern[]]:=Module[{rate,data,channels,size,ptwo,xtra,maxits},
	{rate,data,channels}=Import[f,{{"SampleRate","Data","AudioChannels"}}];
	size=Length[Which[channels==1,data,True,First[data]]];
	ptwo=OptionValue[AudioRead,PowersOfTwo];
	If[!(IntegerQ[ptwo] && ptwo>=0),Message[AudioRead::"badpower"];ptwo=0];
	xtra=Mod[size,2^ptwo];
	data=Which[channels==1,Drop[data,xtra],True,Map[Drop[#,xtra]&,data]];
	size=Length[Which[channels==1,data,True,First[data]]];
	maxits=Last[CheckData[data]];
	
	If[TrueFalse[OptionValue[AudioRead,DisplayInfo]],
		Print["The number of samples returned is ",size," and this length is divisible by at least ",2^maxits,".  The sample rate is ",rate," and the number of channels returned is ",channels,"."];
	];
	Return[{data,rate}];
];

WaveletVectorPlay[v_,t___,OptionsPattern[]]:=Module[{datatype,maxits,its,tree,dsp,e,rng,p,d,rate,iter,n,f,clips,lbls},
	{datatype,maxits}=CheckData[v];
	(* Stereo data *)
	maxits=If[datatype==2&&Length[v]==2,Last[CheckData[Last[v]]],maxits];
	If[!MemberQ[{1,2},datatype],Message[WaveletVectorPlay::"badinput"];Return[{}]];
	
	(*Check the iterations and get the tree-if no tree is given,generate a WaveletTree*)
	its=OptionValue[WaveletVectorPlay,NumIterations];
	{its,tree}=If[t==="",{its,WaveletTree[its]},{CheckTree[t,1], t}];
	If[its==False,Message[WaveletVectorPlay::"badtree"];Return[{}]];
	If[!CheckIterations[its,maxits],Message[WaveletVectorPlay::"baditerations",maxits];its=0];
	
	(* Grab the other option values *)
	{dsp,e,rng,p,d,rate,iter}=OptionValue[WaveletVectorPlay,{DisplayFunction,Epilog,PlayRange,Prolog,SampleDepth,SampleRate,Iteration}];
	(* Handle the singular case where its = 0 *)
	If[its==0||First[tree]==1,Return[Column[{TableForm[{Text[Style["Original Data",FontFamily->"Helvetica"]],ListPlay[v,DisplayFunction->dsp,Epilog->e,PlayRange->rng,Prolog->p,SampleDepth->d,SampleRate->rate]},TableAlignments->Center]},Background->Lighter[Gray],Frame->All]]];
	If[!(iter===All||AllTrue[iter, CheckIterationValue[#,its,tree,1]&]),Message[WaveletVectorPlay::"baditeration"];Return[{}]];
	iter=If[iter===All,WaveletRegionList[tree],iter];

	(* Grab the different clips and new sample rates *)
	n=If[datatype==1,Length[v],Length[First[v]]];
	f[{i_,j_}]:=If[datatype==1,Partition[v,n/2^i][[j]],Map[Partition[#,n/2^i][[j]]&,v]];
	clips=Map[f,iter];
	rate=Round[rate*(If[datatype==1,Map[Length,clips],Map[Length[Last[#]]&,clips]]/n)];
	clips=MapThread[ListPlay[#1,DisplayFunction->dsp,Epilog->e,PlayRange->rng,Prolog->p,SampleDepth->d,SampleRate->#2]&,{clips,rate}];

	(* Create labels *)
	f[{i_,j_}]:=Text[Style["Iteration "<>ToString[i]<>", Region "<>ToString[j],FontFamily->"Helvetica"]];
	lbls=Map[f,iter];

	Return[Column[Map[TableForm[#,TableAlignments->Center]&,Transpose[{lbls,clips}]],Background->Lighter[Gray],Spacings->{0,1},Frame ->All]];
];

(*
WaveletVectorPlay[v_,OptionsPattern[]]:=Module[{dsp,e,rng,p,d,its,iteration,rate,maxits,srate,clip,lbls={},snd={},f},
	
	(* Check the integrity of the input vector *)
	If[!ArrayQ[v]||!AllTrue[Flatten[v],NumericQ],Message[WaveletVectorPlay::"badinput",v];Return[{}]];
	
	(* Read options and determine if they are valid*)
	{dsp,e,rng,p,d,rate,its,iteration}=Map[OptionValue[WaveletVectorPlay,#]&,Keys[Options[WaveletVectorPlay]]];
	
	(* Check the integrity of NumIterations *)
	maxits=If[VectorQ[v],MaxIts[v],MaxIts[First[v]]];
	If[its>maxits || !IntegerQ[its] || its<0,Message[WaveletVectorPlay::"baditerations",maxits];Return[{}]];
	
	(* Convert transform to a list *)
	clip=If[VectorQ[v],WaveletVectorToList[v,NumIterations->its],Map[WaveletVectorToList[#,NumIterations->its]&,v]];
	
	(* Generate different sample rates *)
	srate=Round[rate/2^Range[its]];
	srate=Reverse[Append[srate,Last[srate]]];
	
	Which[
		(* Check and see if we are to plot the low/highpass region of the final iteration *)
		VectorQ[iteration]&&Length[iteration]==2&&First[iteration]==its&&MemberQ[{1,2},Last[iteration]],
			clip=If[VectorQ[v],clip[[Last[iteration]]],Map[clip[[#,Last[iteration]]]&,Range[Length[clip]]]];
			lbls=If[Last[iteration]==1,"Lowpass","Highpass"]<>", Iteration "<>ToString[First[iteration]];
			snd=ListPlay[clip,DisplayFunction->dsp,Epilog->e,PlayRange->rng,Prolog->p,SampleDepth->d,SampleRate->srate[[Last[iteration]]]],
		(* Handle the special case where both parts of the final iteration are requested *)
		IntegerQ[iteration]&&iteration==its,
			clip=If[VectorQ[v],Flatten[Take[clip,2]],Map[Flatten[Take[#,2]]&,clip]];
			lbls="Iteration "<>ToString[iteration];
			snd=ListPlay[clip,DisplayFunction->dsp,Epilog->e,PlayRange->rng,Prolog->p,SampleDepth->d,SampleRate->First[srate]],
		(* Get an iteration besides the last one *)
		IntegerQ[iteration]&&MemberQ[Range[0,its-1],iteration],
			clip=If[VectorQ[v],Reverse[Rest[clip]][[iteration]],Map[Reverse[Rest[#]][[iteration]]&,clip]];
			lbls="Iteration "<>ToString[iteration];
			snd=ListPlay[clip,DisplayFunction->dsp,Epilog->e,PlayRange->rng,Prolog->p,SampleDepth->d,SampleRate->Reverse[Rest[srate]][[iteration]]],
		(* All iterations are requested *)
		iteration==All,
			f[m_]:=Map[Flatten[Take[#1,{m}]]&,clip];
			clip=If[VectorQ[v],clip,Map[f,Range[its+1]]];
			lbls=Prepend[Map[("Highpass "<>ToString[#])&,Range[its]],"Lowpass"];
			snd=MapThread[ListPlay[#1,DisplayFunction->dsp,Epilog->e,PlayRange->rng,Prolog->p,SampleDepth->d,SampleRate->#2]&,{clip,srate}];
	];
	
	If[lbls=={}&&snd=={},Message[WaveletVectorPlay::"baditeration"];Return[{}]];
	
	Return[TableForm[{lbls,snd}]];
];
*)

(* Error Messages *)

AudioInfo::packagenotfound="WaveletWare package not found in any $Path folders.  Check installation instructions."

AudioRead::badpower="Warning :: The power of 2 must be a nonnegative integer - resetting the power to 0."

WaveletVectorPlay::badinput="The input must be a vector consisting of numeric values - WaveletVectorPlay failed."
WaveletVectorPlay::baditeratons="The value of NumIterations must be a nonnegative integer no larger than `1` - WaveletVectorPlay failed."
WaveletVectorPlay::baditeration="The value of Iteration should be All (default), a nonnegative integer 0,...,its, or a list that is either {its,1} (lowpass) or {its,2} (highpass)."


End[]
EndPackage[]