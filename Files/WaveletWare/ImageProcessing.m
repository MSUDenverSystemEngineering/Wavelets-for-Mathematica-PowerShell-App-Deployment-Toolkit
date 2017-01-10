(* Mathematica package *)

BeginPackage["WaveletWare`ImageProcessing`",{"WaveletWare`CommonUsage`","WaveletWare`MiscFunctions`","WaveletWare`OrthogonalTransforms`"}]

(* Exported symbols added here with SymbolName::usage *)

(* Image processing functions *)

ImageFormats::usage="ImageFormats is a WaveletWare package variable that is a list of the Mathematica image raster formats."
ImageInfo::usage="ImageInfo[opts] returns the absolute path name for each image in the WaveletWare package.  It can also supply thumbnail sketches and data about the package images."
ImageRead::usage="ImageRead[file,opts] is a function that reads a Mathematica raster image from disk or a url and converts it to a matrix (grayscale image) or three matrices (color image)."
ImageEnlarge::usage="ImageEnlarge[a,scale] uses the inverse Haar wavelet transform to enlarge the input image a by the given scale factor."
RGBToYCbCr::usage="RGBToYCbCr[{r,g,b},opts] is a function that converts either an r,g,b triple or a list of r,g,b matrices to YCbCr space."
YCbCrToRGB::usage="YCbCrToRGB[{y,cb,cr},opts] is a function that converts either an y,cb,cr triple or a list of y,cb,cr matrices to RGB space."
RGBToHSI::usage="RGBToHSI[{r,g,b},opts] is a function that converts either an r,g,b triple or a list of r,g,b matrices to HSI space."
HSIToRGB::usage="HSIToRGB[{h,s,i},opts] is a function that converts either an h,s,i triple or a list of h,s,i matrices to RGB space."
DCT::usage="DCT[a] computes the discrete cosine transform (II) of an input vector or matrix."
InverseDCT::usage="InverseDCT[a] computes the inverse discrete cosine transform (II) of an input vector or matrix."
CE::usage="CE[a,opts] computes the cumulative energy vector for an input vector or matrix."
PercentCE::usage="PercentCE[v,pct] takes a cumulative energy vector v and a value 0 <= pct <=1 and returns an index k that represents the number of values in v that are less than or equal to pct."
nCE::usage="nCE[v,pct] takes a cumulative energy vector v and a value 0 <= pct <=1 and returns an index k that represents the number of values in v that are less than or equal to pct."
QuantizeCE::usage="QuantizeCE[a,k] takes a vector or matrix a and a positive integer k and returns a vector or matrix where the largest k (absolute) values plus any redundancies are preserved and the remaining values are converted to 0."
Comp::usage="Comp[a,k] takes a vector or matrix a and a positive integer k and returns a vector or matrix where the largest k (absolute) values plus any redundancies are preserved and the remaining values are converted to 0."
Entropy2::usage="Entropy2[v] is equivalent to the Mathematica function call Entropy[2,v] and returns the entropy (base 2) of input vector v."
MSE::usage="MSE[a,b] returns the mean squared error between the two input vectors or matrices."
PSNR::usage="PSNR[a,b] returns the peak signal-to-noise ratio between two input vectors or matrices."
GammaCorrection::usage="GammaCorrection[a,r] takes a grayscale image matrix a and a positive exponent r and returns an r-gamma corrected matrix."
HistogramEQ::usage="HistogramEQ[a] takes a grayscale image matrix and returns a histogram-equalized version of the image matrix."
HistogramMatch::usage="HistogramMatch[target,reference] performs histogram matching of a reference matrix to a target matrix."
AutomatedThreshold::usage="AutomatedThreshold[x,tol] takes a data set x and a stopping tolerance tol and returns a threshold tau for dividing the values in x."
MakeHuffmanCodes::usage = "MakeHuffmanCodes[x] takes a vector, matrix or string and returns the Huffman codes for the input element as well as information about the old and new bitstream lengths."
HuffmanTree::usage="HuffmanTree[codes,opts] takes the Huffman codes generated by MakeHuffmanCodes (first argument) and makes a Huffman tree from them."


(* Option Values for image processing functions *)

ImageForm::usage="ImageForm is an option used by ImageInfo to designate between Color and GrayScale images."
Color::usage="Color is a value for the option ImageForm used by ImageInfo to indicate the images under consideration is an RGB color image."
GrayScale::usage="GrayScale is an option option used by ImageInfo and ImagePlot to indicate the image under consideration is a GrayLevel image."
ShowThumbnails::usage="ShowThumbnails is an option used by ImageInfo to indicate whether or not the function should display available images as thumbnails."
ThumbnailColumns::usage="ThumbnailColumns is the number of columns for the table of thumbnail images displayed by ImageInfo."
ThumbnailSize::usage="ThumbnailSize is the width of each thumbnail image displayed in ImageInfo."
DisplayMode::usage="DisplayMode is an option for HistogramMatch, RGBToYCbCr, YCbCrToRGB, RGBToHSI and HSIToRGB that indicates whether or not the returned values are in the range 0,...,255."
Resize::usage="Resize is an option for ImageRead that can be used to resize the input image."
NodeColor::usage="NodeColor is an option for HuffmanTree that gives the background color for the nodes."
NodeEdgeColor::usage="NodeEdgeColor is an option for HuffmanTree that gives the edge color for the nodes."
NodeEdgeThickness::usage="NodeEdgeThickness is an option for HuffmanTree that sets the thickness of the edge of the nodes."
NodeSize::usage="NodeSize is an option for HuffmanTree that sets the radius of the nodes in the tree."
NodeFontSize::usage="NodeFontSize is an option for HuffmanTree that sets the font size when labeling nodes in the tree."
NodeFontColor::usage="NodeFontColor is an option for HuffmanTree that sets the font color when labeling nodes in the tree."
BranchLength::usage="BranchLength is an option for HuffmanTree that sets the lengths between nodes."
BranchColor::usage="BranchColor is an option for HuffmanTree that sets the color of the lines between nodes."
BranchFontSize::usage="BranchFontSize is an option for HuffmanTree that sets the font size when labeling edges."
BranchFontColor::usage="BranchFontColor is an option for HuffmanTree that sets the font color when labeling edges."
BranchThickness::usage="BranchThickness is an option for HuffmanTree that sets the thickness of the branches."
ShowString::usage="ShowString is an option for HuffmanTree that indicates whether or not string input should be displayed."


(* Options for image processing functions *)

Options[ImageInfo]={ImageForm->All,ShowThumbnails->False,ThumbnailColumns->3,ThumbnailSize->150,DisplayInfo->False};
Options[ImageRead]={DisplayInfo->False,PowersOfTwo->0,Resize->False};
Options[RGBToYCbCr]={DisplayMode->False};
Options[YCbCrToRGB]=Options[RGBToYCbCr];
Options[RGBToHSI]=Options[RGBToYCbCr];
Options[HSIToRGB]=Options[RGBToYCbCr];
Options[CE]={Computation->Numerical};
Options[HistogramMatch]={DisplayMode->False};
Options[HuffmanTree]={NodeColor->Lighter[Gray],NodeEdgeColor->Black,NodeEdgeThickness->Medium,NodeSize->.2,NodeFontSize->12,NodeFontColor->Darker[Gray],BranchLength->.1,BranchColor->Black,BranchFontSize->12,BranchFontColor->Black,BranchThickness->Medium,ShowString->"",ImageSize->Automatic};

Begin["`Private`"]
(* Implementation of the package *)

(* Image processing functions *)

ImageFormats={"GIF","JPEG","TIFF","PNG","WebP","JPEG2000","OpenEXR","ICC","BMP","PICT","WMF","EMF","ACO","XBM","CUR","ICO","ICNS","PCX","PBM","PGM","PPM","PNM","PXR","TGA","SCT","GeoTIFF","FITS","DICOM","HDF","HDF5","NASACDF","NetCDF","RawBitmap"};

ImageInfo[OptionsPattern[]]:=Module[{HomeDirectory=Directory[],ImageDirectory,files,filenames,colorgray,thumbcols,thumbsize,thumbs,labels,pad,f,k,size,type,its},
	
	(*
	WaveletWareFile=Complement[Map[FindFile[#<>"WaveletWare.m"]&,Map[StringJoin[#,"/WaveletWare/"]&,$Path]],{$Failed}];
	If[Length[WaveletWareFile]==0,Message[ImageInfo::"packagenotfound"];Return[{}]];
	PackageDirectory=DirectoryName[First[WaveletWareFile]];
	*)
	(* Next retrieve list of all image files *)
	(*ImageDirectory=PackageDirectory<>"Images/";*)
	ImageDirectory=PackageDirectory[DataType->Images];
	SetDirectory[ImageDirectory];
	files=Map[StringJoin[ImageDirectory,#]&,Union[FileNames[],FileNames["*",{"*"},Infinity]]];
	files=Select[files,MemberQ[ImageFormats,FileFormat[#]]&];
	
	(* Now determine if we we want color, grayscale or both. *)
	colorgray=OptionValue[ImageInfo,ImageForm];
	files=Which[colorgray===Color,Select[files,ImageChannels[Import[#]]==3&],colorgray===GrayScale,Select[files,ImageChannels[Import[#]]==1&],True,files];
	SetDirectory[HomeDirectory];
	
	(* Generate thumbnail display if requested. *)
	If[TrueFalse[OptionValue[ImageInfo,ShowThumbnails]],
		thumbcols=Min[Length[files],OptionValue[ImageInfo,ThumbnailColumns]];
		thumbsize=OptionValue[ImageInfo,ThumbnailSize];
		(* For display purposes, grab just the names of the files. *)
		filenames=Map[FileNameDrop[#,FileNameDepth[ImageDirectory]]&,files];
		pad=Which[Mod[Length[filenames],thumbcols]==0,0,True,thumbcols-Mod[Length[filenames],thumbcols]];
		thumbs=Partition[PadRight[Map[ImageResize[Import[#],thumbsize]&,files],Length[filenames]+pad," "],thumbcols];
		labels=Partition[PadRight[MapThread["Image"<>ToString[#1]<>"="<>#2&,{Range[Length[filenames]],filenames}],Length[filenames]+pad," "],thumbcols];
		f=Map[#->True&,Tuples[{Partition[Range[2*Length[Flatten[thumbs]]/thumbcols],2],Table[{k,k},{k,1,thumbcols}]}]];
		Print[GraphicsGrid[Riffle[thumbs,labels],ImageSize->thumbsize*thumbcols,Frame->{None,None,f}]];
	];
	
	(* Display image info if requested. *)
	
	If[TrueFalse[OptionValue[ImageInfo,DisplayInfo]],
		{size,type}=Transpose[Map[Import[#,{{"ImageSize","ColorSpace"}}]&,files]];
		size=Map[Reverse,size];
		its=Map[Min, Map[IntegerExponent[#, 2] &, size, 1]];
		size=MapThread[ToString[#1]<>"x"<>ToString[#2]&,Transpose[size]];
		filenames=Map[FileNameDrop[#,FileNameDepth[ImageDirectory]]&,files];
		Print["Information for Images in the WaveletWare Package\n",TableForm[Transpose[{Range[Length[filenames]],filenames,size,type,its}],TableHeadings->{None,{"Image Number","Filename","Dimensions","Type","Max Iterations"}},TableAlignments->{Center,Center},TableSpacing->{5,5}]];
	];
	Return[files];
];

ImageRead[file_,OptionsPattern[]]:=Module[{img,rows,cols,type,scale,p,data},
	(* Check validity of given file and read in image if okay. *)
	img=Import[file];
	If[!ImageQ[img],Message[ImageRead::"badformat",file];Return[{}]];
	
	{rows,cols}=ImageDimensions[img];
	type=ImageChannels[img];
	
	(* Resize the image if so instructed *)
	scale=OptionValue[ImageRead,Resize];
	If[scale===False,
		0,
		If[!((IntegerQ[scale] && scale>0) || (VectorQ[scale,IntegerQ] && Length[scale]==2 && Min[scale]>0)),Message[ImageRead::"badinput","Resize"," either False (no resizing) or a positive integer or a list of two positive integers"];Return[{}]];
		If[VectorQ[scale]&&Length[scale]==2,scale=Reverse[scale]];
		img=ImageResize[img,scale];
		{rows,cols}=ImageDimensions[img];
	];
	
	(* Crop image so that the dimensions are divisible by the prescribed power of 2 *)
	p=OptionValue[ImageRead,PowersOfTwo];
	If[!(IntegerQ[p] && p>=0),Message[ImageRead::"badinput","PowersOfTwo"," a nonnegative integer"];Return[{}]]; 
	{rows,cols}=Map[Which[Mod[#, 2^p] == 0, #, True, If[# < 2^p, 2^p, # - Mod[#, 2^p]]]&,{rows,cols}];
	img=ImageCrop[img,{rows,cols}];
	
	(* Extract data into matrix channels and return *)
	data=ImageData[img];
	Which[type==1,data=Round[255*data],type==3,data=Round[255*Map[data[[All,All,#]]&,Range[3]]]]; (* {Map[data[[All,All,#]]&,Range[3]]} *)

	p=If[type==1,MaxIts[data],MaxIts[First[data]]];
	(* Print image info if requested *)
	If[TrueFalse[OptionValue[ImageRead,DisplayInfo]],
		Print["The returned image has dimensions ",cols,"x",rows," and was processed as a ",Which[type==1,"grayscale ",type==3,"color "],"image.  A maximum of ",p," iterations of a wavelet transformation can be performed on this image."]];
	
	Return[data];
];

ImageEnlarge[a_,s_]:=Module[{wt,log,chk,d},
	chk=First[CheckData[a]];
	If[!MemberQ[{2,3},chk],Message[ImageEnlarge::"badinput"];Return[a];];
	d=Which[chk==2,Dimensions[a],True,Dimensions[First[a]]];
	log=Log[2,s];
	If[!(IntegerQ[log]&&Positive[log]),Message[ImageEnlarge::"badfactor"];Return[a];];
	(*If[!MatrixQ[a],Message[ImageEnlarge::"badinput"];Return[a];];*)
	wt=Which[chk==2,PadRight[a,s*d],True,Map[PadRight[#,s*d]&,a]];
	Return[InverseHWT[wt,NumIterations->log,Computation->Symbolic]];
];

RGBToYCbCr[p_,OptionsPattern[]]:=Module[{dmode,m={{.299,.587,.114},{-.299,-.587,.886}/1.772,{.701,-.587,-.114}/1.402},s=DiagonalMatrix[{219,224,224}],t={16,128,128},data},
	(* Check for bad input *)
	If[Length[p]!=3,Message[RGBToYCbCr::"badinput"];Return[{}]];
	(* Get DisplayMode *)
	dmode = TrueFalse[OptionValue[RGBToYCbCr,DisplayMode]];
	
	(* First do the ordered triple case *)
	If[(VectorQ[p] && ArrayQ[p,_,IntegerQ]==True && Max[p]<=255 && Min[p]>=0),
		data=m.p/255.;
		If[!dmode,Return[data],Return[s.data+t]];
	];
	
	(* Now do the list of matrices case *)
	If[(AllTrue[p,MatrixQ] && ArrayQ[p,_,IntegerQ] && AllTrue[Map[Max, p], # <= 255 &] && AllTrue[Map[Min, p], # >= 0 &] && AllTrue[Map[Dimensions, p], # == Drop[Dimensions[p], 1] &]),
		data=Map[m.#&,Transpose[Map[Flatten,p]]/255.];
		data=Return[Map[Partition[#,Last[Dimensions[p]]]&,Transpose[If[!dmode,data,Map[(s.#+t)&,data]]]]];
	];
	
	Message[RGBToYCbCr::"badinput"]; 
	Return[{}];
];

YCbCrToRGB[p_,OptionsPattern[]]:=Module[{dmode,m=Inverse[{{.299,.587,.114},{-.299,-.587,.886}/1.772,{.701,-.587,-.114}/1.402}],s=DiagonalMatrix[1./{219,224,224}],t={16,128,128},data,intervalcheck},
	(* Check for bad input *)
	If[Length[p]!=3,Message[YCbCrToRGB::"badinput"];Return[{}]];
	
	(* Get DisplayMode *)
	dmode = TrueFalse[OptionValue[YCbCrToRGB,DisplayMode]];
		
	intervalcheck[i_,c_]:=AllTrue[Boole[MapThread[IntervalMemberQ[Interval[#1],#2]&,{i,c}]],(#==1)&];
	
	(* First do the ordered triple case *)
	If[VectorQ[p],
		(* Check that the input is in the correct domain, depending on dmode *)
		If[(!dmode && !intervalcheck[{{0, 1}, {-1/2, 1/2}, {-1/2, 1/2}},p]),
			Message[YCbCrToRGB::"baddisplay","False","[0,1],[-1/2,1/2],[-1/2,1/2]"];Return[{}]];
		If[(dmode && !intervalcheck[{{16, 235}, {16,240}, {16, 240}},p]),
			Message[YCbCrToRGB::"baddisplay","True","[16,235],[16,240],[16,240]"];Return[{}]];
	
		Return[If[!dmode,Round[255*m.p],Round[255*m.s.(p-t)]]]
		];
	
	If[AllTrue[p,MatrixQ],
		(* Now check that the input are three matrices of the same dimensions *)
		If[(AllTrue[Map[Dimensions, p],(#==Drop[Dimensions[p],1]&)])===False,Message[YCbCrToRGB::"badinput"];Return[{}]];

		data=Transpose[Map[Flatten,p]];
		(* Check that the input is in the correct domain, depending on dmode *)
		If[(!dmode && MemberQ[Map[intervalcheck[{{0, 1}, {-1/2, 1/2}, {-1/2, 1/2}},#]&,data],False]),
			Message[YCbCrToRGB::"baddisplay","False","[0,1],[-1/2,1/2],[-1/2,1/2]"];Return[{}]];
		If[(dmode && MemberQ[Map[intervalcheck[{{16, 235}, {16,240}, {16,240}},#]&,data],False]),
			Message[YCbCrToRGB::"baddisplay","False","[0,1],[-1/2,1/2],[-1/2,1/2]"];Return[{}]];
				
		data=If[!dmode,Round[Map[255*m.#&,data]],Round[Map[(255*m.s.(#-t))&,data]]];
		Return[Map[Partition[#,Last[Dimensions[p]]]&,Transpose[data]]];
	];
	
	Message[YCbCrToRGB::"badinput"]; 
	Return[{}];

];


RGBToHSI[a_,OptionsPattern[]]:=Module[{sc,r,g,b,intensity,top,bot,t,saturation,hue,rows,cols,p,redgeqb,glessb},
	
	sc=If[TrueFalse[OptionValue[RGBToHSI,DisplayMode]],255,1];
	(* First do the ordered triple case *)
	If[(Length[a]==3 && VectorQ[a] && ArrayQ[a,_,IntegerQ] && Max[a]<=255 && Min[a] >=0),
		intensity = Total[a]/3/255.;
		{r,g,b}=a;
		If[r==g==b, Return[{0., 0., intensity}]];
		top = (2*r - g - b)/2;
  		bot = ((r - g)^2 + (r - b)*(g - b));
  		t = ArcCos[top/Sqrt[bot]];
  		saturation = N[1 - Min[a]/intensity/255];
  		hue = N[If[g >= b, t, 2*Pi - t]/2/Pi];
  		Return[sc*{hue, saturation, intensity}];
  	];
		
	(* Now do the case where the input values are R,G,B matrices. *)
	If[(Length[a]==3 && AllTrue[a,MatrixQ] && AllTrue[IntegerQ,a] && Map[Dimensions,a]==ConstantArray[Dimensions[First[a]],3] && Max[a]<=255 && Min[a]>=0),
		intensity = Total[a]/3/255.;
  		{r, g, b} = a;
  		{rows, cols} = Drop[Dimensions[a], 1];
		top = (2*r - g - b)/2;
		bot = Sqrt[((r - g)^2 + (r - b)*(g - b))];
		p = Position[bot, 0];
		redgeqb = ReplacePart[ConstantArray[1, {rows, cols}], p -> 0];
		bot = ReplacePart[bot, p -> 1000];
		hue = ArcCos[N[top/bot]]*redgeqb;
		glessb = Boole[MapThread[Less[#1, #2] &, {g, b}, 2]];
		hue = 2*Pi*glessb + hue*(-1)^glessb;
		saturation = redgeqb*(1 - Partition[Map[Min, Transpose[Map[Flatten, {r, g, b}]]], cols]/ReplacePart[intensity, Position[intensity, 0.] -> 1]/255.);
		Return[sc*{hue/2/Pi, saturation, intensity}];
	];
	
	Message[RGBToHSI::"badinput"];
	Return[{}];
];

HSIToRGB[a_,OptionsPattern[]]:=Module[{sc,vectorq,matrixq,d,h,s,i,deg,idx,angle,t,u,v,cols,gray,p,s1,s2,s3,r,g,b},
	
	sc=If[TrueFalse[OptionValue[HSIToRGB,DisplayMode]],255,1];
	
	(* Determine whether or not the input is a list of triples or a list of matrices *)
	vectorq=(Length[a]==3 && VectorQ[a] && ArrayQ[a,_,NumericQ] && Max[a/sc]<=1 && Min[a] >=0);
	matrixq=(Length[a]==3 && AllTrue[a,MatrixQ] && AllTrue[Map[Dimensions, a],(#==Drop[Dimensions[a],1]&)] && Max[a/sc]<=1 && Min[a]>=0);
	
	(* Compute some variable for either a list of triples or a list of matrices *)
	If[vectorq || matrixq,
		d = Dimensions[a];
		{h, s, i} = a*{2*Pi,1,1}/sc;
		deg = 180*h/Pi;
		idx = (deg - Mod[deg, 120])/120;
		angle = h - idx*2*Pi/3;
	
		t = (1 - s)/3;
		u = (1 + s*Cos[angle]/Cos[Pi/3 - angle])/3;
		v = 1 - t - u;
	];
	If[vectorq,
		If[h == 0 && s == 0,Return[Round[ConstantArray[255*i, {3}]]]]; 
		Return[Round[255*3*i*Which[(0 <= h && h < 2*Pi/3), {u, v, t}, (2*Pi/3 <= h && h < 4*Pi/3), {t, u, v}, True, {v, t, u}]]];
	];
	
	If[matrixq,
		cols = Last[d];
		gray = ConstantArray[1, Drop[d, 1]];
		p = Intersection[Position[h, 0], Position[s, 0]];
		gray = ReplacePart[gray, p -> 0];
			
		d = Flatten[h];
		s1 = Partition[Boole[Map[GreaterEqual[#, 0] &, d]]*Boole[Map[Less[#, 2*Pi/3] &, d]], cols];
		s2 = Partition[Boole[Map[GreaterEqual[#, 2*Pi/3] &, d]]*Boole[Map[Less[#, 4*Pi/3] &, d]], cols];
		s3 = Partition[Boole[Map[GreaterEqual[#, 4*Pi/3] &, d]], cols]*Partition[Boole[Map[Less[#, 2*Pi] &, d]], cols];
		
		r = (s1*u + s2*t + s3*v)*gray + (1 - gray)/3;
		g = (s1*v + s2*u + s3*t)*gray + (1 - gray)/3;
		b = (s1*t + s2*v + s3*u)*gray + (1 - gray)/3;
	
		Return[Round[255*3*{i*r, i*g, i*b}]];
	];
	
	Message[HSIToRGB::"badinput"];
	Return[{}];
];

DCT[a_]:=Module[{r,c,d1,d2},
	If[VectorQ[a],Return[Join[{1},ConstantArray[Sqrt[2],Length[a]-1]]*FourierDCT[a,2]]];
	
	If[MatrixQ[a],
		{r,c}=Dimensions[a];
		d1=DiagonalMatrix[Join[{1},ConstantArray[Sqrt[2],r-1]]];
		d2=DiagonalMatrix[Join[{1},ConstantArray[Sqrt[2],c-1]]];
		Return[d1.FourierDCT[a,2].d2];
	];
	
	Message[DCT::"badinput"];
	Return[{}];
];

InverseDCT[a_]:=Module[{r,c,d1,d2},
	If[VectorQ[a],
   		Return[FourierDCT[Join[{1},ConstantArray[1/Sqrt[2],Length[a]-1]]*a,3]]
	];
	
	If[MatrixQ[a],
		{r,c}=Dimensions[a];
		d1=DiagonalMatrix[Join[{1},ConstantArray[1/Sqrt[2],r-1]]];
		d2=DiagonalMatrix[Join[{1},ConstantArray[1/Sqrt[2],c-1]]];
		Return[FourierDCT[d1.a.d2,3]];
	];
	
	Message[InverseDCT::"badinput"];
	Return[{}];
];

CE[v_,OptionsPattern[]]:=Module[{w,m,p},
	If[!((VectorQ[v] && Length[v]>0) || (MatrixQ[v] && Min[Dimensions[v]]>0)),Message[CE::"badinput"];Return[{}]];
	
	If[Max[Abs[Flatten[v]]]==0,Message[CE::"zeroinput"];Return[{}]];
	
	w=Accumulate[Sort[Flatten[v]^2, Greater]]/Norm[Flatten[v]]^2;
	m=OptionValue[CE,Computation];
    p=If[m===Symbolic,1,1.];
    Return[p*w];
];

PercentCE[v_,p_]:=Module[{s},
	If[!(VectorQ[v] && Min[v]>0 && Round[Max[v]]==1),Message[PercentCE::"badinput"];Return[{}]];
	If[v!=Sort[v],Message[PercentCE::"badinput"];Return[{}]];
	If[!(NumericQ[p] && p>0 && p<=1),Message[PercentCE::"badpercent"];Return[{}]];
	
	s=Select[v,#<p&];
	Return[Length[s]];
];

nCE[v_,p_]:=PercentCE[v,p];

QuantizeCE[v_,k_]:=Module[{r,a,b},
	If[!((VectorQ[v] && Length[v]>0) || (MatrixQ[v] && Min[Dimensions[v]]>0)),Message[QuantizeCE::"badinput"];Return[{}]];
	If[!(IntegerQ[k] && k>=0 && k<=Length[Flatten[v]]),Message[QuantizeCE::"badinteger"];Return[{}]];
	
	r=If[AllTrue[Flatten[v],ExactNumberQ],1,0];
	
	a=Abs[Sort[Flatten[Abs[v]],Greater][[k]]];
	b=Chop[N[v],N[a]];
	Return[If[r==1,Rationalize[b],b]];	
];

Comp[v_,k_]:=QuantizeCE[v,k];

Entropy2[v_]:=Entropy[2,Flatten[v]];

MSE[a_,b_]:=Module[{},
	If[!(ArrayQ[a]&&ArrayQ[b]&&Dimensions[a]==Dimensions[b]&&Min[Dimensions[a]]>0),Message[MSE::"badinput"];Return[{}]];
	Return[Total[Flatten[(a-b)^2]]/Apply[Times,Dimensions[a]]];
];
	
PSNR[a_,b_] := If[MSE[a,b]==0,Return[Infinity],Return[10.*Log[10, 255^2/MSE[a,b]]]];

GammaCorrection[a_,r_]:=Module[{},
		If[MatrixQ[a,IntegerQ]==False,Message[GammaCorrection::"badinput"];Return[a];];
		If[Min[a]<0 || Max[a]>255,Message[GammaCorrection::"badinput"];Return[a];];
		If[NumericQ[r]==False,Message[GammaCorrection::"badvalue"];Return[a];];
		If[r<=0,Message[GammaCorrection::"badvalue"];Return[a];];
		Return[Round[((a/255)^r)*255]];
];
	
HistogramEQ[a_]:=Module[{r,c,t,val,tot,l=Range[0,255],p,z,pairs,rules},
	If[MatrixQ[a,IntegerQ]==False,Message[HistogramEQ::"badinput"];Return[a];];
	If[Min[a]<0 || Max[a]>255,Message[HistogramEQ::"badinput"];Return[a];];
   	{r,c}=Dimensions[a];
   	t=Split[Sort[Flatten[a]]];
   	{val,tot}={Map[First,t],Map[Length,t]};
   	p=Complement[l,val];
   	z=ConstantArray[0,Length[p]];
   	pairs=Transpose[{Join[p, val], Join[z, tot]}];
   	t=Floor[255*Accumulate[Last[Transpose[SortBy[pairs, First]]]]/r/c];
   	rules=MapThread[#1->#2&,{l,t}];
   	Return[a/.rules];
 ];
 
 HistogramMatch[ref_,target_,OptionsPattern[]]:=Module[{rf,tgt,refcdf,targetcdf,v,p,f,q,rules,display},
 	
 	Map[If[MatrixQ[#,NumberQ]==False,Message[HistogramMatch::"badinput"];Return[ref];]&,{ref,target}];
 	Map[If[Min[#]<0 || Max[#]>1,Message[HistogramMatch::"badinput"];Return[ref];]&,{ref,target}];
 	
 	{rf,tgt}=Map[Round[255*#]&,{ref,target}];
 	{refcdf,targetcdf}=Map[Accumulate[BinCounts[Flatten[#],{0,256,1}]]/Apply[Times,Dimensions[#]]&,{rf,tgt}];
 	v=Map[Abs[refcdf-#]&,targetcdf];
 	p=Map[Flatten,Map[Position[#,Min[#]]&,v]-1];
 	
 	f[l_,r_]:=Module[{m, s, i},
  		m = First[Flatten[Position[l, r]]];
  		s = Abs[r - m];
  		i = First[Flatten[Position[s, Min[s]]]];
  		Return[r[[i]]];
  		];
 		
 	q = Map[Which[Length[#] == 1, First[#], True, f[p, #]] &, p];
	rules = MapThread[#1 -> #2 &, {Range[0, 255], q}];
	display=TrueFalse[OptionValue[HistogramMatch,DisplayMode]];
	v=tgt/.rules;
	Return[If[display,v,v/255.]];	
 ];
 
 AutomatedThreshold[x_,tol_]:=Module[{alpha,tau,S,S1,S2,s1,s2,tau1},
 	
 	If[!(VectorQ[x,NumberQ] || MatrixQ[x,NumberQ]),Message[AutomatedThreshold::"badinput"];Return[0];];
 	If[!(NumberQ[tol]&&Positive[tol]),Message[AutomatedThreshold::"badtolerance"];Return[0];];
 	alpha=2*tol;
 	S=Flatten[N[x]];
 	tau=Mean[{Min[S],Max[S]}];
 	While[alpha>tol,
 		S1=Select[S,#<tau&];
 		S2=Complement[S,S1];
 		{s1,s2}=Map[Mean,{S1,S2}];
 		tau1=(s1+s2)/2;
 		alpha=Abs[tau1-tau];
 		tau=tau1;
 	];
 	Return[tau];
 ];

MakeHuffmanCodes[v_]:=Module[{ch,chars,totals,c,freq,HStep,HList,HList1,tbl,bin,y},
 	
	If[!(MatrixQ[v,(IntegerQ[#]&&NonNegative[#])&]||VectorQ[v,(IntegerQ[#]&&NonNegative[#])&]||StringQ[v]),Message[MakeHuffmanCodes::"badinput"];Return[{}]];

	ch = Characters[Which[ArrayQ[v],FromCharacterCode[Flatten[v]],StringQ[v],v,True,{}]];
	If[Length[ch]==1,Return[{{{ToString[First[ch]],1.,{0}}},8,1}]];

	c=Sort[Map[Reverse,Tally[ch]]];
	chars=Map[Last,c];
	totals=Map[First,c];
	freq=N[totals]/Length[ch];

	HStep[a_]:={Total[First[Transpose[Take[a,2]]]],Take[a,2]};
	HList[a_]:=Sort[Append[Drop[a,2],HStep[a]]];
	HList1[a_]:=Module[{d=Drop[a,2],nl,p},
		nl=Map[First,d];
		p=LengthWhile[y,#<First[HStep[a]]&];
		Return[Insert[d,HStep[a],p]];
	];
	
	tbl=Nest[HList,c,Max[0,Length[c]-2]];

	bin=Map[Take[#,{1,Length[#],2}]&,Map[Flatten[Position[tbl,#]]&,chars]]-1;

	If[ArrayQ[v],chars=Flatten[ToCharacterCode[chars]]];

	Return[{Transpose[{chars,freq,bin}],8*Total[totals],Map[Length,bin].totals}];
];

HuffmanTree[codes_,OptionsPattern[]]:=Module[{l,p,b,nc,nec,net,ns,nfs,nfc,el,ec,efs,efc,bt,lvls,node,pr,lbl,k,idx,tree={},j,plt,ss,str={},imgsize},
	If[Length[Dimensions[codes]]!=2 || Last[Dimensions[codes]]!=3,Message[HuffmanTree::"badinput"];Return[{}];];
	
	{l,p,b}=Transpose[codes];
	(* Check the validity of l,p,b *)
	If[!AllTrue[{l,p},VectorQ],Message[HuffmanTree::"badinput"];Return[{}];];
	If[!(AllTrue[l,StringQ]||AllTrue[l,IntegerQ]),Message[HuffmanTree::"badinput"];Return[{}];];
	If[!(AllTrue[p,NumberQ]&&AllTrue[p,#>=0&]&&AllTrue[p,#<=1&]),Message[HuffmanTree::"badinput"];Return[{}];];
	If[!AllTrue[b,VectorQ] || !AllTrue[Flatten[b],MemberQ[{0,1},#]&],Message[HuffmanTree::"badinput"];Return[{}];];
	
	{nc,nec,net,ns,nfs,nfc,el,ec,efs,efc,bt,ss,imgsize}=OptionValue[HuffmanTree,{NodeColor,NodeEdgeColor,NodeEdgeThickness,NodeSize,NodeFontSize,NodeFontColor,BranchLength,BranchColor,BranchFontSize,BranchFontColor,BranchThickness,ShowString,ImageSize}];
	
	(* Check integrity of options *)
	If[!ColorQ[nc],Message[HuffmanTree::"badcolor","NodeColor"];nc=NodeColor/.Options[HuffmanTree];];
	If[!ColorQ[nec],Message[HuffmanTree::"badcolor","NodeEdgeColor"];nec=NodeEdgeColor/.Options[HuffmanTree];];
	If[!ColorQ[nfc],Message[HuffmanTree::"badcolor","NodeFontColor"];nfc=NodeFontColor/.Options[HuffmanTree];];
	If[!ColorQ[ec],Message[HuffmanTree::"badcolor","BranchColor"];nfc=BranchColor/.Options[HuffmanTree];];
	If[!ColorQ[efc],Message[HuffmanTree::"badcolor","BranchFontColor"];nfc=BranchFontColor/.Options[HuffmanTree];];
	If[!Positive[ns]||!NumericQ[ns],Message[HuffmanTree::"badsize","NodeSize"];ns=NodeSize/.Options[HuffmanTree];];
	If[!Positive[el]||!NumericQ[el],Message[HuffmanTree::"badsize","BranchLength"];el=NodeSize/.Options[HuffmanTree];];
	If[!StringQ[ss],Message[HuffmanTree::"badstring"];ss=""];
	
	str=Style[ss,Gray,Plain,24];
	p=SetPrecision[p,2];
	
	(* Plot the case where the number of characters is 1. *)
	If[Length[l]==1,
		Return[TreePlot[{1->1},Top,1,VertexLabeling->True,VertexRenderingFunction->({nc,EdgeForm[{Thickness[net],nec}],Disk[#,ns],Text[Style[l[[1]],nfc,nfs],#1]}&),PlotLabel->str,ImageSize->imgsize]];
	];

	lvls=Max[Map[Length,b]]+1;
	node=Table[0,{2^lvls-1}];
	pr=node;
	lbl=Table[" ",{2^lvls-1}];

	For[k=1,k<=Length[l],k++,
		idx=2^Length[b[[k]]]+FromDigits[b[[k]],2];
		node[[idx]]=idx;
		pr[[idx]]=p[[k]];
		lbl[[idx]]=ToString[l[[k]]]<>"\n"<>ToString[pr[[idx]]];
	];

	For[j=lvls-1,j>=1,j--,
		For[k=0,k<2^j,k+=2,
			idx=2^j+k;
			If[node[[idx]]!=0,
				node[[idx/2]]=idx/2;
				pr[[idx/2]]=pr[[idx]]+pr[[idx+1]];
				lbl[[idx/2]]=" \n"<>ToString[pr[[idx/2]]];
				tree=Prepend[tree,{idx/2->idx,"0"}];
				tree=Prepend[tree,{idx/2->idx+1,"1"}];
			];
		];
	];
	plt=TreePlot[Sort[tree],Top,1,LayerSizeFunction-> (#^el&),VertexRenderingFunction->({nc,EdgeForm[{Thickness[net],nec}],Disk[#,ns],Text[Style[lbl[[#2]],nfc,nfs],#1]}&),EdgeRenderingFunction -> ({Thickness[bt],Text[Style[#3, efc, FontSize -> efs], 
     Mean[#1], Background -> White], ec, Line[#]} &),PlotLabel->str,ImageSize->imgsize];
	Return[plt];
];


(* Error messages *)

ImageInfo::packagenotfound="WaveletWare package not found in any $Path folders.  Check installation instructions."
ImageRead::badformat="The format `1` is not one of the accepted Mathematica image raster formats - ImageRead failed." 
ImageRead::badinput="The value of `1` should be `2` - ImageRead failed."
ImageEnlarge::badfactor="The scale value must be a power of 2 - ImageEnlarge failed."
ImageEnlarge::badinput="The input must be a matrix - ImageEnlarge failed."
RGBToYCbCr::badinput="The input for RGBToYCbCr is either an r,g,b triple or a list of R,G,B matrices - RGBToYCbCr failed."
YCbCrToRGB::badinput="The input for YCbCrToRGB is either a y,cb,cr triple or a list of Y,Cb,Cr matrices - YCbCrToRGB failed."
YCbCrToRGB::baddisplay="For DisplayMode set to `1` the range of Y,Cb,Cr must be the intervals `2` - YCbCrRGB failed."
RGBToHSI::badinput="The input for RGBToHSI is either an r,g,b triple or a list of R,G,B matrices whose values are all integers 0,...,255 - RGBToHSI failed."
HSIToRGB::badinput="The input for HSIToRGB is either an h,s,i triple or a list of H,S,I matrices - HSIToRGB failed."
DCT::badinput="The input must either be a vector or a matrix - DCT failed."
InverseDCT::badinput="The input must either be a vector or a matrix - InverseDCT failed."
CE::badinput="The input must either be a vector or matrix - CE failed."
CE::zeroinput="The cumulative energy vector cannot be computed for a zero vector or matrix - CE failed."
PercentCE::badinput="The input must be a vector whose values are nondecreasing and contained in the interval (0,1] - PercentCE failed."
PercentCE::badpercent="The value must be a number in the interval (0,1] - PercentCE failed."
QuantizeCE::badinput="The input must either be a vector or matrix - QuantizeCE failed."
QuantizeCE::badinteger="The second input must be a nonnegative integer whose length is at most the total length of the elements in the first input - QuantizeCE failed."
MSE::badinput="The input must be two vectors or matrices of equal dimensions."
GammaCorrection::badinput="The input matrix must be comprised of integers between 0 and 255 inclusive."
GammaCorrection::badvalue="The exponent must be a positive numeric value."
HistogramEQ::badinput="The input matrix must be comprised of integers between 0 and 255 inclusive."
HistogramMatch::badinput="The inputs must be matrices whose values are numbers between 0 and 1 inclusive."
AutomatedThreshold::badinput="The input must be a vector or matrix of numbers - AutomatedThreshold failed."
AutomatedThreshold::badtolerance="The tolerance must be a positive number - AutomatedThreshold failed."
MakeHuffmanCodes::badinput="The input must be an integer-valued vector or matrix or a string - MakeHuffmanCodes failed."
HuffmanTree::badinput="The format of the input must be the same as the first item returned by MakeHuffmanCodes - HuffmanTree failed."
HuffmanTree::badcolor="Warning :: The value given for `1` is not a valid Mathematica color - resetting `1` to its default value."
HuffmanTree::badsize="Warning :: The value given for `1` is not valid - resetting `1` to its default value."
HuffmanTree::badstring="Warning :: The value is not a valid string - no string will be displayed with the plot."

End[]
EndPackage[]