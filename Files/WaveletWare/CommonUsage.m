(* Mathematica package *)

BeginPackage["WaveletWare`CommonUsage`"]
(* Exported symbols added here with SymbolName::usage *) 

(* Common Usage Options:

	Border 
	BorderColor
	Boundary
	ChannelColor
	CompositeTree
	Computation
	DisplayInfo
	IntegerMap
	Iteration
	Numerical
	NumIterations
	PowersOfTwo
	Shift
	Symbolic
*)

If[!ValueQ[DisplayInfo::usage],
	DisplayInfo::usage="DisplayInfo is an option that toggles the display of information in ImageInfo, ImageRead, AudioInfo, AudioRead, DataInfo and SplineFilters."
];

If[!ValueQ[Shift::usage],
	Shift::usage="Shift is an option for discrete wavelet transform functions to indicate if the lowpass/highpass filters should be shifted for the computation."
];

If[!ValueQ[NumIterations::usage],
	NumIterations::usage="NumIterations is an option used by several functions in the WaveletWare package.  It indicates the number of iterations of the wavelet transform performed on the data."
];

If[!ValueQ[ChannelColor::usage],
	ChannelColor::usage="ChannelColor is an option for ImagePlot, WaveletPlot, FullWaveletPlot and WaveletRegionPlot.  It can be used to set ColorFunction to a particular RGB color when a single matrix (with Scaling set to Image) is input in ImagePlot."
];

If[!ValueQ[Border::usage],
	Border::usage="Border is an option for ImagePlot, WaveletPlot, FullWaveletPlot and WaveletRegionPlot.  If set to True (default), a border is placed around the plotted output."
];

If[!ValueQ[BorderColor::usage],
	BorderColor::usage="BorderColor is an option for ImagePlot, WaveletPlot, FullWaveletPlot and WaveletRegionPlot.  It sets the color of the image border."
];

If[!ValueQ[PowersOfTwo::usage],
	PowersOfTwo::usage="PowersOfTwo is an option that can be set in ImageRead or AudioRead so that the dimensions of the imported data is divisible by the given power of 2."
];

If[!ValueQ[Computation::usage],
	Computation::usage="Computation is an option used by wavelet transform functions, their inverses, de-noising functions and the and cumulative energy.  It indicates whether the comptuation should be done numerically or symbolically."
];

If[!ValueQ[Numerical::usage],
	Numerical::usage="Numerical is a symbol in the WaveletWare package that serves as a value of the Computation option for various wavelet transform functions."
];

If[!ValueQ[Symbolic::usage],
	Symbolic::usage="Symbolic is a symbol in the WaveletWare package that serves as a value of the Computation option for various wavelet transform functions."
];

If[!ValueQ[Boundary::usage],
	Boundary::usage="Boundary is an option for several routines that either compute or use the biorthogonal wavelet transform.  If set to True, a symmetric biorthogonal transformation is computed."
];

If[!ValueQ[IntegerMap::usage],
	IntegerMap::usage="IntegerMap is an option used by LeGall wavelet transform functions.  If set to True (default is False) and the input is integer-valued, then the output is integer-valued."
];

If[!ValueQ[Orthogonal::usage],
	Orthogonal::usage="Orthogonal is an option used by the Haar transform routines.  If set to True, the Haar filter is used to computed the Haar transform, else it uses {1/2,1/2} as the lowpass filter."
];

If[!ValueQ[CompositeTree::usage],
	CompositeTree::usage="CompositeTree is an option for WPT that when set to True (default is False) and the input is a list of three matrices, creates a single best basis tree."
];

If[!ValueQ[Iteration::usage],
	Iteration::usage="Iteration is an option for WaveletPlot, WaveletVectorPlot, FullWaveletPlot, FullWaveletVectorPlot and WaveletVectorPlay that can be used to designate as output a particular region of the input data."
];


EndPackage[]
