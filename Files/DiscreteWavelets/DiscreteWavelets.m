(* ::Package:: *)

(* : Title : DiscreteWavelets.m-- a package template*)
(* : Context : DiscreteWavelets`DiscreteWavelets`*)
(* : Author : Patrick J.Van Fleet*)
(* : Summary : This package contains routines that will allow the user to work with 
	 one and two dimensional wavelet transformations contains routines that will allow the user to visualize and/or hear digital images and signals and their wavelet transforms.*)
(* : Copyright : \[Copyright] 2007 by Patrick J.Van Fleet*)
(* : Package Version : 1.0*)
(* : Mathematica Version : 6.0*)
(* : History :*)
(* : Keywords : wavelets, filters, cumulative energy, 
  peak signal to noise ratio, entropy, Huffman codes, color spaces*)
(* : Sources :
      Patrick J.Van Fleet, 
  Discrete Wavelet Transformations : 
      An Elementary Approach With Applications.*)
(* : Warnings :
      This package only works with Mathematica version 6.0 or later!!*)
(* : Limitations :
      < special cases not handled, known problems >*)
(* : Discussion :
      < description of algorithm, information for experts >*)
(* : Requirements : None *)
(* : Examples :
      < sample input that demonstrates the features of this package >*)
      
 BeginPackage["DiscreteWavelets`DiscreteWavelets`"];

(*Functions usage messages*)

CE::usage = "CE returns the cumulative energy vector of the input vector v or \
matrix A."
Comp::usage = "Comp takes a vector v or a matrix A and an integer index k and \
first computes the kth largest element (in absolute value) in v.  
	If we call this value q, Comp then sets to 0 all values in v smaller (in \
absolute value) than q and returns this modified vector.  If k is larger than \
the number of elements in v or A, then k is reset to the number of elements \
in v or A.  If k is 0 or negative, Comp converts all elements to 0."
nCE::usage = "nCE takes a vector built from CE and a percentage p and returns \
an index k so that p% of the energy of v is stored in the first k components \
of v."
Entropy::usage = "Entropy takes a vector v or a matrix A and returns the \
entropy."
MSE::usage="MSE takes as input two matrices and returns their mean squared \
error."
PSNR::usage = "PSNR takes as input two matrices and returns the Peak Signal \
to Noise Ratio between the two matrices.  PSNR calls the MSE function."
RGBToYCbCr::usage="RGBToYCbCr converts from RGB to YCbCr.  The input can be a \
single rgb point {r,g,b}, a list of such points, or a list of three matrices \
(presumably red, green, and blue).  If DisplayMode is set to True, then the \
conversion is stretched to prevent roundoff error.  This is useful for \
displaying the three channels."
YCbCrToRGB::usage="YCbCrToRGB converts from YCbCr to RGB.  The input can be a \
single rgb point {y,u,v}, a list of such points, or a list of three matrices \
(presumably y, u, and v values).  If DisplayMode is set to True, then the \
conversion is stretched to prevent roundoff error.  This is useful for \
displaying the three channels."
MakeHuffmanCodes::usage = "MakeHuffmanCodes takes a string of characters, a \
vector of integers, or a matrix of integers, and returns the codes for the \
elements, the original bitstream length, and the bitstream length of the \
Huffman codes."
FiniteFourier::usage = "FiniteFourier creates a finite Fourier series from \
independent variable w, finite length list v, and a starting index."
Cs::usage = "Cs is simply the Cosine function written in terms of E^(Iw).  \
(See Exercise 4 in Sectin 4.2)."
Sn::usage = "Sn is simply the Sine function written in terms of E^(Iw).  (See \
Exercise 4 in Sectin 4.2)."
Haar::usage = "Haar returns the Haar filter."
Daub::usage = "Daub takes an even positive integer n as an argument and \
returns the Daubechies filter of length n."
Coif::usage="Coif takes a positive integer K as an argument and returns the \
Coiflet filter of length 6*K.  Currently this routines only works for K=1, 2, \
3."
SplineFilters::usage="SplineFilters takes two even integers M and Mt and \
returns the associated biorthogonal spline filter pair.  If only one integer \
is given, SplineFilters returns the Haar filter pair."
CDF97::usage="CDF97 returns the Cohen/Daubechies/Feauveau 9/7 filter pair \
(see Section 10.3)."
LeGall::usage="LeGall returns the LeGall (5,3) biorthogonal filter pair (see \
Section 12.2)."
HWT1D1::usage="HWT1D1 takes a vector v of even length and returns its Haar \
wavelet transform."
IHWT1D1::usage="IHWT1D1 takes a vector v of even length and returns its \
inverse Haar wavelet transform."
HWT1D::usage="HWT1D takes a vector v of length 2^p and an integer i between 0 \
and p inclusive and returns i iterations of the Haar Wavelets Transform."
IHWT1D::usage="IHWT1D takes a vector v of length 2^p and an integer i between \
0 and p inclusive and returns i iterations of the Inverse Haar Wavelet \
Transform."
WT1D1::usage="WT1D1 takes an even-length vector v and a filter h and returns \
one iteration of the wavelet transform.  If a filter h is not input, WT1D1 \
computes one iteration of the Haar wavelet transform."
IWTht::usage="IWTht is an auxiliary routine used by IWT1D1 to compute half of \
the inverse wavelet transform on input vector v."
IWT1D1::usage="IWT1D1 takes an even-length vector v and a filter h and \
returns one iteration of the inverse wavelet transform.  If a filter h is not \
input, IWT1D1 computes one iteration of the inverse Haar wavelet transform."
WT1D::usage="WT1D takes an even-length vector v and a filter h and performs \
an iterated wavelet transform.  The number of iterations is set as a \
directive.  If no directive is given, the routine computes one iteration.  If \
no filter is given, the routine computes the iterated Haar wavelet \
transform."
IWT1D::usage="IWT1D takes an even-length vector v and a filter h and performs \
an iterated inverse wavelet transform.  The number of iterations is set as a \
directive.  If no directive is given, the routine computes one iteration.  If \
no filter is given, the routine computes the iterated inverse Haar wavelet \
transform."
BWT1D1::usage="BWT1D1 takes an even-length vector v and a biorthogonal filter \
pair {h,hw} and returns one iteration of the biorthogonal wavelet transform."
IBWTht::usage="IBWTht is an auxiliary routine used by IBWT1D1 to compute half \
of the inverse wavelet transform on input vector v."
IBWT1D1::usage="IBWT1D1 takes an even-length vector v and a biorthogonal \
filter pair {h,hw} and returns one iteration of the inverse biorthogonal \
wavelet transform."
BWT1D::usage="BWT1D takes an even-length vector v and a filter pair h, hw and \
performs an iterated biorthogonal wavelet transform.  The number of \
iterations is set as a directive.  If no directive is given, the routine \
computes one iteration.  If no filter is given, the routine computes the \
iterated Haar wavelet transform."
IBWT1D::usage="IBWT1D takes an even-length vector v and a filter pair h, hw \
and performs an iterated inverse biorthogonal wavelet transform.  The number \
of iterations is set as a directive.  If no directive is given, the routine \
computes one iteration.  If no filter is given, the routine computes the \
iterated inverse Haar wavelet transform."
LWT1D1::usage="LWT1D1 takes an even-length vector v and employs lifting to \
return one iteration of the biorthogonal wavelet transform using the LeGall \
filter pair."
ILWT1D1::usage="ILWT1D1 takes an even-length vector v and employs lifting to \
return one iteration of the inverse biorthogonal wavelet transform using the \
LeGall filter pair."
LWT1D::usage="LWT1D takes an even-length vector v and employs lifting to \
return an iterated biorthogonal wavelet transform using the LeGall filter \
pair.  The number of iterations is set as a directive.  If no directive is \
given, the routine computes one iteration."
ILWT1D::usage="ILWT1D takes an even-length vector v and employs lifting to \
return an iterated inverse biorthogonal wavelet transform using the LeGall \
filter pair.  The number of iterations is set as a directive.  If no \
directive is given, the routine computes one iteration."
GetCorner::usage="GetCorner takes a matrix a and integers r and c and returns \
the r x c submatrix that resides in the upper left corner of a."
PutCorner::usage="PutCorner takes matrices m and b as input and replaces the \
upper left corner of m with b."
LeftHWT::usage="LeftHWT multiplies input matrix A on the left by the \
transpose of the Haar transform matrix."
RightHWT::usage="RightHWT multiplies the input matrix A on the right by the \
Haar transform matrix."
HWT2D1::usage="HWT2D1 returns one iteration of the two-dimensional Haar \
wavelet transform."
LeftIHWT::usage="LeftIHWT multiplies input matrix A on the left by the \
inverse Haar transform matrix."
IHWT2D1::usage="IHWT2D1 returns one iteration of the two-dimensional inverse \
Haar wavelet transform."
HWT2D::usage="HWT2D returns the iterated two-dimensional Haar wavelet \
transform."
IHWT2D::usage="HWT2D returns the iterated inverse two-dimensional Haar \
wavelet transform."
LeftWT::usage="LeftWT multiplies input matrix A on the left by the transpose \
of the wavelet transform matrix built from input filter h."
RightWT::usage="RighHWT multiplies the input matrix A on the right by the \
wavelet transform matrix built from input filter h."
WT2D1::usage="WT2D1 returns one iteration of the two-dimensional wavelet \
transform built from input filter h."
LeftIWT::usage="LeftIWT multiplies input matrix A on the left by the inverse \
wavelet transform matrix built from input filter h."
IWT2D1::usage="IWT2D1 returns one iterations of the two-dimensional inverse \
wavelet transform built from input filter h."
WT2D::usage="WT2D returns the iterated two-dimensional wavelet transform \
built from input filter h."
IWT2D::usage="IWT2D returns the iterated two-dimensional inverse wavelet \
transform built from input filter h."
BWT2D1::usage="BWT2D1 returns one iteration of the two-dimensional \
biorthogonal wavelet transform built from input filter pair {h, hw}."
IBWT2D1::usage="IBWT2D1 returns one iteration of the two-dimensional inverse \
biorthogonal wavelet transform built from input filter pair {h, hw}."
BWT2D::usage="BWT2D returns the iterated two-dimensional biorthogonal wavelet \
transform built from input filter pair {h,hw}."
IBWT2D::usage="IBWT2D returns the iterated two-dimensional inverse \
biorthogonal wavelet transform built from input filter pair {h,hw}."
LWT2D1::usage="LWT2D1 returns one iteration of the two-dimensional \
biorthogonal wavelet transform via lifting using the LeGall filter pair."
ILWT2D1::usage="ILWT2D1 returns one iteration of the two-dimensional inverse \
biorthogonal wavelet transform via lifting using the LeGall filter pair."
LWT2D::usage="LWT2D returns the iterated two-dimensional biorthogonal wavelet \
transform via lifting built from the LeGall filter pair."
ILWT2D::usage="IBWT2D returns the iterated two-dimensional inverse \
biorthogonal wavelet transform via lifting built from the LeGall filter \
pair."
DCT1D::usage="DCT1D computes the Type II (orthogonal) discrete cosine transform of the \
input vector."
IDCT1D::usage="IDCT1D computes the inverse of the Type II (orthogonal) discrete cosine \
transform of the input vector."
DCT2D::usage="DCT2D computes the two-dimensional Type II (orthogonal) discrete cosine \
transform of the input matrix."
IDCT2D::usage="IDCT2D computes the two-dimensional Type II inverse (orthogonal) discrete \
cosine transform of the input matrix."
WaveletMatrixToList::usage="WaveletMatrixToList takes a matrix produced by a \
2D wavelet transform function and creates a list.  The first component of the \
list is the lowpass portion while the remaining components are the highpass \
components.  Each of the highpass components have three elements."
WaveletListToMatrix::usage="WaveletListToMatrix takes a list of lowpass and \
highpass components (generated presumably by WaveletMatrixToList) and creates \
a matrix such as that returned by a 2D wavelet transform function."
WaveletVectorToList::usage="WaveletVectorToList takes a wavelet vector and \
converts it into a list whose components are the lowpass portion and the \
highpass portions of the transformation."
WaveletListToVector::usage="WaveletListToVector takes a list that has as its \
components the lowpass portion and the highpass portions and simply returns \
the Flattened list.  The command is equivalent to Flatten[v]."
ChopVector::usage="ChopVector will truncate a vector so that its new length \
has a factorization that contains n factors of 2."
GammaCorrection::usage="GammaCorrection takes a grayscale image matrix a and \
a value r and uses r to perfom gamma correction on a."
MakeHistogramEQ::usage="MakeHistogramEQ takes a grayscale image matrix a and \
returns a histogram that shows the distribution of intensities in the image."
HistogramEQ::usage="HistogramEQ takes a grayscale image matrix a and returns \
the histogram equalized image."
MAD::usage="MAD takes as input a matrix or vector and returns the Median Absolute Deviation."
ShrinkageFunction::usage="ShrinkageFunction takes as input a real number t and a tolerance lambda and either returns 0 or shrink the value lambda units closer to the t-axis."
DonohoSure::usage="DonohoSure takes as input a vector v and returns the SUREShrink tolerance lambda."
TestSparseness::usage="Testsparse takes a vector or matrix and determines whether or not it is sparse.  If it is sparse, the routine returns the UniversalThreshold tolerance, and if it is not sparse, the routine returns the SureShrink threshold tolerance."
NoiseEstimate::usage="NoiseEstimate takes a vector or matrix and an orthogonal filter and returns an estimate of the noise present."
UniversalThreshold::usage="UniversalThreshold takes a vector or matrix, an orthogonal filter, and a number of iterations for the wavelet transform and returns the universal threshold tolerance."
WaveletShrinkage::usage="WaveletShrinkage takes a vector or matrix, an orthogonal filter, a tolerance (or list of tolerances), and a number of iterations for the wavelet transform and returns a denoised version of the input vector or matrix.  The module uses Algorithm 9.1 from Section 9.1 of the book."
SureShrink::usage="SureShrink takes a vector or matrix, an orthogonal filter, and a number of iterations for the wavelet transform and uses the Sureshrink tolerance with Algorithm 9.1 to return a denoised version of the input."


ImageList::usage="ImageList gives information about images that come with the DiscreteWavelets package."
ImageNames::usage="ImageNames gives the absolute path name for each image in the DiscreteWavelets package."
ShowThumbnails::usage="ShowThumbnails plots thumbnails of a given list of images."
AudioList::usage="AudioList gives information about audio files that come with the DiscreteWavelets package."
AudioNames::usage="AudioNames gives the absolute path name for each audio file in the DiscreteWavelets package."
DataList::usage="DataList gives information about data files that come with the DiscreteWavelets package."
DataNames::usage="DataNames gives the absolute path name for each data file in the DiscreteWavelets package."
DataList::usage="DataList gives information about data files that come with the DiscreteWavelets package."
DataNames::usage="DataNames gives the absolute path name for each data file in the DiscreteWavelets package."
MaximumIterations::usage="MaxIterations takes a numeric matrix or a list of three numeric matrices and returns the maximum number of iterations of the wavelet transformed that can be performed on the input."
ImageRead::usage="ImageRead takes a filename or a url and returns either a matrix or a list of three matrices that represent the digital image.  If the file does not exist, the routine returns a 2 x 2 zero matrix."
ImagePlot::usage="ImagePlot takes either a matrix or a list of three matrices and displays the image."
CreateImageObject::usage="CreateImageObject take either a matrix or a list of three matrices and returns a graphics object of the image that can be displayed using ImagePlot (or Show)."
LinMap::usage="LinMap scales a list so that its smallest value is 0 and its largest value is mx."
WaveletDensityPlot::usage="WaveletDensityPlot takes either a matrix or a list of three matrices (presumably transformed data) and displays the output."
wtLDP::usage="wtLDP is an auxiliary routine that takes a matrix a (assumed to be a wavelet transform) and uses LinMap on each lowpass/highpass part.  The routine uses LinMap with mx = 255.  The absolute value of the highpass portions are sent to LinMap."
ExtractRegion::usage="ExtractRegion takes either a matrix or a list of three matrices and returns a part of the data.  It also returns the indices of the upper left/lower right coordinates of the extracted region."
ExtractVectorPart::usage="ExtractVectorPart will extract a portion of the vector indicated by opts. It returns the extracted portion of the vector as well as the starting and stopping point."
WaveletVectorPlot::usage="WaveletVectorPart simply plots the transformed vector.  The user can decide to plot the parts they have indicated in different colors as well as add dividing lines."
WaveletVectorPlay::usage="WaveletVectorPlay will play a vector (presumably a wavelet transformation) as a sound file.  The user can decide which parts they want to play."
HuffmanTree::usage="HuffmanTree takes the Huffman codes generated by MakeHuffmanCodes (first argument) and makes a Huffman tree from them."


(* Usage for directives *)
Precision::usage="Default is $MachinePrecision.  Changing this directive to integer n means that NSolve works using Precision n."
NumIterations::usage="NumIterations is the number of iterations you want the transformation routine to perform.  The default is 1.  If set to Max, the transform will do as many iterations as the dimensions of the data will allow."
PowersOfTwo::usage="PowersOfTwo is an option that can be set in ChopVector or ReadImage.  The default value is None so that no data is removed from the input vector/matrix.  If PowersOfTwo is set to n, then the data is truncated so that n factors of 2 appear in the dimensions of the modified data."
PrintInfo::usage="PrintInfo is an option that can be set in ChopVector or ReadImage.  The default value is False.  If it is set to True, the routines print channel information (in the case of images) as well as length/dimension of the data."
DisplayMode::usage="DisplayMode is an option for RGBtoYCbCr and YCbCrtoRGB.  The default value is False in which case the routines do not stretch the channels in order to enhance the display. DisplayMode can be set to a 3x3 matrix as well."
Region::usage="Region is used to specify a particular part of transformed data.  The default value is All.  For 1D data, options are LowPass and HighPass while options for 2D data are Blur, Horizontal, Vertical, or Diagonal."
Iteration::usage="Iteration is a directive that instructs the routine to only retrieve/display the particular iteration.  The default is Max.  This means to display all iterations."
UseColors::usage="A directive for WaveletVectorPlot that indicates that different portions of the transformed data are to be different colors.  Default value is False."
ColorList::usage="ColorList is a directive for WaveletVectorPlot that lists the colors to be used if UseColors is set to True.  ColorList can be a list of elements of the form RGBColor[r,g,b], GrayLevel[h], Hue[h], etc."
DivideLines::usage="DivideLines is a directive for WaveletVectorPlot and WaveletDensityPlot.  Default value is False.  If set to True, the routine draws lines around different parts of the transformed data."
DivideLinesColor::usage="DivideLinesColor is a directive for WaveletVectorPlot and WaveletDensityPlot.  If DivideLines is set to True, then the lines around different parts of the transformed data are colored using this directive.  Default color is red.  Values can be entered using RGBColor, Hue, GrayLevel, or CYMKColor."
DivideLinesThickness::usage="DivideLinesThickness is a directive for WaveletVectorPlot and WaveletDensityPlot.  If DivideLines is set to True then the thickness of the lines around different parts of the transformed data are rendered using this thickness.  Default value is 0.002."
LinearScaling::usage="LinearScaling is a directive for ImagePlot.  Default value is False.  If it is set to true, the absolute value of the input data is linearly scaled so that the min maps to 0 and the max maps to 255.  The value can also be set to LeftWT or RightWT when plotting partial wavelet transforms of images."
ChannelColor::usage="ChannelColor is an alternative to ColorFunction in ImagePlot.  It can be used to change the color with which a grayscale image can be shaded.  The default value is Gray.  Values can be of the form RGBColor[r,g,b], GrayLevel[h], Hue[h], etc."
ImageSize::usage="ImageSize is a standard Graphics options, but in this package it has a new default value of FixEdges.  If ImageSize is set to FixEdges, ImagePlot utilizes Radka Turcajova's fix of the borders that are used in Raster."
Boundary::usage="Boundary is a directive used by biorthogonal wavelet transform modules.  The default value is None and in this case, the wavelet transform is computed normally.  If the directive is set to Reflective, then the module utilizes the symmetry of the biorthogonal filters to better deal with edge conditions.  See Section 11.3 for more details."
IntegerMap::usage="IntegerMap is a directive used by LeGall transform modules.  The default value is False.  If set to true, the module is modified so that the output is integer-valued."
ImageType::usage="ImageType is a directive used by ImageNames, ImageList, and ShowThumbnails to determine the type of image to process.  The default value is All and in this case, the modules process both GrayScale and Color images.  Other options are GrayScale and Color."
ListThumbnails::usage="ListThumbnails is a directive used by ImageNames.  The default value is False and in this case, ImageNames returns complete path and file names to the full-size version of the images.  If set to True, the routine provides complete paths and names for thumbnails of the images."
SlideShow::usage="SlideShow is a directive used by ShowThumbnails.  The default value is False and in this case, the thumbnails are displayed in a rectangular grid.  If the directive is set to True, the thumbnails are displayed in a slide show format."
NodeColor::usage="NodeColor is a directive used by HuffmanTree.  The default value is a dark red.  It can be set to any RGBColor, GrayLevel, Hue, etc. color desired."
NodeEdgeColor::usage="NodeEdgeColor is a directive used by HuffmanTree to set the color for node boundaries.  The default value is Black.  It can be set to any RGBColor, GrayLevel, Hue, etc. color desired."
NodeSize::usage="NodeSize is a directive used by HuffmanTree to set the radius of each node.  The default value is 0.2."
NodeFontSize::usage="NodeFontSize is a directive used by HuffmanTree to set the size of the font used to label each node.  The default value is 12."
NodeFontColor::usage="NodeFontColor is a directive used by HuffmanTree to set the color of the text used to label each node.  The default value is a dark gold.  It can be set to any RGBColor, GrayLevel, Hue, etc. color desired."
EdgeLength::usage="EdgeLength is a directive used by HuffmanTree.  It is a numerical value that is passed to the TreePlot\!\(\*
StyleBox[\" \",\nFontWeight->\"Bold\"]\)directive LayerSizeFunction to help control the proportion sizes of nodes and branches.  The function used by HuffmanTree is #^r, where r is the value of EdgeLength (the default value is 0.1)."
ShowString::usage="ShowString is a directive used by HuffmanTree.  It is set to False, but it can be set to a character string that serves as a label for the rendered Huffman tree."
PrintResult::usage="PrintResult is a directive used by TestSparseness.  The default value is False.  If set to True, the module reports whether or not the input was sparse and which threshold value was returned."


(* Options *)

Options[RGBToYCbCr] = Options[YCbCrToRGB] = {DisplayMode->False};
Options[Daub] = Options[Coif] = {Precision -> $MachinePrecision};
Options[SplineFilters]={PrintInfo->False};
Options[HWT1D] = Options[IHWT1D] = Options[WT1D] = Options[IWT1D] = {NumIterations->1};
Options[BWT1D1] = Options[IBWT1D1] = {Boundary->None};
Options[BWT1D] = Options[IBWT1D] = {NumIterations->1,Boundary->None};
Options[LWT1D1] = Options[ILWT1D1] = {IntegerMap->False};
Options[LWT1D] = Options[ILWT1D] = {NumIterations->1,IntegerMap->False};
Options[HWT2D] = Options[IHWT2D] = Options[WT2D] = Options[IWT2D] = {NumIterations->1};
Options[BWT2D1] = Options[IBWT2D1] = {Boundary->None};
Options[BWT2D] = Options[IBWT2D] = {NumIterations->1,Boundary->None};
Options[LWT2D1] = Options[ILWT2D1] = {IntegerMap->False};
Options[LWT2D] = Options[ILWT2D] = {NumIterations->1,IntegerMap->False};
Options[WaveletMatrixToList] = Options[WaveletListToMatrix] = {NumIterations->1};
Options[ChopVector]={PowersOfTwo->None,PrintInfo->False};
Options[WaveletVectorToList]={NumIterations->1};
Options[ImageList]={ImageType->All};
Options[ImageNames]={ImageType->All,ListThumbnails->False};
Options[ShowThumbnails]={ImageType->All,SlideShow->False};
Options[ImageRead]={PowersOfTwo->None,PrintInfo->False};
Options[ImagePlot]=Join[Join[Options[Graphics],Options[Raster]],{LinearScaling->False,ChannelColor->None}];
SetOptions[ImagePlot,AspectRatio->Automatic,ImageSize->FixEdges];
Options[WaveletDensityPlot]=Join[Options[ImagePlot],{Region->All,Iteration->All,NumIterations->1,DivideLines->True,DivideLinesColor->RGBColor[1,1,1],DivideLinesThickness->.002}];
Options[ExtractRegion]={NumIterations->1,Iteration->All,Region->All};
Options[ExtractVectorPart]={Region->All,Iteration->All,NumIterations->1};
Options[WaveletVectorPlot]=
    Join[Options[Graphics],{NumIterations->1,Iteration->All,
        Region->All,PointSize->0.004,UseColors->False,
        ColorList->{RGBColor[.25,0,.5],RGBColor[.75,0,.5],
            RGBColor[.5,.5,0],RGBColor[.25,.5,.5],RGBColor[.5,.25,.25],
            RGBColor[.5,.5,.75]},DivideLines->True,
        DivideLinesColor->RGBColor[1,0,0],
        DivideLinesThickness->.004}];
SetOptions[WaveletVectorPlot,Axes->True,AspectRatio->1/GoldenRatio];
Options[WaveletVectorPlay]=Join[Options[ListPlay],{NumIterations->1,Iteration->1,Region->All,SampleRate->11025}];
Options[HuffmanTree]={NodeColor->RGBColor[178./255,142./255,46./255],NodeEdgeColor->Black,NodeSize->.2,NodeFontSize->12,NodeFontColor->RGBColor[138./255,30./255,2./255],EdgeLength->.1,ShowString->""};
Options[TestSparseness]={PrintResult->False};
Options[UniversalThreshold]={NumIterations->1};
Options[WaveletShrinkage]={NumIterations->1};
Options[SureShrink]={NumIterations->1};


        
Begin["`Private`"]

Off[General::spell];
Off[General::spell1];



CE[v_] := Module[{w},
	If[Length[v]==0,Message[CE::"zerolength"];w={},
    w=Drop[FoldList[Plus, 0, Sort[Abs[Flatten[v]], Greater]^2], 1]/
      Norm[Flatten[v]]^2];
    Return[w];
      ];
      
Comp[v_,k_]:=Module[{cols,w,j,n},
	If[Length[v]==0,Message[Comp::"zerolength"];Return[{}]];
	If[!IntegerQ[k],Message[Comp::"noninteger"];Return[{}]];
	n=Which[ArrayDepth[v]==1,Length[v],True,Apply[Times,Dimensions[v]]];
	j=If[k>n,n,k];
	w=If[k<=0,Table[0,{Length[Flatten[v]]}],Chop[Flatten[N[v]], \
Sort[Flatten[Abs[v]],Greater][[j]]]];
			Which[ArrayDepth[v]==1,Return[w],True,cols=Dimensions[v][[2]];Return[\
Partition[w,cols]]]];

nCE[v_, p_]:= Module[{},
	If[ArrayDepth[v]!=1,Message[nCE::"nonvector"];Return[0]];
	Return[Length[Select[v, (# <= p) &]]]];

Entropy[v_]:=Module[{c,s,t,n},
				If[Length[v]==0,Message[Entropy::"nonvectormatrix"];Return[0]];
				n=If[ArrayDepth[v]==1,Length[v],Apply[Times,Dimensions[v]]];
				c=Map[Length,Split[Sort[Flatten[v]]]];
				s=c/n;t=Log[2,n/c];
    			Return[s.t]];
      
MSE[a_,b_]:=Module[{},
	If[(ArrayDepth[a] || \
ArrayDepth[b])!=2,Message[MSE::"nonmatrix"];Return[-1]];
	If[Dimensions[a]!=Dimensions[b],Message[MSE::"baddimensions"];Return[-1]];
	Return[Total[Flatten[(a-b)^2]]/Apply[Times,Dimensions[a]]];
	];

PSNR[a_, b_] := Module[{},
	If[MSE[a,b]==-1,Message[PSNR::"badmatrix"];Return[0]];
	If[MSE[a,b]==0,Return[Infinity]];
	Return[10.*Log[10, 255^2/MSE[a,b]]];
];  



RGBToYCbCr[c_,opts___]:=
    Module[{k,cols,t,d,
        a={{.299,.587,.114},{-.299,-.587,.886}/1.772,{.701,-.587,-.114}/1.402}/255.,
        r,g,b,v={219,224,224},w={16,128,128},m,dflag=0},
      m=DisplayMode/.{opts}/.Options[RGBToYCbCr];
      If[m===True,dflag=1];
      d=Depth[c];
      If[d==2,t=a.c;If[dflag==0,Return[t],Return[v*t+w]]];
      If[d==3,t=Transpose[a.Transpose[c]];If[dflag==0,Return[t],Return[\
Transpose[v*Transpose[t]+w]]]];
      If[d==4,{r,g,b}=c;cols=Dimensions[r][[2]];
        t=
          Partition[Flatten[Transpose[{Flatten[r],Flatten[g],Flatten[b]}]],
            3];
        t=Transpose[a.Transpose[t]];
        If[dflag==1,t=Transpose[v*Transpose[t]+w]];
        Return[
          Table[Partition[Take[Flatten[t],{k,3*Length[t],3}],cols],{k,1,
              3}]];];
      ];


YCbCrToRGB[c_,opts___]:=
    Module[{k,cols,t,d,a=Inverse[{{.299,.587,.114},{-.299,-.587,.886}/1.772,{.701,-.587,-.114}/1.402}/255.],
    r,g,b,m,dflag=0,v={219,224,224},w={16,128,128}},
      m=DisplayMode/.{opts}/.Options[YCbCrToRGB];
      If[m===True,dflag=1];
      d=Depth[c];
      If[d==2,
      	If[dflag==0,Return[a.c],
      		tmp=(c-w)/v;Return[a.tmp]]];
      If[d==3,
      	If[dflag==0,Return[Transpose[a.Transpose[c]]],
      		tmp=Transpose[(Transpose[c]-w)/v];Return[Transpose[a.Transpose[tmp]]]];
      Return[Transpose[a.Transpose[c]]];];
      If[d==4,{r,g,b}=c;cols=Dimensions[r][[2]];
        t=Partition[Flatten[Transpose[{Flatten[r],Flatten[g],Flatten[b]}]],3];
        If[dflag==0,t=Transpose[a.Transpose[t]],tmp=Transpose[(Transpose[t]-w)/v];t=Transpose[a.Transpose[tmp]]];
        Return[
          Table[Partition[Take[Flatten[t],{k,3*Length[t],3}],cols],{k,1,
              3}]];];
      ];


MakeHuffmanCodes[v_]:=
    Module[{data,alphabet,queue,getFreq,orderFunct,left,right,newnode,codes,cnt},
    If[Head[v]==List,
    	If[MemberQ[Union[Map[IntegerQ,Flatten[v]]],False]==True,Message[MakeHuffmanCodes::"badinput"];Return[{0,0,0}]];
    ];
    data=Characters[Switch[Head[v],String,v,_,FromCharacterCode[Flatten[v]]]];
      alphabet=Union[data];
      queue=Map[{#,Count[data,#]}&,alphabet];
      orderFunct[a_,b_]:=
        Switch[Head[a[[2]]],Integer,a[[2]]<b[[2]],String,
          StringLength[a[[2]]]<=StringLength[b[[2]]],_,
          Print["nope: ",Head[a[[2]]]];];
      queue=Sort[queue,orderFunct];
      While[Length[queue]>1,left=queue[[1]];
        right=queue[[2]];
        newnode={{left,right},left[[2]]+right[[2]]};
        queue=Sort[Append[Drop[queue,2],newnode],orderFunct];];
      getCodes[node_,codeString_]:=Module[{left,right},left=node[[1]];
          Switch[Head[left],List,right=left[[2]];left=left[[1]];
            
            Union[getCodes[left,codeString<>"0"],
              getCodes[right,codeString<>"1"]],Symbol,{{left,codeString}},
            String,{{left,codeString}}]];
      codes=Transpose[Sort[getCodes[queue[[1]],""],orderFunct]];
      codes=Append[codes,Map[ToCharacterCode,codes[[2]]]-48];
      cnt=Map[Count[data,#]&,codes[[1]]];
      codes[[2]]=N[cnt/Length[data]];
      If[Head[v]==List,
        codes[[1]]=Flatten[ToCharacterCode[codes[[1]]]]];
      Return[{Transpose[codes],Length[data]*8,cnt.Map[Length,codes[[3]]]}];
]; 


FiniteFourier[w_,v_,idx_:0]:=Module[{k},
	If[!IntegerQ[idx],Message[FiniteFourier::"badindex"];Return[0]];
	Return[v.Table[E^(I*k*w),{k,idx,idx+Length[v]-1}]];
	];
	
Cs[w_]:=(E^(I*w)+E^(-I*w))/2;

Sn[w_]:=(E^(I*w)-E^(-I*w))/(2*I);

Haar[___]:=Sqrt[2]*{1,1}/2;

Daub[n_Integer, opts___] := 
	Module[{prec, j, m, z, eqs, k, s, coeffs, rtwo, rthree, b, c, x, sls},  
    If[!IntegerQ[n],Message[Daub::"ninteger",n];Return[{0,0}]];
    If[(n<=0)||(Mod[n,2]!=0),Message[Daub::"nposeven",n];Return[{0,0}];];
    
   (* In the case where n = 2, 4, 6, go ahead and return the symbolic values \
of the filter.*)
   rtwo = Sqrt[2];
   If[n == 2, Return[{rtwo/2, rtwo/2}]];
   rthree = Sqrt[3];
   If[n == 4, Return[{1 + rthree, 3 + rthree, 3 - rthree, 1 - \
rthree}/(4*rtwo)]];
   b = 2*Sqrt[5]; c = Sqrt[5/2 + Sqrt[10]];
   If[n == 6, Return[{rtwo + b + 2*c, 5*rtwo + b + 6*c, 10*rtwo - 2*b + 4*c, \
10*rtwo - 2*b - 4*c, 5*rtwo + b - 6*c, rtwo + b - 2*c}/32]];

  (* Create an array of unknowns, build the system of equations and then \
solve via NSolve.  Note that the user can enter the precision used by NSolve. \
*)
  
  z = Array[h, n];
  eqs = Flatten[{Sum[z[[k]]^2, {k, 1, n}] == 1, 
  Table[Sum[z[[k + 1]]*z[[k + 1 - 2*m]], {k, 2*m, n - 1}] == 0, {m, 1, n/2 - \
1}], 
  Sum[z[[k]], {k, 1, n}] == Sqrt[2], 
  Sum[(-1)^(k + 1)*z[[k]], {k, 1, n}] == 0, 
  Table[Sum[(-1)^(k - 1)*(k - 1)^m*z[[k]], {k, 2, n}] == 0, {m, 1, n/2 - \
1}]}];
  
  prec = Precision /. {opts} /. Options[Daub];
  s = NSolve[eqs, z, prec];
  (* Now extract the coefficients and only take those with real values. *)
  
  coeffs = Table[z /. s[[k]], {k, 1, Length[s]}];
  coeffs = Table[Select[coeffs[[k]], (Conjugate[#] == #) &], {k, 1, \
Length[coeffs]}];
  coeffs = Select[coeffs, Length[#] > 0 &];
  
  (* Finally, form the trig polynomial using each solution and loop through \
the roots of the polynomial until we find the one 
   whose maximum value (in modulus) is 1.  Note that the trig polynomial uses \
a positive k so we need to reverse the coefficients we return.  *)
   
  For[m = 1, m <= Floor[n/2], m++,
        s = NSolve[Sum[coeffs[[m, k + 1]]*x^k, {k, 0, n - 1}] == 0, x];
  		p = Max[Abs[Table[x /. s[[k]], {k, 1, Length[s]}]]];
  		If[p == 1.0, Return[Reverse[coeffs[[m]]]]];
        ];
  Return[m];
  ];
  
Coif[K_,opts___]:=Module[{a0=(Sqrt[7]-1)/2,a1=(3-Sqrt[7])/2,h,z,H,k,m,deriv0,\
derivPi,orth,eqs,prec,s,sln},
  	If[!IntegerQ[K],Message[Coif::"noninteger",K];Return[{0,0}]];
  	If[(K<1) || (K>3),Message[Coif::"badinteger",K];Return[{0,0}]];
  	If[K==1,Return[Sqrt[2]*{-a0,4-a1,8+2*a0,4+2*a1,-a0,-a1}/16]];
  	z=Array[h,6*K,-2*K];
	H[w_]:=z.Table[E^(I*k*w),{k,-2*K,4*K-1}];
	deriv0=Table[Derivative[m][H][0]==Sqrt[2]*DiscreteDelta[m],{m,0,2*K-1}];
	derivPi=Table[Derivative[m][H][Pi]==0,{m,0,2*K-1}];
	orth=Table[Sum[h[k]*h[k+2m],{k,-2*K,4*K-1-2*m}]==DiscreteDelta[m],{m,0,3*K-1}\
];
	eqs=Flatten[{deriv0,derivPi,orth}];
	prec = Precision /. {opts} /. Options[Coif];
	s=NSolve[eqs,z,prec];
	sln=ReplaceAll[z,s];
	If[K==2,Return[First[sln]]];
	If[K==3,Return[sln[[4]]]];
  ];
  
SplineFilters[M_,Mw_,opts___]:=
  Module[{n,nw,t,w,m,idx,init,term,f,fhat,ph,y,z},
  If[Head[M]==Rule || Head[Mw]==Rule,Return[{Haar[],Haar[]}];];
  If[!IntegerQ[M] || \
!IntegerQ[Mw],Message[SplineFilters::"badintegers"];Return[{{0,0},{0,0}}];];
  If[M<0 || Mw <0, \
Message[SplineFilters::"badintegers"];Return[{{0,0},{0,0}}];];
  {n,nw}=({M,Mw}-Mod[M,2])/2;
  If[OddQ[Abs[M-Mw]], \
Message[SplineFilters::"badintegers"];Return[{{0,0},{0,0}}];];
  idx=Mod[M+1,2];
  {init,term}={-2*n-nw+idx,2*n+nw+(-1)^Mod[M+1,2]};
  t[w_]=E^(I*w*Mod[M,2]/2)*Cs[w/2]^M*
      Sum[Binomial[n+nw-idx+m,m]*Sn[w/2]^(2*m),{m,0,n+nw-idx}];
  f=Sqrt[2]*Table[Coefficient[Expand[t[w]],E^(I*w),k],{k,init,term}];
  fhat=Sqrt[2]*Table[Binomial[Mw,m],{m,0,Mw}]/2^Mw;
  ph=PrintInfo/.{opts}/.Options[SplineFilters];
  Clear[y,z];
  y=Table["h["<>ToString[k]<>"]",{k,init,term}];
  z=Table["ht["<>ToString[k]<>"]",{k,-Mw/2+Mod[Mw,2]/2,Mw/2+Mod[Mw,2]/2}];
    If[ph===True,
        Print["Spline Filter Information: \n\nLength of the first filter = ",
          2*M+Mw-1," and it is of the form ",y,
          ".\n\nLength of the second filter = ",Mw+1,
          " and it is of the form ",z,".\n"]];
    Return[{f,fhat}];
    ];
    
    SplineFilters[M_]:={Haar[],Haar[]};

CDF97[___]:=Module[{t,rts,r,k,p,pw,cs,sn,H,Hw,a,y,z},
      	rts=NSolve[1+4*t+10*t^2+20*t^3==0,t];
      	r=Table[t/.rts[[k]],{k,1,3}];
      	p[t_]:=a*(t-r[[1]]);
      	pw[t_]:=20*(t-r[[2]])*(t-r[[3]])/a;
      	H[w_]:=Sqrt[2]*Cs[w/2]^4*p[Sn[w/2]^2];
      	Hw[w_]:=Sqrt[2]*Cs[w/2]^4*pw[Sn[w/2]^2];
      	Clear[a];
      	a=a/.Flatten[NSolve[H[0]==Sqrt[2],a]];
      	y=Chop[Table[Coefficient[Expand[H[w]],E^(I*w),k],{k,-3,3}]];
      	z=Chop[Table[Coefficient[Expand[Hw[w]],E^(I*w),k],{k,-4,4}]];
      	Return[{z,y}];
      ];

LeGall[___]:=Reverse[SplineFilters[2,2]]*{Sqrt[2],Sqrt[2]/2}
  
HWT1D1[v_]:=Module[{x},
  	If[!EvenQ[Length[v]],Message[HWT1D1::"evenlength"];Return[v]];
  	x=Partition[v,2];
  	Return[Sqrt[2]*Join[x.{1,1},x.{-1,1}]/2];
  	];
  	
IHWT1D1[v_]:=Module[{x,y},
  	If[!EvenQ[Length[v]],Message[IHWT1D1::"evenlength"];Return[v]];
  	x=Transpose[Partition[v,Length[v]/2]];
  	y=Sqrt[2]*Flatten[Transpose[{x.{1,-1},x.{1,1}}]]/2;
  	Return[y];
  	];
  	
HWT1D[v_, opts___]:=Module[{n,hp,a,b,j,maxits,its},
      (* Make sure the length of the vector and figure out how many \
iterations are desired.*)
      n=Length[v];
      
      If[ArrayDepth[v]!=1,Message[HWT1D::"nonvector"];Return[v];];
      If[Mod[n,2]==1,Message[HWT1D::"oddlength"];Return[v];];
      
      maxits=FactorInteger[n][[1,2]];
      its=NumIterations/.{opts}/.Options[HWT1D];
      If[!IntegerQ[its],Message[HWT1D::"badits",maxits];its=1];
      If[its<0,Message[HWT1D::"badits",maxits];its=1];
      
      If[its>maxits,Message[HWT1D::"maxits",its,maxits];its=maxits];
      	
      (* hp will hold the highpass portions - we append it as we loop through \
the iterations.  a holds the current lowpass portion.  
      We prepend the final a to hp and return this. *)
    
      hp={};
      a=v;
      For[j=1,j<=its,j++,
        b=HWT1D1[a];
        a=Take[b,n/(2^j)];
        hp=Join[Drop[b,n/(2^j)],hp];
        ];
      Return[Join[a,hp]];
      ];
      
IHWT1D[v_,opts___]:=Module[{n,a,j,b,z,s,its},
      
(* Make sure the length of the vector and figure out how many iterations are \
desired.*)
      
	n=Length[v];
      
    If[ArrayDepth[v]!=1,Message[IHWT1D::"nonvector"];Return[v];];
    If[Mod[n,2]==1,Message[IHWT1D::"oddlength"];Return[v];];
      
    maxits=FactorInteger[n][[1,2]];
    its=NumIterations/.{opts}/.Options[IHWT1D];
    If[!IntegerQ[its],Message[IHWT1D::"badits",maxits];its=1];
    If[its<0,Message[IHWT1D::"badits",maxits];its=1];
      
    If[its>maxits,Message[IHWT1D::"maxits",its,maxits];its=maxits];

      
    (* a starts out holding the transformed vector.  At each piece we replace \
the smallest lowpass highpass tandem with the next level's lowpass portion.  
            When we finish, a holds the inversely transformed data. *)
      
      a=v;
      For[j=1,j<=its,j++,
        y=Take[a,n/(2^(its-j))];
        b=IHWT1D1[y];
        z=Drop[a,n/(2^(its-j))];
        a=Join[b,z];
        ];
      Return[a];
      ];
      
WT1D1[v_,h_]:=Module[{hh=Reverse[h],g,x},
		If[OddQ[Length[h]],Message[WT1D1::"badfilterlength"];Return[v]];
		If[!EvenQ[Length[v]],Message[WT1D1::"evenlength"];Return[v]];
      	g=h*Table[(-1)^k,{k,1,Length[h]}];
      	x=Partition[v,Length[h],2,{1,2}];
      	Return[Join[x.hh,x.g]];
      	];
      	
WT1D1[v_]:=HWT1D1[v];
      
IWTht[v_,h_]:=Module[{eh,oh,z,L},
      	L=Length[h];
	    If[L==2,Return[Flatten[Transpose[{h[[1]]*v,h[[2]]*v}]]],
	    	{oh,eh}={Take[Reverse[h],{1,L,2}],Take[Reverse[h],{2,L,2}]};
        	z=Partition[Join[Take[v,-(L/2-1)],v],L/2,1];
        	Return[Flatten[Transpose[{z.eh,z.oh}]]];
        ];
      ];
      
IWT1D1[v_,h_]:=Module[{hh=Reverse[h],t,b,g},
		If[OddQ[Length[h]],Message[IWT1D1::"badfilterlength"];Return[v]];
      	If[!EvenQ[Length[v]],Message[IWT1D1::"evenlength"];Return[v]];
      	{t,b}=Partition[v,Length[v]/2];
	    g=Table[(-1)^k,{k,1,Length[h]}]*h;
	    Return[IWTht[t,hh]+IWTht[b,g]];
	  ];
	  
IWT1D1[v_]:=IHWT1D1[v];
	  
WT1D[v_,h_,opts___]:=Module[{n,hp,a,b,j,maxits,its},
      (* Make sure the length of the vector and figure out how many \
iterations are desired.*)
      
      If[Head[h]==Rule,Return[HWT1D[v,h]]];
      n=Length[v];
      If[ArrayDepth[v]!=1,Message[WT1D::"nonvector"];Return[v];];
      If[Mod[n,2]==1,Message[WT1D::"oddlength"];Return[v];];
      
      maxits=FactorInteger[n][[1,2]];
      its=NumIterations/.{opts}/.Options[WT1D];
      If[!IntegerQ[its],Message[WT1D::"badits",maxits];its=1];
      If[its<0,Message[WT1D::"badits",maxits];its=1];
      
      If[its>maxits,Message[WT1D::"maxits",its,maxits];its=maxits];
      	
      (* hp will hold the highpass portions - 
          we append it as we loop through the iterations.  
              a holds the current lowpass portion.  
              We prepend the final a to hp and return this. *)
    
      hp={};
      a=v;
      For[j=1,j<=its,j++,
        b=WT1D1[a,h];
        a=Take[b,n/(2^j)];
        hp=Join[Drop[b,n/(2^j)],hp];
        ];
      Return[Join[a,hp]];
      ];
      
IWT1D[v_,h_,opts___]:=Module[{n,a,j,b,z,s,its},
      
      If[Head[h]==Rule,Return[IHWT1D[v,h]]];
      n=Length[v];
      If[ArrayDepth[v]!=1,Message[IWT1D::"nonvector"];Return[v];];
      If[Mod[n,2]==1,Message[IWT1D::"oddlength"];Return[v];];
      
      maxits=FactorInteger[n][[1,2]];
      its=NumIterations/.{opts}/.Options[IWT1D];
      If[!IntegerQ[its],Message[IWT1D::"badits",maxits];its=1];
      If[its<0,Message[IWT1D::"badits",maxits];its=1];
      
      If[its>maxits,Message[IWT1D::"maxits",its,maxits];its=maxits];
      
      (* a starts out holding the transformed vector.  
            At each piece we replace the smallest lowpass/
            highpass tandem with the next level's lowpass portion.  
            When we finish, a holds the inversely transformed data. *)
      
      a=v;
      For[j=1,j<=its,j++,
        y=Take[a,n/(2^(its-j))];
        b=IWT1D1[y,h];
        z=Drop[a,n/(2^(its-j))];
        a=Join[b,z];
        ];
      Return[a];
      ];
     
BWT1D1[v_, \
opts___]:=Module[{},Message[BWT1D1::"nofilter"];Return[HWT1D1[v]]];

BWT1D1[v_, {h__, hw___}, \
opts___]:=Module[{b,Lh,Lhw,Endh,Endhw,p,gw,xtop,xbot,vv,r,s},
	  If[ArrayDepth[v]!=1,Message[BWT1D1::"nonvector"];Return[v];];
      If[Mod[Length[v],2]==1,Message[BWT1D1::"oddlength"];Return[v];];
      If[{Length[{h}],Length[{hw}]}!={1,1} || \
Length[Flatten[{h,hw}]]==2,Message[BWT1D1::"onefilter"];Return[WT1D1[v,\
Flatten[{h,hw}]]]];
      
      
      {Lh,Lhw}=Map[Length,{h,hw}];
      If[OddQ[Abs[Lh-Lhw]],Message[BWT1D1::"badfilterlengths"];Return[v];];
      
      b=Boundary/.{opts}/.Options[BWT1D1];
      
      vv=If[ToString[b]==="Reflective",
          {r,s}={1,Length[v]}+Mod[Lh,2]*{1,-1};Join[v,Take[Reverse[v],{r,s}]],\

          v];
      
      {Endh,Endhw}=Map[((#-Mod[#,2])/2)&,{Lh,Lhw}];
      p={0,Lh-1}+Mod[Endh+1,2];
      gw=Reverse[h]*Table[(-1)^k,{k,First[p],Last[p]}];
      
      xtop=Partition[RotateRight[vv,Endhw-Mod[Lh+1,2]],Lhw,2,{1,2}];
      
      xbot=Partition[RotateRight[vv,Endh-1],Lh,2,{1,2}];
      
      {xtop,xbot}=
        If[ToString[b]==="Reflective",{Take[xtop,Length[v]/2],
            Take[xbot,Length[v]/2]},{xtop,xbot}];
      Return[Join[xtop.hw,xbot.gw]];
      ];
      
IBWTht[v_,f_,i_]:=Module[{Lf,Stopf,oddf,evenf,p,RotateEven,RotateOdd,vEven,\
vOdd,t},
    Lf=Length[f];
    Stopf=(Lf-Mod[Lf,2])/2;
      
    Which[OddQ[Lf],p=If[OddQ[Stopf],{1,2},{2,1}],True,p=If[OddQ[Stopf],{2,1},{\
1,2}]];
    {oddf,evenf}={Take[f,{First[p],Lf,2}],Take[f,{Last[p],Lf,2}]};
      
      
    (* This If statement handles the case of the highpass odd length filter \
*)
    If[(i==1 && OddQ[Lf]),{oddf,evenf}={evenf,oddf};
    Stopf=Stopf+1];
      
    RotateEven=Floor[Stopf/2];
    RotateOdd=Floor[(Stopf-1)/2];
      
    vEven=Partition[RotateRight[v,RotateEven],Length[evenf],1,{1,1}];
      
    vOdd=Partition[RotateRight[v,RotateOdd],Length[oddf],1,{1,1}];
      
    t={vEven.Reverse[evenf],vOdd.Reverse[oddf]};
    Return[Flatten[Transpose[t]]];
];


IBWT1D1[v_, opts___]:=Module[{},Message[IBWT1D1::"nofilter"];Return[IHWT1D1[v]]];

IBWT1D1[v_, {h__,hw___}, opts___]:=Module[{b,lp,hp,n,vv,Lhw,hwStop,p,k,g,y},
	If[ArrayDepth[v]!=1,Message[BWT1D1::"nonvector"];Return[v];];
    If[Mod[Length[v],2]==1,Message[BWT1D1::"oddlength"];Return[v];];
    If[{Length[{h}],Length[{hw}]}!={1,1} || \
Length[Flatten[{h,hw}]]==2,Message[IBWT1D1::"onefilter"];Return[IWT1D1[v,\
Flatten[{h,hw}]]]];

    b=Boundary/.{opts}/.Options[IBWT1D1];
      
    Lhw=Length[hw];
    vv=If[ToString[b]==="Reflective",{lp,hp}=Partition[v,Length[v]/2];
  		lp=Join[lp,Drop[Reverse[lp],-Mod[Lhw,2]]];
       	hp=Join[hp,(-1)^(Lhw+1)*Drop[Reverse[hp],Mod[Lhw,2]]];
       	Join[lp,hp],v];
      
    hwStop=(Lhw-Mod[Lhw,2])/2;
    p={0,Lhw-1}+Mod[hwStop+1,2];
    g=hw*Table[(-1)^k,{k,First[p],Last[p]}];
    y=IBWTht[Take[vv,Length[vv]/2],h,0]+IBWTht[Drop[vv,Length[vv]/2],g,1];
    If[ToString[b]==="Reflective",y=Take[y,Length[v]]];
    Return[y];
];

BWT1D[v_, {h__,hw___}, opts___]:=Module[{n,hp,a,b,j,maxits,its},
(*Make sure the length of the vector and figure out how many iterations are \
desired.*)

	n=Length[v];
	If[ArrayDepth[v]!=1,Message[BWT1D::"nonvector"];Return[v];];
    If[Mod[n,2]==1,Message[BWT1D::"oddlength"];Return[v];];
    If[{Length[{h}],Length[{hw}]}!={1,1} || \
Length[Flatten[{h,hw}]]==2,Message[BWT1D::"onefilter"];Return[WT1D[v,Flatten[{\
h,hw}],opts]]];
      
    maxits=FactorInteger[n][[1,2]];
    its=NumIterations/.{opts}/.Options[BWT1D];
    If[!IntegerQ[its],Message[BWT1D::"badits",maxits];its=1];
    If[its<0,Message[BWT1D::"badits",maxits];its=1];
      
    If[its>maxits,Message[BWT1D::"maxits",its,maxits];its=maxits];

    (*hp will hold the highpass portions-we append it as we loop through the \
iterations.a holds the current lowpass portion.
      		We prepend the final a to hp and return this.*)
      		
    hp={};
    a=v;
    For[j=1,j<=its,j++,b=BWT1D1[a,{h,hw},opts];
    	a=Take[b,n/(2^j)];
        hp=Join[Drop[b,n/(2^j)],hp];
    ];
    Return[Join[a,hp]];
];

BWT1D[v_, \
opts___]:=Module[{},Message[BWT1D::"nofilter"];Return[HWT1D[v,opts]]];

IBWT1D[v_,{h__,hw___},opts___]:=Module[{n,a,j,b,z,s,its},

	n=Length[v];
	If[ArrayDepth[v]!=1,Message[IBWT1D::"nonvector"];Return[v];];
    If[Mod[n,2]==1,Message[IBWT1D::"oddlength"];Return[v];];
    If[{Length[{h}],Length[{hw}]}!={1,1} || \
Length[Flatten[{h,hw}]]==2,Message[IBWT1D::"onefilter"];Return[IWT1D[v,\
Flatten[{h,hw}],opts]]];
      
    maxits=FactorInteger[n][[1,2]];
    its=NumIterations/.{opts}/.Options[IBWT1D];
    If[!IntegerQ[its],Message[IBWT1D::"badits",maxits];its=1];
    If[its<0,Message[IBWT1D::"badits",maxits];its=1];
      
    If[its>maxits,Message[IBWT1D::"maxits",its,maxits];its=maxits];

    (*a starts out holding the transformed vector.At each piece we replace \
the smallest lowpass
            highpass tandem with the next level's lowpass portion.When we \
finish, a holds the inversely transformed data.*)
    a=v;
    For[j=1,j<=its,j++,y=Take[a,n/(2^(its-j))];
    	b=IBWT1D1[y,{h,hw},opts];
        z=Drop[a,n/(2^(its-j))];
        a=Join[b,z];
    ];
    Return[a];
];

IBWT1D[v_, \
opts___]:=Module[{},Message[IBWT1D::"nofilter"];Return[IHWT1D[v,opts]]];

LWT1D1[v_,opts___]:=Module[{e,o,s,d,f,t},
	If[ArrayDepth[v]!=1,Message[LWT1D1::"nonvector"];Return[v];];	
	If[!EvenQ[Length[v]],Message[LWT1D1::"evenlength"];Return[v]];
  	f=IntegerMap/.{opts}/.Options[LWT1D1];
  	If[ToString[f]=="True" && \
MemberQ[Map[IntegerQ,v],False],Message[LWT1D1::"noninteger"];Return[v]];
   	e=Take[v,{2,Length[v],2}];
   	o=Take[v,{1,Length[v],2}];
   	o=Append[o,Last[o]];
   	t=Map[Total,Partition[o,2,1]]/2;
   	d=e-If[f===True,Floor[t],t];
   	d=Prepend[d,First[d]];
   	t=Map[Total,Partition[d,2,1]]/4;
   	s=Drop[o,-1]+If[f===True,Floor[t+1/2],t];
   	Return[Join[s,Drop[d,1]]];
];


ILWT1D1[v_,opts___]:=Module[{s,d,e,o,f,t},
	If[ArrayDepth[v]!=1,Message[LWT1D1::"nonvector"];Return[v];];	
	If[!EvenQ[Length[v]],Message[ILWT1D1::"evenlength"];Return[v]];
    f=IntegerMap/.{opts}/.Options[ILWT1D1];
    If[ToString[f]=="True" && \
MemberQ[Map[IntegerQ,v],False],Message[ILWT1D1::"noninteger"];Return[v]];
    {s,d}=Partition[v,Length[v]/2];
    d=Prepend[d,First[d]];
    t=Map[Total,Partition[d,2,1]]/4;
    o=s-If[f===True,Floor[t+1/2],t];
    o=Append[o,Last[o]];
    t=Map[Total,Partition[o,2,1]]/2;
    e=Drop[d,1]+If[f===True,Floor[t],t];
    Return[Flatten[Transpose[{Drop[o,-1],e}]]];
];     

LWT1D[v_, opts___]:=Module[{n,hp,a,b,j,maxits,its},
	(* Make sure the length of the vector and figure out how many iterations are \
desired.*)
      
    n=Length[v];
    If[ArrayDepth[v]!=1,Message[LWT1D::"nonvector"];Return[v];];
    If[Mod[n,2]==1,Message[LWT1D::"oddlength"];Return[v];];
      
    maxits=FactorInteger[n][[1,2]];
    its=NumIterations/.{opts}/.Options[LWT1D];
    If[!IntegerQ[its],Message[LWT1D::"badits",maxits];its=1];
    If[its<0,Message[LWT1D::"badits",maxits];its=1];
      
    If[its>maxits,Message[LWT1D::"maxits",its,maxits];its=maxits];
      	
    (* hp will hold the highpass portions - we append it as we loop through \
the iterations.  a holds the current lowpass portion.  
    	We prepend the final a to hp and return this. *)
    
     hp={};
     a=v;
     For[j=1,j<=its,j++,
     	b=LWT1D1[a,opts];
        a=Take[b,n/(2^j)];
        hp=Join[Drop[b,n/(2^j)],hp];
     ];
     Return[Join[a,hp]];
];

ILWT1D[v_,opts___]:=Module[{n,a,j,b,z,s,its},
      
	(* Make sure the length of the vector and figure out how many iterations are \
desired.*)
      
	n=Length[v];
      
    If[ArrayDepth[v]!=1,Message[ILWT1D::"nonvector"];Return[v];];
    If[Mod[n,2]==1,Message[ILWT1D::"oddlength"];Return[v];];
      
    maxits=FactorInteger[n][[1,2]];
    its=NumIterations/.{opts}/.Options[ILWT1D];
    If[!IntegerQ[its],Message[ILWT1D::"badits",maxits];its=1];
    If[its<0,Message[ILWT1D::"badits",maxits];its=1];
      
    If[its>maxits,Message[ILWT1D::"maxits",its,maxits];its=maxits];

      
    (* a starts out holding the transformed vector.  At each piece we replace \
the smallest lowpass highpass tandem with the next level's lowpass portion.  
            When we finish, a holds the inversely transformed data. *)
      
      a=v;
      For[j=1,j<=its,j++,
        y=Take[a,n/(2^(its-j))];
        b=ILWT1D1[y,opts];
        z=Drop[a,n/(2^(its-j))];
        a=Join[b,z];
        ];
      Return[a];
];

GetCorner[m_,r_,c_]:=Module[{rows,cols},

	If[MatrixQ[m,NumericQ]==False,Message[GetCorner::"badinput"];Return[m];];
	If[Map[IntegerQ,{r,c}]!={True,True},Message[GetCorner::"integerdimensions"];\
Return[m];];
	{rows,cols}=Dimensions[m];
	If[r>rows || \
c>cols,Message[GetCorner::"baddimensions",r,c,rows,cols];Return[m];];
	
	Return[Take[m,r,c]];
];

PutCorner[m_,b_]:=Module[{mat,br,bc,mr,mc,mt,bt,tp,i},
	If[MatrixQ[m,NumericQ]==False,Message[PutCorner::"badinput"];Return[m];];
	If[MatrixQ[b,NumericQ]==False,Message[PutCorner::"badinput"];Return[b];];
	{br,bc}=Dimensions[b];
	{mr,mc}=Dimensions[m];
	If[Map[IntegerQ,{br,bc}]!={True,True},Message[PutCorner::"integerdimensions"];Return[m];];
	If[Map[IntegerQ,{mr,mc}]!={True,True},Message[PutCorner::"integerdimensions"];Return[m];];
	
	If[{br,bc}=={mr,mc},Return[b]];
	If[br==mr,Return[Join[b,Take[m,mr,bc-mc],2]];];
	If[bc==mc,Return[Join[b,Take[m,br-mr]]];];
	
	tp =Join[Join[b,Take[m,br,bc-mc],2],Take[m,br-mr]];
	Return[tp];
];


LeftHWT[a_]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[LeftHWT::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[LeftHWT::"baddimensions"];Return[a];];
	
	Return[Transpose[Map[HWT1D1,Transpose[a]]]];
];

RightHWT[a_]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[RightHWT::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[RightHWT::"baddimensions"];Return[a];];
	
	Return[Map[HWT1D1,a]];
];


HWT2D1[a_]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[HWT2D1::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[HWT2D1::"baddimensions"];Return[a];];
	
	Return[Transpose[LeftHWT[Transpose[LeftHWT[a]]]]];
];
	
LeftIHWT[a_]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[LeftIHWT::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[LeftIHWT::"baddimensions"];Return[a];];
	
	Return[Transpose[Map[IHWT1D1,Transpose[a]]]];
];

IHWT2D1[a_]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[IHWT2D1::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[IHWT2D1::"baddimensions"];Return[a];];
	
	Return[Transpose[LeftIHWT[Transpose[LeftIHWT[a]]]]];
];

HWT2D[a_,opts___]:=Module[{rows,cols,its,maxits,b,j,d,z},

	If[MatrixQ[a,NumericQ]==False,Message[HWT2D::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[HWT2D::"baddimensions"];Return[a];];
	{rows,cols}=Dimensions[a];

	(*Compute the maximum number of iterations.*)

	maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[HWT2D];
	If[!IntegerQ[its],Message[HWT2D::"badits",maxits];its=1];
	If[its<0,Message[HWT2D::"badits",maxits];its=1];
	If[its>maxits,Message[HWT2D::"maxits",its,maxits];its=maxits];
	
	b=a;
	
    For[j=0,j<its,j++,
    	d=GetCorner[b,rows/(2^j),cols/(2^j)];
        z=HWT2D1[d];
        b=PutCorner[b,z];
    ];
    Return[b];
];

IHWT2D[a_,opts___]:=Module[{rows,cols,its,maxits,b,j,d,z},

	If[MatrixQ[a,NumericQ]==False,Message[HWT2D::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[HWT2D::"baddimensions"];Return[a];];
	{rows,cols}=Dimensions[a];

	(*Compute the maximum number of iterations.*)

	maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[HWT2D];
	If[!IntegerQ[its],Message[HWT2D::"badits",maxits];its=1];
	If[its<0,Message[HWT2D::"badits",maxits];its=1];
	If[its>maxits,Message[HWT2D::"maxits",its,maxits];its=maxits];
	
	b=a;
      For[j=1,j<=its,j++,
        d=GetCorner[b,rows/(2^(its-j)),cols/(2^(its-j))];
        z=IHWT2D1[d];
        b=PutCorner[b,z];];
      Return[b];
];



LeftWT[a_,h___]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[LeftWT::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[LeftWT::"baddimensions"];Return[a];];
	If[VectorQ[h,NumericQ]==False,Message[LeftWT::"nofilter"];Return[LeftHWT[N[a]]]];
	If[OddQ[Length[h]],Message[LeftWT::"badfilterlength"];Return[a]];

	Return[Transpose[Map[WT1D1[#,h]&,Transpose[a]]]];
];

RightWT[a_,h___]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[RightHWT::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[RightHWT::"baddimensions"];Return[a];];
	If[VectorQ[h,NumericQ]==False,Message[RightWT::"nofilter"];Return[RightHWT[N[a]]]];
	If[OddQ[Length[h]],Message[RightWT::"badfilterlength"];Return[a]];
	
	Return[Map[WT1D1[#,h]&,a]];
];

WT2D1[a_,h___]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[WT2D1::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[WT2D1::"baddimensions"];Return[a];];
	If[VectorQ[h,NumericQ]==False,Message[WT2D1::"nofilter"];Return[HWT2D1[N[a]]]];
	If[OddQ[Length[h]],Message[RightWT::"badfilterlength"];Return[a]];
	
	Return[Transpose[LeftWT[Transpose[LeftWT[a,h]],h]]];
];

LeftIWT[a_,h___]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[LeftIWT::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[LeftIWT::"baddimensions"];Return[a];];
	If[VectorQ[h,NumericQ]==False,Message[LeftIWT::"nofilter"];Return[LeftIHWT[N[a]]]];
	If[OddQ[Length[h]],Message[LeftIWT::"badfilterlength"];Return[a]];
	
	Return[Transpose[Map[IWT1D1[#,h]&,Transpose[a]]]];
];

IWT2D1[a_,h___]:=Module[{},

	If[MatrixQ[a,NumericQ]==False,Message[IWT2D1::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[IWT2D1::"baddimensions"];Return[a];];
	If[VectorQ[h,NumericQ]==False,Message[IWT2D1::"nofilter"];Return[IHWT2D1[N[a]]]];
	If[OddQ[Length[h]],Message[IWT2D1::"badfilterlength"];Return[a]];
	
	Return[Transpose[LeftIWT[Transpose[LeftIWT[a,h]],h]]];
];

WT2D[a_,h__,opts___]:=Module[{rows,cols,its,maxits,b,j,d,z},

	If[MatrixQ[a,NumericQ]==False,Message[WT2D::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[WT2D::"baddimensions"];Return[a];];
	If[VectorQ[h,NumericQ]==False,Message[WT2D::"nofilter"];Return[HWT2D[N[a],opts]]];
	If[OddQ[Length[h]],Message[WT2D::"badfilterlength"];Return[a]];
	{rows,cols}=Dimensions[a];

	(*Compute the maximum number of iterations.*)

	maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[WT2D];
	If[!IntegerQ[its],Message[WT2D::"badits",maxits];its=1];
	If[its<0,Message[WT2D::"badits",maxits];its=1];
	If[its>maxits,Message[WT2D::"maxits",its,maxits];its=maxits];
	
	b=a;
	
    For[j=0,j<its,j++,
    	d=GetCorner[b,rows/(2^j),cols/(2^j)];
        z=WT2D1[d,h];
        b=PutCorner[b,z];
    ];
    Return[b];
	
];

IWT2D[a_,h__,opts___]:=Module[{rows,cols,its,maxits,b,j,d,z},

	If[MatrixQ[a,NumericQ]==False,Message[IWT2D::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[IWT2D::"baddimensions"];Return[a];];
	If[VectorQ[h,NumericQ]==False,Message[IWT2D::"nofilter"];Return[IHWT2D[N[a],opts]]];
	If[OddQ[Length[h]],Message[IWT2D::"badfilterlength"];Return[a]];
	{rows,cols}=Dimensions[a];

	(*Compute the maximum number of iterations.*)

	maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[IWT2D];
	If[!IntegerQ[its],Message[IWT2D::"badits",maxits];its=1];
	If[its<0,Message[IWT2D::"badits",maxits];its=1];
	If[its>maxits,Message[IWT2D::"maxits",its,maxits];its=maxits];
	
	b=a;
    For[j=1,j<=its,j++,
    	d=GetCorner[b,rows/(2^(its-j)),cols/(2^(its-j))];
        z=IWT2D1[d,h];
        b=PutCorner[b,z];
    ];

    Return[b];
];


BWT2D1[a_, opts___]:=Module[{},Message[BWT2D1::"nofilter"];Return[HWT2D1[N[a]]]];

BWT2D1[a_,{h__,hw___},opts___]:=Module[{b,c},
	If[MatrixQ[a,NumericQ]==False,Message[BWT2D1::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[BWT2D1::"baddimensions"];Return[a];];
    If[{Length[{h}],Length[{hw}]}!={1,1} || Length[Flatten[{h,hw}]]==2,Message[BWT1D1::"onefilter"];Return[WT2D1[a,Flatten[{h,hw}]]]];

    b=Transpose[Map[BWT1D1[#,{h,hw}]&,Transpose[a]]];
    c=Map[BWT1D1[#,{h,hw},opts]&,b];
    Return[c];
];

IBWT2D1[a_, opts___]:=Module[{},Message[IBWT2D1::"nofilter"];Return[IHWT2D1[N[a]]]];

IBWT2D1[a_,{h__,hw___},opts___]:=Module[{b,c},
	If[MatrixQ[a,NumericQ]==False,Message[IBWT2D1::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[IBWT2D1::"baddimensions"]; Return[a];];
    If[{Length[{h}],Length[{hw}]}!={1,1} || Length[Flatten[{h,hw}]]==2,Message[IBWT2D1::"onefilter"];Return[IWT2D1[a,Flatten[{h,hw}]]]];

    b=Transpose[Map[IBWT1D1[#,{h,hw}]&,Transpose[a]]];
    c=Map[IBWT1D1[#,{h,hw},opts]&,b];
    Return[c];
];

BWT2D[a_,opts___]:=Module[{},Message[BWT2D::"nofilter"];Return[HWT2D[N[a],opts]]];

BWT2D[a_,{h__,hw___},opts___]:=Module[{rows,cols,maxits,its,b,j,d,z},
	If[MatrixQ[a,NumericQ]==False,Message[BWT2D::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[BWT2D::"baddimensions"]; Return[a];];
    If[{Length[{h}],Length[{hw}]}!={1,1} || Length[Flatten[{h,hw}]]==2,Message[BWT2D::"onefilter"];Return[WT2D[a,Flatten[{h,hw}],opts]]];
	{rows,cols}=Dimensions[a];
	
	(*Compute the maximum number of iterations.*)

	maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[BWT2D];
	If[!IntegerQ[its],Message[BWT2D::"badits",maxits];its=1];
	If[its<0,Message[BWT2D::"badits",maxits];its=1];
	If[its>maxits,Message[BWT2D::"maxits",its,maxits];its=maxits];

    b=a;
	For[j=0,j<its,j++,
		d=GetCorner[b,rows/(2^j),cols/(2^j)];
  		z=BWT2D1[d,{h,hw},opts];
    	b=PutCorner[b,z];
    ];
	Return[b];
];

IBWT2D[a_,opts___]:=Module[{},Message[IBWT2D::"nofilter"];Return[IHWT2D[N[a],opts]]];

IBWT2D[a_,{h__,hw___},opts___]:=Module[{rows,cols,maxits,its,b,j,d,z},
	If[MatrixQ[a,NumericQ]==False,Message[IBWT2D::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[IBWT2D::"baddimensions"]; Return[a];];
    If[{Length[{h}],Length[{hw}]}!={1,1} || Length[Flatten[{h,hw}]]==2,Message[IBWT2D::"onefilter"];Return[IWT2D[a,Flatten[{h,hw}],opts]]];
	{rows,cols}=Dimensions[a];
	
	(*Compute the maximum number of iterations.*)

	maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[IBWT2D];
	If[!IntegerQ[its],Message[IBWT2D::"badits",maxits];its=1];
	If[its<0,Message[IBWT2D::"badits",maxits];its=1];
	If[its>maxits,Message[IBWT2D::"maxits",its,maxits];its=maxits];

    b=a;
	For[j=1,j<=its,j++,
		d=GetCorner[b,rows/(2^(its-j)),cols/(2^(its-j))];
	    z=IBWT2D1[d,{h,hw},opts];
    	b=PutCorner[b,z];
    ];
	Return[b];
];


LWT2D1[a_,opts___]:=Module[{b,c},
	If[MatrixQ[a,NumericQ]==False,Message[LWT2D1::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[LWT2D1::"baddimensions"];Return[a];];
    
    b=Transpose[Map[LWT1D1[#,opts]&,Transpose[a]]];
    c=Map[LWT1D1[#,opts]&,b];
    Return[c];
];

ILWT2D1[a_,opts___]:=Module[{b,c},
	If[MatrixQ[a,NumericQ]==False,Message[ILWT2D1::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[ILWT2D1::"baddimensions"];Return[a];];
    
    b=Transpose[Map[ILWT1D1[#,opts]&,Transpose[a]]];
    c=Map[ILWT1D1[#,opts]&,b];
    Return[c];
];

LWT2D[a_,opts___]:=Module[{rows,cols,maxits,its,b,j,d,z},
	If[MatrixQ[a,NumericQ]==False,Message[LWT2D::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[LWT2D::"baddimensions"];Return[a];];
    {rows,cols}=Dimensions[a];
	
	(*Compute the maximum number of iterations.*)

	maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[LWT2D];
	If[!IntegerQ[its],Message[LWT2D::"badits",maxits];its=1];
	If[its<0,Message[LWT2D::"badits",maxits];its=1];
	If[its>maxits,Message[LWT2D::"maxits",its,maxits];its=maxits];

    b=a;
	For[j=0,j<its,j++,
		d=GetCorner[b,rows/(2^j),cols/(2^j)];
  		z=LWT2D1[d,opts];
    	b=PutCorner[b,z];
    ];
	Return[b];
];

ILWT2D[a_,opts___]:=Module[{rows,cols,maxits,its,b,j,d,z},
	If[MatrixQ[a,NumericQ]==False,Message[ILWT2D::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[ILWT2D::"baddimensions"]; Return[a];];
    {rows,cols}=Dimensions[a];
	
	(*Compute the maximum number of iterations.*)

	maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[ILWT2D];
	If[!IntegerQ[its],Message[ILWT2D::"badits",maxits];its=1];
	If[its<0,Message[ILWT2D::"badits",maxits];its=1];
	If[its>maxits,Message[ILWT2D::"maxits",its,maxits];its=maxits];

    b=a;
	For[j=1,j<=its,j++,
		d=GetCorner[b,rows/(2^(its-j)),cols/(2^(its-j))];
	    z=ILWT2D1[d,opts];
    	b=PutCorner[b,z];
    ];
	Return[b];
];

DCT1D[v_]:=Module[{d},
	If[VectorQ[v,NumericQ]==False,Message[DCT1D::"badinput"];Return[v];];
   	d=Join[{1},Table[Sqrt[2],{Length[v]-1}]];
	Return[d*FourierDCT[v,2]];
];

IDCT1D[v_]:=Module[{n,W,j,k},
	If[VectorQ[v,NumericQ]==False,Message[IDCT1D::"badinput"];Return[v];];
   	d=Join[{1},Table[1/Sqrt[2],{Length[v]-1}]];
	Return[FourierDCT[d*v,3]];
];

DCT2D[a_]:=Module[{r,c,d1,d2},
	If[MatrixQ[a,NumericQ]==False,Message[DCT2D::"badinput"];Return[a];];
	{r,c}=Dimensions[a];
	d1=DiagonalMatrix[Join[{1},Table[Sqrt[2],{r-1}]]];
	d2=DiagonalMatrix[Join[{1},Table[Sqrt[2],{c-1}]]];
	Return[d1.FourierDCT[a,2].d2];
];
	
IDCT2D[a_]:=Module[{r,c,d1,d2},
	If[MatrixQ[a,NumericQ]==False,Message[IDCT2D::"badinput"];Return[a];];
	{r,c}=Dimensions[a];
	d1=DiagonalMatrix[Join[{1},Table[1/Sqrt[2],{r-1}]]];
	d2=DiagonalMatrix[Join[{1},Table[1/Sqrt[2],{c-1}]]];
	Return[FourierDCT[d1.a.d2,3]];
];



WaveletMatrixToList[a_,opts___]:=Module[{rows,cols,maxits,its,hp,b,k},
	If[MatrixQ[a,NumericQ]==False,Message[WaveletMatrixToList::"badinput"];Return[a];];
	If[Map[EvenQ,Dimensions[a]]!={True,True},Message[WaveletMatrixToList::"baddimensions"];Return[a];];
    {rows,cols}=Dimensions[a];
    
    maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[WaveletMatrixToList];
	If[!IntegerQ[its],Message[WaveletMatrixToList::"badits",maxits];its=1];
	If[its<0,Message[WaveletMatrixToList::"badits",maxits];its=1];
	If[its>maxits,Message[WaveletMatrixToList::"maxits",its,maxits];its=maxits];

    hp=Table[Table[{},{3}],{its}];
    For[k=1,k<=its,k++,
    	b=Partition[a,{rows/2^k,cols/2^k}];
       	hp[[its-k+1,1]]=b[[1,2]];
       	hp[[its-k+1,2]]=b[[2,1]];
       	hp[[its-k+1,3]]=b[[2,2]];
       	b=b[[1,1]];
    ];
    Return[Join[{b},hp]];
];
	
WaveletListToMatrix[a_,opts___]:=Module[{rows,cols,maxits,its,k,m,tmp},
	If[MatrixQ[First[a],NumericQ]==False,Message[WaveletListToMatrix::"badinput"];Return[a];];
	If[Union[Map[MatrixQ[#,NumericQ]&,Last[a]]]!={True},Message[WaveletListToMatrix::"badinput"];Return[a];];
	{rows,cols}=Dimensions[Last[a][[1]]]*2;
	
    maxits=Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];

	its=NumIterations/.{opts}/.Options[WaveletListToMatrix];
	If[!IntegerQ[its],Message[WaveletListToMatrix::"badits",maxits];its=1];
	If[its<0,Message[WaveletListToMatrix::"badits",maxits];its=1];
	If[its>maxits,Message[WaveletListToMatrix::"maxits",its,maxits];its=maxits];
	
	(*tmp=Flatten[Union[Prepend[Map[Dimensions,Map[First,Drop[a,1]]]*Reverse[\
Table[2^k,{k,1,its}]],Dimensions[First[wtl]]*2^its]]];
	If[Length[tmp]!=2,Message[WaveletListToMatrix::"baddimensions"];Return[];];*)
      
    m=First[a];
    For[k=1,k<=its,k++,
    	m=Join[Join[Transpose[Join[Transpose[m],Transpose[a[[k+1,1]]]]]],Join[Transpose[Join[Transpose[a[[k+1,2]]],Transpose[a[[k+1,3]]]]]]]];
 	Return[m];
 ];  
 
 ChopVector[v_,opts___]:=Module[{n,nt,chop,ph},
      n=Length[v];
      
      ph=PrintInfo/.{opts}/.Options[ChopVector];
      
      nt=PowersOfTwo/.{opts}/.Options[ChopVector];
      
      If[nt===None,
      	If[ph===True,Print["The length of the vector is ",n,"."]];
      	Return[v]];
      If[2^nt>n,Message[ChopVector::"maxnt",nt,n];Return[v];];
      
      chop=Mod[n,2^nt];
      
      If[ph===True,Print["The length of the vector is ",n-chop,"."]];
      
      Return[Drop[v,-chop]];
      ];


      
GammaCorrection[a_,r_:1]:=Module[{},
		If[MatrixQ[a,IntegerQ]==False,Message[GammaCorrection::"badinput"];Return[a];];
		If[Min[a]<0 || Max[a]>255,Message[GammaCorrection::"badrange"];Return[a];];
		If[NumericQ[r]==False,Message[GammaCorrection::"badvalue"];Return[a];];
		If[r<=0,Message[GammaCorrection::"positivevalue"];Return[a];];
		Return[Round[((a/255)^(r))*255]];
	];
	
	
MakeHistogramEQ[a_]:=Module[{p,lvl,cnt,k,n,m,i},
	If[MatrixQ[a,IntegerQ]==False,Message[MakeHistogramEQ::"badinput"];Return[a];];
	If[Min[a]<0 || Max[a]>255,Message[MakeHistogramEQ::"badrange"];Return[a];];
	p=Split[Sort[Flatten[a]]];
    {lvl,cnt}=Transpose[Table[{p[[k,1]],Length[p[[k]]]},{k,1,Length[p]}]];
    i=1;
    n=Table[0,{256}];
    For[m=0,m<=255,m++,
 		If[m==lvl[[i]],n[[m+1]]=cnt[[i]];i=i+1];
        If[i>Length[p],i=0];
    ];
    Return[n];
];

HistogramEQ[a_]:=Module[{n,s},
	If[MatrixQ[a,IntegerQ]==False,Message[MakeHistogramEQ::"badinput"];Return[a];];
	If[Min[a]<0 || Max[a]>255,Message[MakeHistogramEQ::"badrange"];Return[a];];
	n=MakeHistogramEQ[a];
    s=N[Drop[FoldList[Plus,0,n],1]/Apply[Times,Dimensions[a]]];
    Return[Map[Ceiling[255*s[[#+1]]]&,a]];
];

WaveletVectorToList[v_,opts___]:=Module[{n,u,k,maxits},
	n=Length[v];
      
   	If[ArrayDepth[v]!=1,Message[WaveletVectorToList::"nonvector"];Return[v];];

    If[Mod[n,2]==1,Message[WaveletVectorToList::"oddlength"];Return[v];];
    maxits=FactorInteger[n][[1,2]];
 	its=NumIterations/.{opts}/.Options[WaveletVectorToList];
    If[its>maxits,Message[WaveletVectorToList::"maxits",maxits];
            Return[];];
    If[its===Max,its=maxits];
    u={};
	For[k=1,k<=its,k++,
    	u=Prepend[u,Take[v,{n/2^k+1,n/2^(k-1)}]];];
	Return[Prepend[u,Take[v,n/2^its]]];
];

WaveletListToVector[v_]:=Flatten[v];

MAD[v_]:=Median[Abs[Flatten[v]-Median[Flatten[v]]]];

ShrinkageFunction[t_, l_:0]:=Sign[t]*Which[Abs[t]>=l,Abs[t]-l,True,0];

DonohoSure[v_]:=Module[{n,a,b,c,k,s,r,idx},
	n=Length[v];
	a=Sort[v^2];
    b=Drop[FoldList[Plus,0,a],1];
    c=Reverse[Table[k,{k,0,n-1}]];
    s=b+c*a;
    r=n-2*Table[k,{k,1,n}]+s;
    idx=Flatten[Position[r,Min[r]]][[1]];
    Return[Sqrt[a[[idx]]]];
];

TestSparseness[v_,opts___]:=Module[{vflag,w,n,s,u,tol,prt},
	w=Flatten[v];
	vflag=If[VectorQ[w,NumericQ]==True,1,0];
	If[vflag==0,Message[TestSparseness::badinput];Return[v]];
	n=Length[w];
	s=Total[w^2]/n-1;
	u=3*Log[2,n]/(2*Sqrt[n]);
	tol=If[s>u,MAD[w]*Sqrt[2*Log[n]]/.6745,DonohoSure[w]];
	prt=PrintResult/.{opts}/.Options[TestSparseness];
	If[prt==True,
		If[s<=u,Print["The input is sparse - returing the universal threshold value"],Print["The input is not sparse - returning the SureShrink threshold value."]];
	];
	Return[tol];
];

NoiseEstimate[a_,h_]:=Module[{vflag,mflag,wt,wtlist,hp,sigma},
	vflag=If[VectorQ[a,NumericQ]==True,1,0];
	mflag=If[MatrixQ[a,NumericQ]==True,1,0];
	If[vflag==0 && mflag==0,Message[NoiseEstimate::badinput];Return[0]];
	wt=If[vflag==1,WT1D1[a,h],WT2D1[a,h]];
	wtlist=If[vflag==1,WaveletVectorToList[wt,NumIterations->1],WaveletMatrixToList[wt,NumIterations->1]];
	hp=Flatten[Last[wtlist]];
	sigma=MAD[hp]/.6745;
	Return[sigma];
];

UniversalThreshold[a_,h_,opts___]:=Module[{vflag,mflag,its,maxits,sigma,len},
	vflag=If[VectorQ[a,NumericQ]==True,1,0];
	mflag=If[MatrixQ[a,NumericQ]==True,1,0];
	If[vflag==0 && mflag==0,Message[UniversalThreshold::badinput];Return[0]];
	maxits=If[vflag==1,FactorInteger[Length[a]][[1,2]],Min[Map[Last,Map[First,Map[FactorInteger,Dimensions[A]]]]]];
	its=NumIterations/.{opts}/.Options[UniversalThreshold];
	If[its==Max ,its=maxits];
	If[its>maxits,Print["NumIterations is a larger value than allowed."];Return[0]];
	sigma=NoiseEstimate[a,h];
	len=Length[Flatten[a]]-If[vflag==1,Length[Flatten[a]]/2^its,Apply[Times,Dimensions[a]/2^its]];
	Return[sigma*Sqrt[2*Log[len]]];
];

WaveletShrinkage[a_,h_,lambda_,opts___]:=Module[{vflag,mflag,maxits,its,wt,wtlist,hp,dep,iwt},
	vflag=If[VectorQ[a,NumericQ]==True,1,0];
	mflag=If[MatrixQ[a,NumericQ]==True,1,0];
	If[vflag==0 && mflag==0,Message[WaveletShrinkage::badinput];Return[0]];
	maxits=If[vflag==1,FactorInteger[Length[a]][[1,2]],Min[Map[Last,Map[First,Map[FactorInteger,Dimensions[a]]]]]];
	its=NumIterations/.{opts}/.Options[WaveletShrinkage];
	If[its==Max ,its=maxits];
	If[its>maxits,Message[WaveletShrinkage::numiterations];Return[0]];
	If[Length[lambda]==0 && lambda < 0,Message[WaveletShrinkage::lambdapositive];Return[a]];
	If[Length[lambda]>0 && Union[Map[NonNegative,Flatten[lambda]]]!={True},Message[WaveletShrinkage::lambdavectorpositive];Return[a]];
	If[vflag==1 && Length[lambda]!=0 && Length[lambda]!=its,Message[WaveletShrinkage::lambdalength,its];Return[a]];
	If[mflag==1 && Length[lambda]!=0 && Dimensions[lambda]!={its,3},Message[WaveletShrinkage::lambdamatrix,its];Return[a]];

	wt=If[vflag==1,WT1D[a,h,NumIterations->its],WT2D[a,h,NumIterations->its]];
	wtlist=If[vflag==1,WaveletVectorToList[wt,NumIterations->its],WaveletMatrixToList[wt,NumIterations->its]];
	hp=Drop[wtlist,1];
	dep=If[vflag==1,2,4];
	hp=If[Length[lambda]==0,Map[ShrinkageFunction[#,lambda]&,hp,{dep}],
	If[vflag==1,Table[Map[ShrinkageFunction[#,lambda[[k]]]&,hp[[k]]],{k,1,its}],Table[Table[Map[ShrinkageFunction[#,lambda[[j,k]]]&,hp[[j,k]],{2}],{k,1,3}],{j,1,its}]]];
	wtlist=Prepend[hp,First[wtlist]];
	wt=If[vflag==1,WaveletListToVector[wtlist],WaveletListToMatrix[wtlist,NumIterations->its]];
	iwt=If[vflag==1,IWT1D[wt,h,NumIterations->its],IWT2D[wt,h,NumIterations->its]];
	Return[iwt];
];

SureShrink[a_,h_,opts___]:=Module[{vflag,mflag,maxits,its,wt,wtlist,sigma,hp,u,s,p,j,k,lambda,iwt},
	vflag=If[VectorQ[a,NumericQ]==True,1,0];
	mflag=If[MatrixQ[a,NumericQ]==True,1,0];
	If[vflag==0 && mflag==0,Message[SureShrink::badinput];Return[a]];
	maxits=If[vflag==1,FactorInteger[Length[a]][[1,2]],Min[Map[Last,Map[First,Map[FactorInteger,Dimensions[a]]]]]];
	its=NumIterations/.{opts}/.Options[SureShrink];
	If[its==Max ,its=maxits];
	If[its>maxits,Message[SureShrink::maxits];Return[a]];
	wt=If[vflag==1,WT1D[a,h,NumIterations->its],WT2D[a,h,NumIterations->its]];
	wtlist=If[vflag==1,WaveletVectorToList[wt,NumIterations->its],WaveletMatrixToList[wt,NumIterations->its]];
	sigma=MAD[Flatten[Last[wtlist]]]/.6745;
	hp=Drop[wtlist,1]/sigma;
	If[mflag==1,p=Map[Flatten,hp,{2}]];
	u=N[If[vflag==1,Map[(3*Log[2,Length[#]]/(2*Sqrt[Length[#]]))&,hp],Map[(3*Log[2,Length[#]]/(2*Sqrt[Length[#]]))&,p,{2}]]];
	s=N[If[vflag==1,Map[(Total[#^2]/Length[#]-1)&,hp],Map[(Total[#^2]/Length[#]-1)&,p,{2}]]];
	lambda=If[vflag==1,Table[If[s[[k]]<=u[[k]],sigma*Sqrt[2*Log[Length[hp[[k]]]]],DonohoSure[hp[[k]]]],{k,1,its}],Table[Table[If[s[[j,k]]<=u[[j,k]],sigma*Sqrt[2*Log[Length[p[[j,k]]]]],DonohoSure[p[[j,k]]]],{k,1,3}],{j,1,its}]];
	hp=If[vflag==1,Table[Map[ShrinkageFunction[#,lambda[[k]]]&,hp[[k]]],{k,1,its}],Table[Table[Map[ShrinkageFunction[#,lambda[[j,k]]]&,hp[[j,k]],{2}],{k,1,3}],{j,1,its}]];
	wtlist=Prepend[hp*sigma,First[wtlist]];
	wt=If[vflag==1,WaveletListToVector[wtlist],WaveletListToMatrix[wtlist,NumIterations->its]];
	iwt=If[vflag==1,IWT1D[wt,h,NumIterations->its],IWT2D[wt,h,NumIterations->its]];
	Return[iwt];
];


ImageList[opts___]:=Module[{fname,f,basedir=$UserBaseDirectory,type,glist,clist,lst},
	type=ToString[ImageType/.{opts}/.Options[ImageList]];
	If[type!="GrayScale" && type!="Color", type="All"];
	
	f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
	If[Length[FileNames[f]]==0,basedir=$BaseDirectory];
	(*If[OpenRead[f]===$Failed,basedir=$UserBaseDirectory];*)
	Print["The base directory for the images is ",basedir," and the images courtesy of Dr. Radka Turcajova."];
	Print["The naming convention for the thumbnail image is to add a _small to the name.  For example, the thumbnail image for benches.png is benches_small.png."];	
	Print["The number in parentheses represents the dimensions and max iterations for the thumbnail image.\n"];

	If[type=="GrayScale" || type=="All",
		fname=basedir<>"/Applications/DiscreteWavelets/Images/GrayScale/Images.txt";
		f=OpenRead[fname];
		glist=ReadList[f];
		glist=Map[Insert[#,"Gray",2]&,glist];
		Close[f];
	];
	
	If[type=="Color" || type=="All",
		fname=basedir<>"/Applications/DiscreteWavelets/Images/Color/Images.txt";
		f=OpenRead[fname];
		clist=ReadList[f];
		clist=Map[Insert[#,"Color",2]&,clist];
		Close[f];
	];
	
	Switch[type,"GrayScale", \
lst=glist,"Color",lst=clist,"All",glist=Append[glist,{"","","","",""}];lst=\
Join[glist,clist]];
	Return[TableForm[lst,TableHeadings->{None,{"File Name", "Type", "Rows", \
"Columns", "Max Iterations"}},TableAlignments->Center]];
];


ImageNames[opts___]:=Module[{type,basedir=$UserBaseDirectory,f,fname,glist,clist,dir,names,thumblist,thumbnails,tglist,gtable},
	type=ToString[ImageType/.{opts}/.Options[ImageNames]];
	If[type!="GrayScale" && type!="Color", type="All"];
	thumblist=ToString[ListThumbnails/.{opts}/.Options[ImageNames]];
		
	f=basedir<>"\\Applications\\DiscreteWavelets\\DiscreteWavelets.m";
	If[Length[FileNames[f]]==0,basedir=$BaseDirectory];
	
	If[type=="GrayScale" || type=="All",
		dir=basedir<>"/Applications/DiscreteWavelets/Images/GrayScale/";
		fname=dir<>"Images.txt";
		f=OpenRead[fname];
		glist=ReadList[f];
		names=Map[First,glist];
		glist=Map[StringSplit,Map[ToString,names]];
		glist=Map[StringJoin[dir,#]&,glist];
		Close[f];
			
		If[thumblist=="True",glist=Map[StringInsert[#,"_small",-5]&,glist]];
	];
	
	If[type=="Color" || type=="All",
		dir=basedir<>"/Applications/DiscreteWavelets/Images/Color/";
		fname=dir<>"Images.txt";
		f=OpenRead[fname];
		clist=ReadList[f];
		names=Map[First,clist];
		clist=Map[StringSplit,Map[ToString,Map[First,clist]]];
		clist=Map[StringJoin[dir,#]&,clist];
		Close[f];
				
		If[thumblist=="True",clist=Map[StringInsert[#,"_small",-5]&,clist]];
	];
	Switch[type,"GrayScale",Return[glist],"Color",Return[clist],"All",Return[{glist,clist}]];
];


ShowThumbnails[opts___]:=Module[{type,lbls,thumbs,thumbnails,gtable,ctable,imgs,k,tmp},
	type=ToString[ImageType/.{opts}/.Options[ShowThumbnails]];
	If[type!="GrayScale" && type!="Color", type="All"];
	If[type=="GrayScale" || type=="All",
		lbls=Map[Last,Map[StringSplit[#,"/"]&,ImageNames[ImageType->"GrayScale"]]];
		lbls=Table["gray"<>ToString[k]<>" = "<>lbls[[k]],{k,1,Length[lbls]}];
		thumbs=ImageNames[ImageType->"GrayScale",ListThumbnails->True];
		thumbnails=Map[ImageRead,thumbs];
		gtable=Table[CreateImageObject[N[thumbnails[[k]]],PlotLabel->lbls[[k]]],{k,1,Length[thumbnails]}];
	];
	If[type=="Color" || type=="All",
		lbls=Map[Last,Map[StringSplit[#,"/"]&,ImageNames[ImageType->"Color"]]];
		lbls=Table["color"<>ToString[k]<>" = "<>lbls[[k]],{k,1,Length[lbls]}];
		thumbs=ImageNames[ImageType->"Color",ListThumbnails->True];
		thumbnails=Map[ImageRead,thumbs];
		ctable=Table[CreateImageObject[N[thumbnails[[k]]],PlotLabel->lbls[[k]]],{k,1,Length[thumbnails]}];
	];
	imgs=Switch[type,"GrayScale",gtable,"Color",ctable,"All",Flatten[{gtable,ctable}]];
	type=SlideShow/.{opts}/.Options[ShowThumbnails];
	tmp=If[type==True,SlideView[imgs],GraphicsGrid[Partition[imgs,3],ImageSize->768,Frame->All]];
	(*Return[GraphicsGrid[Partition[imgs,3],ImageSize->768,Frame->All]];*)
	Return[tmp];
];


	
AudioList[]:=Module[{fname,f,basedir=$UserBaseDirectory,alist},
	f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
	If[Length[FileNames[f]]==0,basedir=$BaseDirectory];
      (*f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
      If[OpenRead[f]===$Failed,basedir=$UserBaseDirectory];*)
      fname=basedir<>"/Applications/DiscreteWavelets/Sounds/sounds.txt";
      Print["The base directory for the audio files is ",basedir,"."];
      f=OpenRead[fname];
      alist=ReadList[f];
      Return[
        TableForm[alist,
          TableHeadings->{None,{"File Name","Type","Duration (s)",
                "Channels","Sample Rate","Samples"}},
          TableAlignments->Center]];
];	

AudioNames[]:=Module[{type,basedir=$UserBaseDirectory,f,fname,alist,names},
	f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
	If[Length[FileNames[f]]==0,basedir=$BaseDirectory];
      (*f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
      If[OpenRead[f]===$Failed,basedir=$UserBaseDirectory];*)
      dir=basedir<>"/Applications/DiscreteWavelets/Sounds/";
      fname=dir<>"sounds.txt";
      f=OpenRead[fname];
      alist=ReadList[f];
      names=Map[First,alist];
      alist=Map[StringSplit,Map[ToString,names]];
      alist=Map[StringJoin[dir,#]&,alist];
      Close[f];
      Return[alist];
      ];

DataList[]:=Module[{fname,f,basedir=$UserBaseDirectory,alist},
	f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
	If[Length[FileNames[f]]==0,basedir=$BaseDirectory];
      (*f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
      If[OpenRead[f]===$Failed,basedir=$UserBaseDirectory];*)
      fname=basedir<>"/Applications/DiscreteWavelets/Data/data.txt";
      Print["The base directory for the audio files is ",basedir,"."];
      f=OpenRead[fname];
      alist=ReadList[f];
      Return[
        TableForm[alist,
          TableHeadings->{None,{"File Name","Length","Description"}},
          TableAlignments->Center]];
];	

DataNames[]:=Module[{type,basedir=$UserBaseDirectory,f,fname,alist,names},
	f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
	If[Length[FileNames[f]]==0,basedir=$BaseDirectory];
      (*f=basedir<>"/Applications/DiscreteWavelets/DiscreteWavelets.m";
      If[OpenRead[f]===$Failed,basedir=$UserBaseDirectory];*)
      dir=basedir<>"/Applications/DiscreteWavelets/Data/";
      fname=dir<>"data.txt";
      f=OpenRead[fname];
      alist=ReadList[f];
      names=Map[First,alist];
      alist=Map[StringSplit,Map[ToString,names]];
      alist=Map[StringJoin[dir,#]&,alist];
      Close[f];
      Return[alist];
      ];

MaximumIterations[a___]:=Module[{cflag,rows,cols},
    If[MatrixQ[a,NumericQ]==True,cflag=0;];
    If[Map[MatrixQ[#,NumericQ]&,a]=={True,True,True},cflag=1;];
    If[cflag==1 && \
Length[Flatten[Union[Map[Dimensions,a]]]]!=2,Message[MaximumIterations::"differentdimensions"];Return[0];];
    If[cflag==-1,Message[MaximumIterations::"badinput"];Return[0];];
    {rows,cols}=If[cflag==0,Dimensions[a],Dimensions[First[a]]];
    If[Union[Map[EvenQ,{rows,cols}]]!={True},Return[0];];

    Return[Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]]];
];


ImageRead[file_,opts___]:=Module[{f,r,c,maxnt,chopr,chopc,img,red,green,blue,\
cflag,str,ph},
	(*f=Import[file];
	If[!(Head[f]===Graphics),Message[ImageRead::"badinput",file];Return[{{0,0},{\
0,0}}]];
    img=f[[1,1]];*)
	img=Import[file,"Data"];
	cflag=Length[Dimensions[img]];
	If[cflag<3 && MatrixQ[img,NumericQ]===False,Message[ImageRead::"badinput",file];Return[{{0,0},{\
0,0}}]];
	If[cflag>=3 && Cases[Map[MatrixQ[#,NumericQ]&,img],False]!={},Message[ImageRead::"badinput",file];Return[{{0,0},{\
0,0}}]];
    cflag=Length[Dimensions[img]];
	If[cflag==3,
		{red,green,blue}={img[[All,All,1]],img[[All,All,2]],img[[All,All,3]]};
		str="color",
		red=img;
		str="grayscale"];
    (*If[cflag==3,
        red=Reverse[img[[All,All,1]]];
        green=Reverse[img[[All,All,2]]];
        blue=Reverse[img[[All,All,3]]];
        str="color",
        red=Reverse[img];
        str="grayscale"];*)
    {r,c}=Dimensions[red];
                  
    nt=PowersOfTwo/.{opts}/.Options[ImageRead];
    If[!IntegerQ[nt],nt=0];
    ph=PrintInfo/.{opts}/.Options[ImageRead];
       
    If[(nt===None)&&(ph===True),Print["The dimensions of the image are ",r," \
x ",c," and the file was processed as a ",str," image.  A total of \
",MaximumIterations[red]," iterations of a wavelet transform can be performed \
on the image."]]; 
    If[(nt===None)&&(cflag==3),Return[{red,green,blue}]];
    If[(nt===None)&&(cflag!=3),Return[red]];
    If[2^nt>Min[r,c],Message[ImageRead::"maxnt",nt,r,c];Return[];];
      
    chopr=Mod[r,2^nt];
    chopc=Mod[c,2^nt];
    If[!(nt===None) && Norm[{chopr,chopc}]!=0,Print["Chopping ",chopc," \
columns off the right and ",chopr," rows off the bottom of the image.  The \
dimensions are now divisible by 2^",nt,"."]];
    red=Take[red,{1,r-chopr},{1,c-chopc}];
      
    If[ph===True,Print["The dimensions of the image are ",r-chopr," x \
",c-chopc," and the file was processed as a ",str," image.  A total of \
",MaximumIterations[red]," iterations of a wavelet transform can be performed \
on the image."]];
      
    If[cflag==3,
        green=Take[green,{1,r-chopr},{1,c-chopc}];
        blue=Take[blue,{1,r-chopr},{1,c-chopc}];
        Return[{red,green,blue}],Return[red]];
];

ReadImage[file_,opts___]:=ImageRead[file,opts];


ImagePlot[a_,opts___]:=Module[{},
	If[MatrixQ[a,NumericQ]==True,cflag=0;];
    If[Map[MatrixQ[#,NumericQ]&,a]=={True,True,True},cflag=1;];
    If[cflag==1 && Length[Flatten[Union[Map[Dimensions,a]]]]!=2,Message[CreateImageObject::"differentdimensions"];Return[];];
    If[cflag==-1,Message[CreateImageObject::"badinput"];Return[];];
    
	Return[Show[CreateImageObject[a,opts]]];
];


CreateImageObject[a_,opts___]:=Module[{m,n,bm,bn,is,clr,c,obj,scaling,red,green,blue,rst,cflag},

	If[MatrixQ[a,NumericQ]==True,cflag=0;];
    If[Map[MatrixQ[#,NumericQ]&,a]=={True,True,True},cflag=1;];
    If[cflag==1 && Length[Flatten[Union[Map[Dimensions,a]]]]!=2,Message[CreateImageObject::"differentdimensions"];Return[];];
    If[cflag==-1,Message[CreateImageObject::"badinput"];Return[];];

	If[cflag==0,red=a,{red,green,blue}=a];
    
    
    (* Use Radka Turcajova's fix on the image border if so requested. *)
    
    is=ImageSize/.{opts}/.Options[ImagePlot];
    {m,n}=Dimensions[red];
    If[is===FixEdges,
      bm=Round[(m+1)/40]; bm=2*bm + 1 + Sign[m- 40*bm +1];
      bn=Round[(n+1)/40]; bn=2*bn + 1 + Sign[n- 40*bn+1];is={n+bn,m+bm}];
     
    
    (* Now determine the color to use for ColorFunction.  
          The user can enter their own color function or use ChannelColor.  
          ChannelColor only works with RGBColor. *)
    
    If[cflag==0,
      clr=ChannelColor/.{opts}/.Options[ImagePlot];
      If[(clr[[0]]===CMYKColor)||(clr[[0]]===GrayLevel)||(clr[[0]]===Hue),
      Message[CreateImageObject::"onlyRGB"];Return[];];
      c={1.,1.,1.};
      If[clr=!=None,c[[1]]=clr[[1]];c[[2]]=clr[[2]];c[[3]]=clr[[3]]]
      ];

    (* Apply LinMap if requested. *)
    
    scaling=LinearScaling/.{opts}/.Options[ImagePlot];
    If[scaling===LeftWT && \
OddQ[m]==True,Message[CreateImageObject::"LeftWT"];scaling=False];
    If[scaling===LeftWT && \
OddQ[n]==True,Message[CreateImageObject::"RightWT"];scaling=False];
    
    If[scaling===True,
        	red=LinMap[Abs[red],255];
        	If[cflag==1,green=LinMap[Abs[green],255];
          blue=LinMap[Abs[blue],255];];
        ];
        
    If[scaling===LeftWT,
    	red = \
Join[LinMap[Abs[Take[red,m/2]],255],LinMap[Abs[Drop[red,m/2]],255]];
    	If[cflag==1,
    		green = \
Join[LinMap[Abs[Take[green,m/2]],255],LinMap[Abs[Drop[green,m/2]],255]];
    		blue = \
Join[LinMap[Abs[Take[blue,m/2]],255],LinMap[Abs[Drop[blue,m/2]],255]];
    	];
    ];

    If[scaling===RightWT,
    	red = \
Transpose[Join[LinMap[Abs[Take[Transpose[red],n/2]],255],LinMap[Abs[Drop[\
Transpose[red],n/2]],255]]];
    	If[cflag==1,
    		green = \
Transpose[Join[LinMap[Abs[Take[Transpose[green],n/2]],255],LinMap[Abs[Drop[\
Transpose[green],n/2]],255]]];
    		blue = \
Transpose[Join[LinMap[Abs[Take[Transpose[blue],n/2]],255],LinMap[Abs[Drop[\
Transpose[blue],n/2]],255]]];
    	];
    ];
    
    If[cflag==0,
        rst=Raster[Reverse[red],{{0,0},Reverse[Dimensions[red]]},{0,255},
            ColorFunction->(RGBColor[c[[1]]*#,c[[2]]*#,c[[3]]*#1]&),
            "ColorFunctionScaling"->ColorFunctionScaling/.{opts}/.Options[ImagePlot]],
        
         rst= Raster[Partition[Transpose[{Flatten[Reverse[red]],Flatten[Reverse[green]],Flatten[Reverse[blue]]}],n]/256];];
    
    obj=Graphics[rst,
    	"AlignmentPoint"->AlignmentPoint/.{opts}/.Options[ImagePlot],
          "AspectRatio"->AspectRatio/.{opts}/.Options[ImagePlot],
          "Axes"->Axes/.{opts}/.Options[ImagePlot],
          "AxesLabel"->AxesLabel/.{opts}/.Options[ImagePlot],
          "AxesOrigin"->AxesOrigin/.{opts}/.Options[ImagePlot],
          "AxesStyle"->AxesStyle/.{opts}/.Options[ImagePlot],
          "Background"->Background/.{opts}/.Options[ImagePlot],
          "BaselinePosition"->BaselinePosition/.{opts}/.Options[ImagePlot],
          "BaseStyle"->BaseStyle/.{opts}/.Options[ImagePlot],
          "ColorOutput"->ColorOutput/.{opts}/.Options[ImagePlot],
          "ContentSelectable"->ContentSelectable/.{opts}/.Options[ImagePlot],
          "DisplayFunction":>DisplayFunction/.{opts}/.Options[ImagePlot],
          "Epilog"->Epilog/.{opts}/.Options[ImagePlot],
          "FormatType":>FormatType/.{opts}/.Options[ImagePlot],
          "Frame"->Frame/.{opts}/.Options[ImagePlot],
          "FrameLabel"->FrameLabel/.{opts}/.Options[ImagePlot],
          "FrameStyle"->FrameStyle/.{opts}/.Options[ImagePlot],
          "FrameTicks"->FrameTicks/.{opts}/.Options[ImagePlot],
          "FrameTicksStyle"->FrameTicksStyle/.{opts}/.Options[ImagePlot],
          "GridLines"->GridLines/.{opts}/.Options[ImagePlot],
          "GridLinesStyle"->GridLinesStyle/.{opts}/.Options[ImagePlot],
          "ImageMargins"->ImageMargins/.{opts}/.Options[ImagePlot],
          "ImagePadding"->ImagePadding/.{opts}/.Options[ImagePlot],
          ImageSize->is,
          "LabelStyle"->LabelStyle/.{opts}/.Options[ImagePlot],
          "Method"->Method/.{opts}/.Options[ImagePlot],
          "PlotLabel"->PlotLabel/.{opts}/.Options[ImagePlot],
          "PlotRange"->PlotRange/.{opts}/.Options[ImagePlot],
          "PlotRangeClipping"->PlotRangeClipping/.{opts}/.Options[ImagePlot],
          "PlotRangePadding"->PlotRangePadding/.{opts}/.Options[ImagePlot],
          "PlotRegion"->PlotRegion/.{opts}/.Options[ImagePlot],
          "PreserveImageOptions"->PreserveImageOptions/.{opts}/.Options[ImagePlot],
          "Prolog"->Prolog/.{opts}/.Options[ImagePlot],
          "RotateLabel"->RotateLabel/.{opts}/.Options[ImagePlot],
          "Ticks"->Ticks/.{opts}/.Options[ImagePlot],
          "TicksStyle"->TicksStyle/.{opts}/.Options[ImagePlot]];
    Return[obj];
];


	
LinMap[A_,mx_]:=Module[{dx},
	If[Min[A]==Max[A],Return[A],dx=mx/(Max[A]-Min[A]);
    If[mx<=1,Return[dx*(A-Min[A])], Return[Round[dx*(A-Min[A])]]]];];

WaveletDensityPlot[a_,opts___]:=Module[{startr,startc,endr,endc,its,dividelines,thk,clr,bx,hlns,clns,
        rows,cols,k,obj,region,cflag=-1,maxits},
        
    If[MatrixQ[a,NumericQ]==True,cflag=0;];
    If[Map[MatrixQ[#,NumericQ]&,a]=={True,True,True},cflag=1;];
    If[cflag==1 && Length[Flatten[Union[Map[Dimensions,a]]]]!=2,Message[WaveletDensityPlot::"differentdimensions"];Return[];];
    If[cflag==-1,Message[WaveletDensityPlot::"badinput"];Return[];];
    {rows,cols}=If[cflag==0,Dimensions[a],Dimensions[First[a]]];
    
    (* At this point, we know we have either one numeric matrix or three numeric matrices all of the same dimensions.  Now we need to
    	check the evenness of the dimensions.*)
    
    If[Union[Map[EvenQ,{rows,cols}]]!={True},Message[WaveletDensityPlot::"baddimensions"];Return[];];
    
    (* Compute the maximum number of iterations. *)
    maxits = Min[Map[Last,Map[First,Map[FactorInteger,{rows,cols}]]]];
          
	its=NumIterations/.{opts}/.Options[WaveletDensityPlot];
    If[!IntegerQ[its],Message[WaveletDensityPlot::"badits",maxits];its=1];
    If[its<0,Message[WaveletDensityPlot::"badits",maxits];its=1];
    If[its>maxits,Message[WaveletDensityPlot::"maxits",its,maxits];its=maxits];
     
    If[cflag==0,
        region=ExtractRegion[wtLDP[N[a],its],opts],
        region=
          ExtractRegion[{wtLDP[N[a[[1]]],its],wtLDP[N[a[[2]]],its],
              wtLDP[N[a[[3]]],its]},opts]
        ];
      
      {startr,endr}=region[[2]];
      {startc,endc}=region[[3]];
      
      obj=CreateImageObject[region[[1]],opts];
      
      divlines=DivideLines/.{opts}/.Options[WaveletDensityPlot];
      reg=Region/.{opts}/.Options[WaveletDensityPlot];
      
      If[(divlines===True)&&(reg===All),
        thk=DivideLinesThickness/.{opts}/.Options[WaveletDensityPlot];
        clr=DivideLinesColor/.{opts}/.Options[WaveletDensityPlot];
        iter=Log[2,rows/endr]+1;
        rows=endr-startr+1;
        cols=endc-startc+1;
        bx=
          Graphics[{Thickness[thk],clr,
              Line[{{0,0},{cols,0},{cols,rows},{0,rows},{0,0}}]}];
        hlns=
          Table[Graphics[{Thickness[thk],clr,
                Line[{{0,(2^k-1)*rows/2^k},{cols/2^(k-1),(2^k-1)*
                        rows/2^k}}]}],{k,1,its-iter+1}];
        	clns=
          Table[Graphics[{Thickness[thk],clr,
                Line[{{cols/2^(k+1),(2^k-1)*rows/2^k},{cols/2^(k+1),
                      rows}}]}],{k,0,its-iter}];
        Return[Show[{obj,bx,hlns,clns}]];
        ];
	Return[Show[obj]];
];

wtLDP[a_,its_]:=Module[{m,mat,k,i,r,c},{r,c}=Dimensions[a];
      m=Table[{0},{k,1,its}];
      mat=Table[{0},{k,1,its}];m[[1]]=Partition[a,{r/2,c/2}];
      For[i=2,i<=its,m[[i]]=Partition[m[[i-1,1,1]],{r/2^i,c/2^i}];
        i++];
      m[[its,1,1]]=LinMap[m[[its,1,1]],255];
      For[i=1,i<=its,m[[i,1,2]]=LinMap[Abs[m[[i,1,2]]],255];
        m[[i,2,1]]=LinMap[Abs[m[[i,2,1]]],255];
        m[[i,2,2]]=LinMap[Abs[m[[i,2,2]]],255];i++];
      mat[[1]]=
        Join[Table[
            Flatten[Append[m[[its,1,1,i]],m[[its,1,2,i]]]],{i,1,r/2^its}],
          Table[Flatten[Append[m[[its,2,1,i]],m[[its,2,2,i]]]],{i,1,
              r/2^its}]];
      For[k=2,k<=its,
        mat[[k]]=
          Join[Table[
              Flatten[Append[mat[[k-1,i]],m[[its-k+1,1,2,i]]]],{i,1,
                r/2^(its-k+1)}],
            Table[Flatten[Append[m[[its-k+1,2,1,i]],m[[its-k+1,2,2,i]]]],{i,1,\
r/2^(its-k+1)}]];
        k++];
      Return[mat[[its]]]
];

ExtractRegion[a_,opts___]:=Module[{startr,startc,endr,endc,iter,reg,its,k,cflag=-1,maxits},
    startr=startc=1;
    If[MatrixQ[a,NumericQ]==True,cflag=0;];
    If[Map[MatrixQ[#,NumericQ]&,a]=={True,True,True},cflag=1;];
    If[cflag==1 && Length[Flatten[Union[Map[Dimensions,a]]]]!=2,Message[ExtractRegion::"differentdimensions"];Return[];];
    If[cflag==-1,Message[ExtractRegion::"badinput"];Return[];];

	{endr,endc}=If[cflag==0,Dimensions[a],Dimensions[First[a]]];
    
    (* At this point, we know we have either one numeric matrix or three numeric matrices all of the same dimensions.  Now we need to
    	check the evenness of the dimensions.*)
    
    If[Union[Map[EvenQ,{endr,endc}]]!={True},Message[ExtractRegion::"baddimensions"];Return[];];
    
      
    (* Compute the maximum number of iterations. *)
    maxits = Min[Map[Last,Map[First,Map[FactorInteger,{endr,endc}]]]];
          
	its=NumIterations/.{opts}/.Options[ExtractRegion];
    If[!IntegerQ[its],Message[ExtractRegion::"badits",maxits];its=1];
    If[its<0,Message[ExtractRegion::"badits",maxits];its=1];
    If[its>maxits,Message[ExtractRegion::"maxits",its,maxits];its=maxits];


    iter=Iteration/.{opts}/.Options[ExtractRegion];
    If[iter>maxits, iter=maxits];
    reg=SymbolName[Region/.{opts}/.Options[ExtractRegion]];
      
    (* Now take care of some special cases.  If Region is Blur then we want the last iteration, so iter=its.  
       The same is the case if iter is All and Region is not All. *)
  
      If[(reg=="Blur"||(iter===All&&reg!="All")),iter=its];
      If[iter===All,iter=1];
      
      (* Once we have iter, we can set the end of rows and the end of columns. *)
      
      endr=endr/2^(iter-1);
      endc=endc/2^(iter-1);
      
      (* start and stop points now depend on specific regions that were picked and easy to figure out. *)
      
      If[reg=="Blur"||reg=="Vertical",endr=endr/2];
      If[reg=="Horizontal"||reg=="Diagonal", startr=IntegerPart[(startr+endr)/2+1]];
      If[reg=="Blur"||reg=="Horizontal", endc=endc/2];
      If[reg=="Vertical"||reg=="Diagonal", startc=IntegerPart[(startc+endc)/2+1]];
      
      (* Return a list consisting of the extracted matric, the starting/stopping row, and the starting/stopping column. *)
      
      If[cflag==1,
        Return[{Table[
              Take[a[[k]],{startr,endr},{startc,endc}],{k,1,3}],{startr,
              endr},{startc,endc}}],
        Return[{Take[a,{startr,endr},{startc,endc}],{startr,endr},{startc,
              endc}}]];
 ];  
 
 ExtractVectorPart[v_,opts___]:=Module[{n,region,nits,iter,start,stop},
      n=Length[v];
      start=1;
      stop=n;
      region=Region/.{opts}/.Options[ExtractVectorPart];
      
      nits=NumIterations/.{opts}/.Options[ExtractVectorPart];
      iter=Iteration/.{opts}/.Options[ExtractVectorPart];
      If[nits===Max,nits=FactorInteger[n][[1,2]]];
      If[iter===All,iter=1]; 
      If[ToString[region]=="LowPass",stop=n/2^nits];
      
      If[ToString[region]=="HighPass",start=n/2^iter+1;stop=start+n/2^iter-1];

      If[ToString[region]=="All",stop=n/2^(iter-1)]; 
      Return[Take[v,{start,stop}]];
      ];
      
WaveletVectorPlot[v_,opts___]:=Module[{nits,iter,region,parts,usecolors,color,lc=0,w,m,data,pts,k,init,
        lens,pthk,clrs,clr,thk,divlines,lns,bx,mx,mn,xv,n,plt},
	n=Length[v];
    nits=NumIterations/.{opts}/.Options[WaveletVectorPlot];
      
    If[nits===Max,nits=FactorInteger[n][[1,2]]];
      
    iter=Iteration/.{opts}/.Options[WaveletVectorPlot];
    If[iter===All,iter=1;parts=nits+1];
      
    region=Region/.{opts}/.Options[WaveletVectorPlot];
    (* Figure out how many parts there are to color. *)
      
    If[region===All,parts=nits-iter+2,parts=1];
    
    pthk=PointSize/.{opts}/.Options[WaveletVectorPlot];
    usecolors=UseColors/.{opts}/.Options[WaveletVectorPlot];
      
    If[usecolors===True,color=ColorList/.{opts}/.Options[WaveletVectorPlot];
    lc=Length[color];
    (* Make sure we have enough colors!! *)
        
    While[lc<parts,color=Join[color,color];lc=Length[color]]];
    (* Now extract the portion of the vector we want to plot and then make a table of ordered pairs out of it.  
    	We will use these ordered pairs to build a Graphics object. *)
   
    w=ExtractVectorPart[v,opts];
    m=Length[w];
    data=Point/@Transpose[{Table[k,{k,1,m}],w}];
      
    (* No colors - just return the graphics object. *)
      
    If[usecolors===False,pts=Transpose[{Table[PointSize[pthk],{m}],data}]];
      
    (* If region is not All, then the user only wanted one part - 
          assign the graphics object the first color in the list and return. *)
      
    If[(usecolors===True)&&(region=!=All),
    	pts=Transpose[{Table[PointSize[pthk],{m}],Table[color[[1]],{k,1,m}],data}]];
      
    (* If the region is everything, then we need to figure out the lengths of each part.  
            That's what we are doing here. *)
      
    If[(usecolors===True)&&(region===All),
        init={m/2^(parts-1),m/2^(parts-1)};
        lens=Join[init,Table[m/2^(parts-k+1),{k,3,parts}]];
        
    (* Now join a color to each bunch of points. *)
        
    clrs= Join[Table[color[[1]],{k,1,lens[[1]]}],
            Table[color[[2]],{k,1,lens[[2]]}]];
    For[k=3,k<=parts,k++,
    	clrs=Join[clrs,Table[color[[k]],{j,1,lens[[k]]}]]];
        
        pts=Transpose[{Table[PointSize[pthk],{m}],clrs,data}];
    ];
      
      
    (* Add the divide lines if so desired. *)
      
    divlines=DivideLines/.{opts}/.Options[WaveletVectorPlot];
    lns={};
    bx={};
    If[(divlines===True)&&(region===All),
        init={m/2^(parts-1),m/2^(parts-1)};
        lens=Join[init,Table[m/2^(parts-k+1),{k,3,parts}]];
    	thk=DivideLinesThickness/.{opts}/.Options[WaveletVectorPlot];
    	clr=DivideLinesColor/.{opts}/.Options[WaveletVectorPlot];
    	mn=Min[v];
    	mx=Max[v];
    	xv=Drop[FoldList[Plus,1,lens],1];
    	lns=Table[Graphics[{clr,Thickness[thk],
                Line[{{xv[[k]],mn},{xv[[k]],mx}}]}],{k,1,Length[lens]-1}];
    	bx=Graphics[{Thickness[thk],clr,
              Line[{{0,mn},{m,mn},{m,mx},{0,mx},{0,mn}}]}];
    ];
      
    plt = Show[{bx,lns,Graphics[pts]},
    		"AlignmentPoint"->AlignmentPoint/.{opts}/.Options[WaveletVectorPlot],
            "AspectRatio"->AspectRatio/.{opts}/.Options[
                WaveletVectorPlot],
            "Axes"->Axes/.{opts}/.Options[WaveletVectorPlot],
            "AxesLabel"->AxesLabel/.{opts}/.Options[WaveletVectorPlot],
            "AxesOrigin"->AxesOrigin/.{opts}/.Options[WaveletVectorPlot],
            "AxesStyle"->AxesStyle/.{opts}/.Options[WaveletVectorPlot],
            "Background"->Background/.{opts}/.Options[WaveletVectorPlot],
            "BaselinePosition"->BaselinePosition/.{opts}/.Options[WaveletVectorPlot],
            "BaseStyle"->BaseStyle/.{opts}/.Options[WaveletVectorPlot],
            "ColorOutput"->ColorOutput/.{opts}/.Options[WaveletVectorPlot],
            "ContentSelectable"->ContentSelectable/.{opts}/.Options[WaveletVectorPlot],
            "DisplayFunction":>DisplayFunction/.{opts}/.Options[WaveletVectorPlot],
            "Epilog"->Epilog/.{opts}/.Options[WaveletVectorPlot],
            "FormatType":>FormatType/.{opts}/.Options[WaveletVectorPlot],
            "Frame"->Frame/.{opts}/.Options[WaveletVectorPlot],
            "FrameLabel"->FrameLabel/.{opts}/.Options[WaveletVectorPlot],
            "FrameStyle"->FrameStyle/.{opts}/.Options[WaveletVectorPlot],
            "FrameTicks"->FrameTicks/.{opts}/.Options[WaveletVectorPlot],
            "FrameTicksStyle"->FrameTicksStyle/.{opts}/.Options[WaveletVectorPlot],
            "GridLines"->GridLines/.{opts}/.Options[WaveletVectorPlot],
            "GridLinesStyle"->GridLinesStyle/.{opts}/.Options[WaveletVectorPlot],
            "ImageMargins"->ImageMargins/.{opts}/.Options[WaveletVectorPlot],
            "ImagePadding"->ImagePadding/.{opts}/.Options[WaveletVectorPlot],
            "ImageSize"->ImageSize/.{opts}/.Options[WaveletVectorPlot],
            "LabelStyle"->LabelStyle/.{opts}/.Options[WaveletVectorPlot],
            "Method"->Method/.{opts}/.Options[WaveletVectorPlot],
            "PlotLabel"->PlotLabel/.{opts}/.Options[WaveletVectorPlot],
            "PlotRange"->PlotRange/.{opts}/.Options[WaveletVectorPlot],
            "PlotRangeClipping"->PlotRangeClipping/.{opts}/.Options[WaveletVectorPlot],
            "PlotRangePadding"->PlotRangePadding/.{opts}/.Options[WaveletVectorPlot],
            "PlotRegion"->PlotRegion/.{opts}/.Options[WaveletVectorPlot],
            "PreserveImageOptions"->PreserveImageOptions/.{opts}/.Options[WaveletVectorPlot],
            "Prolog"->Prolog/.{opts}/.Options[WaveletVectorPlot],
            "RotateLabel"->RotateLabel/.{opts}/.Options[WaveletVectorPlot],
            "Ticks"->Ticks/.{opts}/.Options[WaveletVectorPlot],
            "TicksStyle"->TicksStyle/.{opts}/.Options[WaveletVectorPlot]];
         Return[plt];
];

WaveletVectorPlay[v_,opts___]:=Module[{nits,iter,region,wtlist,sr,srlist,k,lbl,clips,snip,output},

    nits=NumIterations/.{opts}/.Options[WaveletVectorPlay];
    	If[nits===Max,nits=FactorInteger[Length[v]][[1,2]]];
      	
    iter=Iteration/.{opts}/.Options[WaveletVectorPlay];
      
    region=ToString[Region/.{opts}/.Options[WaveletVectorPlay]];
    wtlist=WaveletVectorToList[v,NumIterations->nits];        
    sr=SampleRate/.{opts}/.Options[WaveletVectorPlay];
srlist=Table[Round[sr/2^k],{k,1,nits}];
srlist=Prepend[Reverse[srlist],Last[srlist]];
     
    lbl=Table[Text[Style["Highpass Iteration "<>ToString[k],Gray,Plain,18]],{k,Length[wtlist]-1,1,-1}];
lbl=Prepend[lbl,Text[Style["Lowpass",Gray,Plain,18]]];

If[region=="All",
clips=Table[ ListPlay[wtlist[[k]],
    	"DisplayFunction"->DisplayFunction/.{opts}/.Options[WaveletVectorPlay],
        "Epilog"->Epilog/.{opts}/.Options[WaveletVectorPlay],
        "PlayRange"->PlayRange/.{opts}/.Options[WaveletVectorPlay],
        "Prolog"->Prolog/.{opts}/.Options[WaveletVectorPlay],
        "SampleDepth"->SampleDepth/.{opts}/.Options[WaveletVectorPlay],
      	"SampleRate"->srlist[[k]]],{k,1,Length[wtlist]}],
If[region=="LowPass",snip=First[wtlist];sr=First[srlist];lbl=First[lbl];];
If[region=="HighPass",snip=Reverse[Drop[wtlist,1]][[iter]];sr=Reverse[Drop[srlist,1]][[iter]];lbl=Reverse[Drop[lbl,1]][[iter]];];
clips=ListPlay[snip,
    	"DisplayFunction"->DisplayFunction/.{opts}/.Options[WaveletVectorPlay],
        "Epilog"->Epilog/.{opts}/.Options[WaveletVectorPlay],
        "PlayRange"->PlayRange/.{opts}/.Options[WaveletVectorPlay],
        "Prolog"->Prolog/.{opts}/.Options[WaveletVectorPlay],
        "SampleDepth"->SampleDepth/.{opts}/.Options[WaveletVectorPlay],
      	"SampleRate"->sr];
];
If[region=="All",output=TableForm[Flatten[Transpose[{clips,lbl}]]],output=TableForm[{clips,lbl}]];
Return[output];
];

HuffmanTree[codes_,opts___]:=Module[{l,p,b,nc,nec,ns,nfs,nfc,el,lvls,node,pr,lbl,k,idx,tree={},j,plt,ss,str={}},
	If[Dimensions[codes]!=2,Message[HuffmanTree::"baddimensions"];Return[0];];
	{l,p,b}=Transpose[codes];
	If[Length[l]==1,b={{0}}];
	If[Last[Dimensions[codes]]!=3,Message[HuffmanTree::"badinput"];Return[Dimensions[{l,p,b}]];];
	If[Length[Union[Map[Length,{l,p,b}]]]!=1,Message[HuffmanTree::"listlengths"];Return[0];];
	If[ToString[VectorQ[p,NumericQ]]!=True,Message[HuffmanTree::"probabilities"];Return[0];];
	If[Union[Map[MemberQ[{0,1},#]&,Union[Flatten[b]]]]!={True},Message[HuffmanTree::"binary"];Return[0];];
	nc=NodeColor/.{opts}/.Options[HuffmanTree];
	nec=NodeEdgeColor/.{opts}/.Options[HuffmanTree];
	ns=NodeSize/.{opts}/.Options[HuffmanTree];
	nfs=NodeFontSize/.{opts}/.Options[HuffmanTree];
	nfc=NodeFontColor/.{opts}/.Options[HuffmanTree];
	el=EdgeLength/.{opts}/.Options[HuffmanTree];
	ss=ShowString/.{opts}/.Options[HuffmanTree];
	str=Style[ss,Gray,Plain,24];
	p=SetPrecision[p,2];
	If[Length[l]==1,
Return[TreePlot[{1->1},Top,1,VertexLabeling->True,VertexRenderingFunction->({nc,EdgeForm[nec],Disk[#,ns],Text[Style[l[[1]],nfc,nfs],#1]}&),PlotLabel->str]]];
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
	plt=TreePlot[Sort[tree],Top,1,LayerSizeFunction-> (#^el&),VertexLabeling->True,VertexRenderingFunction->({nc,EdgeForm[nec],Disk[#,ns],Text[Style[lbl[[#2]],nfc,nfs],#1]}&),PlotLabel->str];
	Return[plt];
];


(* Messages *)

CE::zerolength="The input is not a list - cannot compute the cumulative energy."
Comp::zerolength="The input is not a list - returning an empty list."
Comp::noninteger="The second argument of Comp is not an integer - returning an empty list."
nCE::nonvector="The input is not a vector - returning 0."
Entropy::nonvectormatrix="The input must be a vector or a matrix - returning 0."
MSE::nonmatrix="At least one of the inputs is not a matrix - returning -1."
MSE::baddimensions="The dimensions of the input matrices are not the same - returning -1."
PSNR::badmatrix="Either one of the entries is not a matrix or the dimensions of the input matrices are not the same - returning 0."
Daub::ninteger="The argument `1` must be an integer - returning the filter {0,0}."
Daub::nposeven="The argument `1` is not a positive even integer - returning the filter {0,0}."
Coif::noninteger="The argument `1` must be either 1, 2, or 3 - returning the filter {0,0}."
Coif::badinteger="The argument `1` must be either 1, 2, or 3 - returning the filter {0,0}."
SplineFilters::badintegers="The input values must both be even integers - returning the filters {0,0} and {0,0}."
HWT1D1::evenlength="The length of the input vector v must be even - returning the input vector."
IHWT1D1::evenlength="The length of the input vector v must be even - returning the input vector."
HWT1D::nonvector="The input v is not a vector - returning the input."
HWT1D::oddlength="The input vector v is of odd length - returning the input."
HWT1D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
HWT1D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
IHWT1D::nonvector="The input v is not a vector - returning the input."
IHWT1D::oddlength="The input vector v is of odd length - returning the input."
IHWT1D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
IHWT1D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
FiniteFourier::badindex="The starting index must be an integer - returning 0."
MakeHuffmanCodes::badinput="The input values must be integer-valued - returning {0,0,0}."
WT1D1::evenlength="The length of the input vector v must be even - returning the input vector."
WT1D1::badfilterlength="The length of the filter must be even - returning the input vector."
IWT1D1::evenlength="The length of the input vector v must be even - returning the input vector."
IWT1D1::badfilterlength="The length of the filter must be even - returning the input vector."
WT1D::nonvector="The input v is not a vector - returning the input."
WT1D::oddlength="The input vector v is of odd length - returning the input."
WT1D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
WT1D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
IWT1D::nonvector="The input v is not a vector - returning the input."
IWT1D::oddlength="The input vector v is of odd length - returning the input."
IWT1D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
IWT1D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
BWT1D1::nonvector="The input v is not a vector - returning the input."
BWT1D1::oddlength="The input vector v is of odd length - returning the input."
BWT1D1::onefilter="Only one input filter detected - returning the discrete wavelet transform using that filter."
BWT1D1::nofilter="No filter detected - returning the Haar wavelet transform of the input vector."
BWT1D1::badfilterlengths="The lengths of the filters must either both be even or both be odd."
IBWT1D1::nonvector="The input v is not a vector - returning the input."
IBWT1D1::oddlength="The input vector v is of odd length - returning the input."
IBWT1D1::onefilter="Only one input filter detected - returning the discrete wavelet transform using that filter."
IBWT1D1::nofilter="No filter detected - returning the Haar wavelet transform of the input vector."
IBWT1D1::badfilterlengths="The lengths of the filters must either both be even or both be odd."
BWT1D::nonvector="The input v is not a vector - returning the input."
BWT1D::oddlength="The input vector v is of odd length - returning the input."
BWT1D::onefilter="Only one input filter detected - returning the discrete wavelet transform using that filter."
BWT1D::nofilter="No filter detected - returning the Haar wavelet transform of the input vector."
BWT1D::badfilterlengths="The lengths of the filters must either both be even or both be odd."
BWT1D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
BWT1D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
IBWT1D::nonvector="The input v is not a vector - returning the input."
IBWT1D::oddlength="The input vector v is of odd length - returning the input."
IBWT1D::onefilter="Only one input filter detected - returning the discrete wavelet transform using that filter."
IBWT1D::nofilter="No filter detected - returning the Haar wavelet transform of the input vector."
IBWT1D::badfilterlengths="The lengths of the filters must either both be even or both be odd."
IBWT1D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
IBWT1D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
LWT1D1::evenlength="The length of the input vector v must be even - returning the input vector."
LWT1D1::nonvector="The input v is not a vector - returning the input."
LWT1D1::noninteger="When IntegerMap is set to true, the input must be integer-valued.  Returing the input vector."
ILWT1D1::evenlength="The length of the input vector v must be even - returning the input vector."
ILWT1D1::nonvector="The input v is not a vector - returning the input."
ILWT1D1::noninteger="When IntegerMap is set to true, the input must be integer-valued.  Returing the input vector."
LWT1D::nonvector="The input v is not a vector - returning the input."
LWT1D::oddlength="The input vector v is of odd length - returning the input."
LWT1D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
LWT1D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
ILWT1D::nonvector="The input v is not a vector - returning the input."
ILWT1D::oddlength="The input vector v is of odd length - returning the input."
ILWT1D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
ILWT1D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
GetCorner::badinput="The input must be a numeric matrix.  Returning the input."
GetCorner::integerdimensions="The dimensions of the corner to select must be integers.  Returning the input."
GetCorner::baddimensions="One of the input dimensions `1`x`2` is larger than the dimensions `3`x`4` of the input matrix.  Returning the input."
PutCorner::badinput="The input must be a numeric matrix.  Returning the input."
PutCorner::integerdimensions="The dimensions of the corner to select must be integers.  Returning the input."
PutCorner::baddimensions="One of the input dimensions `1`x`2` is larger than the dimensions `3`x`4` of the input matrix.  Returning the input."
LeftHWT::badinput="The input must be a numeric matrix.  Returning the input."
LeftHWT::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
RightHWT::badinput="The input must be a numeric matrix.  Returning the input."
RightHWT::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
HWT2D1::badinput="The input must be a numeric matrix.  Returning the input."
HWT2D1::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
LeftIHWT::badinput="The input must be a numeric matrix.  Returning the input."
LeftIHWT::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
IHWT2D1::badinput="The input must be a numeric matrix.  Returning the input."
IHWT2D1::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
HWT2D::badinput="The input must be a numeric matrix.  Returning the input."
HWT2D::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
HWT2D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
HWT2D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
IHWT2D::badinput="The input must be a numeric matrix.  Returning the input."
IHWT2D::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
IHWT2D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
IHWT2D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
LeftWT::badinput="The input must be a numeric matrix.  Returning the input."
LeftWT::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
LeftWT::badfilterlength="The length of the filter must be even - returning the input matrix."
LeftWT::nofilter="No filter given - using the Haar filter for the computations."
RightWT::badinput="The input must be a numeric matrix.  Returning the input."
RightWT::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
RightWT::badfilterlength="The length of the filter must be even - returning the input matrix."
RightWT::nofilter="No filter given - using the Haar filter for the computations."
WT2D1::badinput="The input must be a numeric matrix.  Returning the input."
WT2D1::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
WT2D1::badfilterlength="The length of the filter must be even - returning the input matrix."
WT2D1::nofilter="No filter given - using the Haar filter to compute the transform."
LeftIWT::badinput="The input must be a numeric matrix.  Returning the input."
LeftIWT::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
LeftIWT::badfilterlength="The length of the filter must be even - returning the input matrix."
LeftIWT::nofilter="No filter given - using the Haar filter for the computations."
IWT2D1::badinput="The input must be a numeric matrix.  Returning the input."
IWT2D1::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
IWT2D1::badfilterlength="The length of the filter must be even - returning the input matrix."
IWT2D1::nofilter="No filter given - using the Haar filter to compute the inverse transform."
WT2D::badinput="The input must be a numeric matrix.  Returning the input."
WT2D::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
WT2D::badfilterlength="The length of the filter must be even - returning the input matrix."
WT2D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
WT2D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
WT2D::nofilter="No filter given - using the Haar filter to compute the transform."
IWT2D::badinput="The input must be a numeric matrix.  Returning the input."
IWT2D::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
IWT2D::badfilterlength="The length of the filter must be even - returning the input matrix."
IWT2D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
IWT2D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
IWT2D::nofilter="No filter given - using the Haar filter to compute the inverse transform."
BWT2D1::nofilter="No filter given - using the Haar filter to compute the inverse transform."
BWT2D1::badinput="The input must be a numeric matrix.  Returning the input."
BWT2D1::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
BWT2D1::onefilter="Only one input filter detected - returning the wavelet transform using that filter."
IBWT2D1::nofilter="No filter given - using the Haar filter to compute the inverse transform."
IBWT2D1::badinput="The input must be a numeric matrix.  Returning the input."
IBWT2D1::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
IBWT2D1::onefilter="Only one input filter detected - returning the inverse wavelet transform using that filter."
BWT2D::nofilter="No filter given - using the Haar filter to compute the transform."
BWT2D::badinput="The input must be a numeric matrix.  Returning the input."
BWT2D::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
BWT2D::onefilter="Only one input filter detected - returning the wavelet transform using that filter."
BWT2D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
BWT2D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
IBWT2D::nofilter="No filter given - using the Haar filter to compute the inverse transform."
IBWT2D::badinput="The input must be a numeric matrix.  Returning the input."
IBWT2D::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
IBWT2D::onefilter="Only one input filter detected - returning the inverse wavelet transform using that filter."
IBWT2D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
IBWT2D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
LWT2D1::badinput="The input must be a numeric matrix.  Returning the input."
LWT2D1::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
ILWT2D1::badinput="The input must be a numeric matrix.  Returning the input."
ILWT2D1::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
LWT2D::badinput="The input must be a numeric matrix.  Returning the input."
LWT2D::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
LWT2D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
LWT2D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
ILWT2D::badinput="The input must be a numeric matrix.  Returning the input."
ILWT2D::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
ILWT2D::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
ILWT2D::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
DCT1D::badinput="The input must be a numeric vector.  Returning the input."
IDCT1D::badinput="The input must be a numeric vector.  Returning the input."
DCT2D::badinput="The input must be a numeric matrix.  Returning the input."
IDCT2D::badinput="The input must be a numeric matrix.  Returning the input."
WaveletMatrixToList::badinput="The input must be a numeric matrix.  Returning the input."
WaveletMatrixToList::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
WaveletMatrixToList::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
WaveletMatrixToList::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
WaveletListToMatrix::badinput="The dimensions of each portion of the transform is not in the correct form.  Returning the input."
WaveletListToMatrix::baddimensions="The input is a matrix with at least one odd dimension.  Returning the input."
WaveletListToMatrix::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
WaveletListToMatrix::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
GammaCorrection::badinput="The input must be a numeric matrix consisting of integers.  Returning the input."
GammaCorrection::badvalue="The correction constant must be a numerical value.  Returning the input."
GammaCorrection::positivevalue="The correction constant must be positive.  Returning the input."
GammaCorrection::badrange="The input matrix must be comprised of integers in [0, 255]."
MakeHistogramEQ::badinput="The input must be a numeric matrix consisting of integers.  Returning the input."
MakeHistogramEQ::badrange="The input matrix must be comprised of integers in [0, 255]."
HistogramEQ::badinput="The input must be a numeric matrix consisting of integers.  Returning the input."
HistogramEQ::badrange="The input matrix must be comprised of integers in [0, 255]."
WaveletVectorToList::nonvector="The input v is not a vector - returning the input."
WaveletVectorToList::oddlength="The input vector v is of odd length - returning the input."
WaveletVectorToList::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
WaveletListToVector::nonvector="The input v is not a vector - returning the input."
ChopVector::maxnt="The power of two `1` exceeds the length of the vector `2` - returning the original input."

MaximumIterations::badinput="The input is not a numeric matrix or list of three numeric matrices - returning 0."
MaximumIterations::baddimensions="One of the dimensions of an input matrix is odd - returning 0"
ImageRead::badinput="The input file `1` was not found or is not a valid image.  Returning the 2x2 zero matrix."
ImageRead::maxnt="2^`1` exceeds at least one of the image dimensions (`2`,`3`)."
ImagePlot::baddimensions="One of the dimensions of an input matrix is odd.  WaveletDensityPlot failed."
ImagePlot::differentdimensions="The three input matrices do not have the same dimensions.  WaveletDensityPlot failed."
CreateImageObject::baddimensions="One of the dimensions of an input matrix is odd.  WaveletDensityPlot failed."
CreateImageObject::differentdimensions="The three input matrices do not have the same dimensions.  WaveletDensityPlot failed."
CreateImageObject::onlyRGB="You must use RGBColor or a color from AllColors for the ChannelColor.  If you want to use this  directive, use it with ColorFunction and omit ChannelColor."
CreateImageObject::LeftWT="The number of rows in the input matrix is odd.  Setting LinearScaling to False."
CreateImageObject::RightWT="The number of columns in the input matrix is odd.  Setting LinearScaling to False."
WaveletDensityPlot::badinput="The input is not a numeric matrix or list of three numeric matrices.  WaveletDensityPlot failed."
WaveletDensityPlot::baddimensions="One of the dimensions of an input matrix is odd.  WaveletDensityPlot failed."
WaveletDensityPlot::differentdimensions="The three input matrices do not have the same dimensions.  WaveletDensityPlot failed."
WaveletDensityPlot::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
WaveletDensityPlot::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
ExtractRegion::badinput="The input is not a numeric matrix or list of three numeric matrices.  WaveletDensityPlot failed."
ExtractRegion::baddimensions="One of the dimensions of an input matrix is odd.  WaveletDensityPlot failed."
ExtractRegion::differentdimensions="The three input matrices do not have the same dimensions.  WaveletDensityPlot failed."
ExtractRegion::maxits="The NumIterations `1` is greater than the maximum number of iterations that can be performed on the vector - setting the its = `2`."
ExtractRegion::badits="NumIterations value must be an integer between 0 and `1` inclusive - setting NumIterations to 1."
HuffmanTree::baddimensions="Bad input length - dimensions longer than 2."
HuffmanTree::badinput="Bad input length - second dimension is not 3."
HuffmanTree::listlengths="List lengths are different."
HuffmanTree::probabilities="The probability list does not consist of numeric values."
HuffmanTree::binary="The code list does not consist of 0s and 1s."
TestSparseness::badinput="Input must be either a vector or a list of vectors comprised of numerical values."
NoiseEstimate::badinput="Input must be either a matrix or vector of numerical values."
UniversalThreshold::badinput="Input must be either a matrix or vector of numerical values."
WaveletShrinkage::badinput="Input must be either a matrix or vector of numerical values."
WaveletShrinage::numiterations="NumIterations is a larger value than allowed."
WaveletShrinkage::lambdapositive="lambda must be a nonnegative number - returing the input."
WaveletShrinkage::lambdavectorpositive="lambda must consist of nonnegative numbers - returing the input."
WaveletShrinkage::lambdalength="lambda must be a list containing `1` nonnegative numbers - returning input."
WaveletShrinkage::lambdamatrix="lambda must be a matrix of dimension `1` x 3 consisting of nonnegative numbers - returning input."
SureShrink::badinput="Input must be either a matrix or vector of numerical values."
SureShrink::maxits="NumIterations is a larger value than allowed."
End[]

EndPackage[];

