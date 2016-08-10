(* ::Package:: *)

(***************************************************************************************)
(*Article: Boundary crossing of Oreder Statistic Point Processes;
Title: BoundaryCrossingOSPPFunctionV2.nb;
Author: Pierre-Olivier Goffard
*)
(* Notations *)
(*
- For all the point processes;
\[Nu] = Time transformation;
\[Mu] = Mean of the Point process

- For Mixed Poisson Process;
W = Distribution of the intensity of the mixed Poisson process(possibly degenerate in the case of the homogeneous Poisson process for instance);
\[Lambda] = Intensity for an homogeneous Poisson process;
m = Mean parameter of the gamma prior distribution associated to a Polya-Lundberg process;
r = Shape parameter of the gamma prior distribution associated to a Polya--Lunndberg process;

- For mixed sample process;
Z = Distribution of the size of the initial population (degenerate in the case of the linear death process);
\[Gamma] = Limit of the expectation of the point process when t goes to infinity;

- For all boundaries;
h = Shape of the boundary;
Inverseh = generalized inverse of the shape function of the boundary;

- For the lower boundary;
\[Alpha] = The lower boundary starts at the point (0,-\[Alpha]);
n = Level at which the meeting occurs (Continuous crossing, note that it is equivalent to specify the time Inverseh[n+L] at which the crossing happens);

- For the upper boundary;
\[Beta] = The upper boundary starts at the point (0,\[Beta]);
t = Time horizon for which crossing (non continuous because it happens during a jump) has not occured yet; 

- Miscelaneous;
G0[n] = Abel-Gontcharov polynomials of degree n evaluated at 0;
A0[n] = Appell polynomials of degree n evaluated at 0;
A1[n] = Appell polynomials of degree n evaluated at 1;
d[n] = Sheffer type polynomials for rectangle probabilities;
\[Alpha][n] = Instant at which the lower boundary reach n (intger valued);
\[Beta][n] = Instant at which the upper boundary reach n (intger valued);
F_t = Disribution function of the jumping time oof the OSPP between 0 and t;
MP = Mixed Poisson process;
MPTT = Mixed Poisson process with a Time Transformation;
MSP = Mixed Sample Process
*)
PDFMPTT[W_,\[Nu]_,n_,t_]:=(D[MomentGeneratingFunction[W,\[Nu][t]*(z-1)],{z,n}]/.z->0)/n!;
(*PDF of the continuous crossing of a Mixed Poisson process with a time transformation and a Lower boundary*)
PDFMeetingLevelMPTTLowerBoundary[W_,\[Nu]_,n_,Inverseh_,\[Alpha]_]:=Module[{\[Alpha]n,F\[Alpha]n,G0,z},
\[Alpha]n=Function[{k},Inverseh[k+\[Alpha]]];
F\[Alpha]n=Function[{k},\[Nu][\[Alpha]n[k]]/\[Nu][\[Alpha]n[n]]];
G0[0]=1;
Do[
G0[k]=-Sum[
Binomial[k,j]*(F\[Alpha]n[j])^(k-j)*G0[j]
,{j,0,k-1}];
,{k,1,n}];
(D[MomentGeneratingFunction[W,\[Nu][\[Alpha]n[n]]*(z-1)],{z,n}]/.z->0)/n!*(-1)^(n)*G0[n]
];
(*PDF of the meeting time of a Mixed Poisson process and a Lower linear boundary*)
PDFMeetingLevelMPLowerLinearBoundary[W_,n_,c_,\[Alpha]_]:=Module[{\[Alpha]n,Fln,G0,z},
\[Alpha]n=Function[{k},(k+\[Alpha])/c];
Fln=Function[{k},\[Alpha]n[k]/\[Alpha]n[n]];
(D[MomentGeneratingFunction[W,\[Alpha]n[n]*(z-1)],{z,n}]/.z->0)/n!*Fln[\[Alpha]n[0]];
];
(*PDF of a Mixed Sample process*)
PDFMSP[Z_,\[Mu]_,\[Gamma]_,n_,t_]:=(D[MomentGeneratingFunction[Z,Log[1-\[Mu][t]/\[Gamma]+\[Mu][t]/\[Gamma]*z]],{z,n}]/.z->0)/n!;
(*PDF of the continuous crossing of a Mixed Sample process and a Lower boundary*)
PDFMeetingLevelMSPLowerBoundary[Z_,\[Mu]_,\[Gamma]_,n_,Inverseh_,\[Alpha]_]:=Module[{\[Alpha]n,z,Fln,G0},
(*Function that gives the times at which the lower boundary is integer valued*)
\[Alpha]n=Function[{k},Inverseh[k+\[Alpha]]];
(*Distribution Function of the sample from which the order statistics are computed*) 
Fln=Function[{k},\[Mu][\[Alpha]n[k]]/\[Mu][\[Alpha]n[n]]];
G0[0]=1;
Do[
G0[k]=-Sum[
Binomial[k,j]*(Fln[j])^(k-j)*G0[j]
,{j,0,k-1}];
,{k,1,n}];
(D[MomentGeneratingFunction[Z,Log[1-\[Mu][\[Alpha]n[n]]/\[Gamma]+\[Mu][\[Alpha]n[n]]/\[Gamma]*z]],{z,n}]/.z->0)/n!*(-1)^(n)*G0[n]
];
(*Survival function of the crossing time of a Mixed Poisson process with a time transformation with an upper boundary*)
SurvivalCrossingTimeMPTTUpperBoundary[W_,\[Nu]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{\[Beta]n,Ft,A0,A1},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
Ft=Function[{k},\[Nu][\[Beta]n[k]]/\[Nu][t]];
A0[0]=1;A1[0]=1;
Do[
A0[k]=-Sum[
Binomial[k,j]*(Ft[k])^(j)*A0[k-j]
,{j,1,k}];
A1[k]=Sum[
Binomial[k,j]*A0[k-j]
,{j,0,k}];
,{k,1,Floor[\[Beta]+h[t]]}];
Sum[
(D[MomentGeneratingFunction[W,\[Nu][t]*(z-1)],{z,n}]/.z->0)/n!*A1[n]
,{n,0,Floor[h[t]+\[Beta]]}]
];
(*Survival function of the crossing time of a Mixed Poisson process and a linear upper boundary where U=0*)
SurvivalCrossingTimeMPLinearUpperBoundaryU0[W_,t_,c_]:=Module[{z},
Sum[
(1-n/c/t)(D[MomentGeneratingFunction[W,t*(z-1)],{z,n}]/.z->0)/n!
,{n,0,Floor[c*t]}];
];
(*Survival function of the crossing time of a Mixed Sample Process and an upper boundary*)
SurvivalFunctionCrossingTimeMSPUpperBoundary[Z_,\[Mu]_,\[Gamma]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{\[Beta]n,Ft,A0,A1,z},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
Ft=Function[{s},\[Mu][Min[s,t]]/\[Mu][t]];
A0[0]=1;A1[0]=1;
Do[
A0[k]=-Sum[
Binomial[k,j]*(Ft[\[Beta]n[k]])^(j)*A0[k-j]
,{j,1,k}];
A1[k]=Sum[
Binomial[k,j]*A0[k-j]
,{j,0,k}];
,{k,1,Floor[\[Beta]+h[t]]}];
Sum[(D[MomentGeneratingFunction[Z,Log[1-\[Mu][t]/\[Gamma]+\[Mu][t]/\[Gamma]*z]],{z,n}]/.z->0)/n!*A1[n],{n,0,Floor[\[Beta]+h[t]]}]
];
(*Two sided probabilities for the Mixed Poisson process with a time transformation*)
(*Probability of exiting through the lower boundary at level n without having crossed the upper boundary*)
PDFTwoSidedMPTTExitByLowerBoundary[W_,\[Nu]_,n_,Inverseh_,\[Alpha]_,\[Beta]_]:=Module[{\[Alpha]n,\[Beta]n,Fln,dn,z},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
\[Alpha]n=Function[{k},Inverseh[k+\[Alpha]]];
Fln=Function[{s},\[Nu][s]/\[Nu][\[Alpha]n[n]]];
dn[0]=1; 
Do[
dn[k]=-Sum[
Binomial[k,j]*(Max[Fln[\[Alpha]n[j]]-Fln[\[Beta]n[k]],0])^(k-j)*(-1)^(k-j)*dn[j]
,{j,0,k-1}];
,{k,1,n}];
(D[MomentGeneratingFunction[W,\[Nu][\[Alpha]n[n]]*(z-1)],{z,n}]/.z->0)/n!*dn[n]
];
(*Probability of not exiting the region between the lower boundary and the upper boundary before time t*)
SurvivalTwoSidedMPTTNotExiting[W_,\[Nu]_,t_,h_,Inverseh_,\[Alpha]_,\[Beta]_]:=Module[{\[Alpha]n,\[Beta]n,Ft,dn,z},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
\[Alpha]n=Function[{k},Inverseh[k+\[Alpha]]];
Ft=Function[{s},\[Nu][Min[s,t]]/\[Nu][t]];
dn[0]=1;
Do[
dn[k]=-Sum[
Binomial[k,j]*(Max[Ft[\[Alpha]n[j]]-Ft[\[Beta]n[k]],0])^(k-j)*(-1)^(k-j)*dn[j]
,{j,0,k-1}];
,{k,1,Floor[h[t]+\[Beta]]}];
Sum[
(D[MomentGeneratingFunction[W,\[Nu][t]*(z-1)],{z,n}]/.z->0)/n!*dn[n]
,{n,Max[Floor[h[t]-\[Alpha]],0],Floor[h[t]+\[Beta]]}]
];
(*Two sided probabilities for the Mixed Sample Process*)
(*Probability of exiting through the lower boundary at level n without having crossed the upper boundary*)
PDFTwoSidedMSPExitByLowerBoundary[Z_,\[Mu]_,\[Gamma]_,n_,Inverseh_,\[Alpha]_,\[Beta]_]:=Module[{\[Alpha]n,\[Beta]n,Fln,dn,z},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
\[Alpha]n=Function[{k},Inverseh[k+\[Alpha]]];
Fln=Function[{s},\[Mu][s]/\[Mu][\[Alpha]n[n]]];
dn[0]=1;
Do[
dn[k]=-Sum[
Binomial[k,j]*(Max[Fln[\[Alpha]n[j]]-Fln[\[Beta]n[k]],0])^(k-j)*(-1)^(k-j)*dn[j]
,{j,0,k-1}];
,{k,1,n}];
(D[MomentGeneratingFunction[Z,Log[1-\[Mu][\[Alpha]n[n]]/\[Gamma]+\[Mu][\[Alpha]n[n]]/\[Gamma]*z]],{z,n}]/.z->0)/n!*dn[n]

];
(*Probability of not exiting the region between the lower boundary and the upper boundary before time t*)
SurvivalTwoSidedMPSNotExiting[Z_,\[Mu]_,\[Gamma]_,t_,h_,Inverseh_,\[Alpha]_,\[Beta]_]:=Module[{\[Alpha]n,\[Beta]n,Ft,dn,z},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
\[Alpha]n=Function[{k},Inverseh[k+\[Alpha]]];
Ft=Function[{s},\[Mu][Min[s,t]]/\[Mu][t]];
dn[0]=1;
Do[
dn[k]=-Sum[
Binomial[k,j]*(Max[Ft[\[Alpha]n[j]]-Ft[\[Beta]n[k]],0])^(k-j)*(-1)^(k-j)*dn[j]
,{j,0,k-1}];
,{k,1,Floor[h[t]+\[Beta]]}];
Sum[
(D[MomentGeneratingFunction[Z,Log[1-\[Mu][t]/\[Gamma]+\[Mu][t]/\[Gamma]*z]],{z,n}]/.z->0)/n!*dn[n]
,{n,Max[Floor[h[t]-\[Alpha]],0],Floor[h[t]+\[Beta]]}]
];


(************************************************************)
(*Complementary Functions*)
(*Plot the PDF, Survival Functions and CDF of discrete distribution*)
PlotDiscreteDistribution[Expression_,LabelSize_,Image_,n_]:=DiscretePlot[Expression[k],{k,0,n},PlotRange->All,ExtentSize->Right,PlotMarkers->Point,AxesOrigin->{0,0},LabelStyle->Directive[Black,LabelSize],ImageSize->Image,PlotStyle->Blue,FillingStyle->Blue];
TableDiscreteDistribution[Expression_,n_]:=Grid[Table[{k,N[Expression[k]]},{k,0,n}]];
(*Function to produce a sample of R replication of a general Mixed Poisson process*)
SampleMPTT[W_,\[Nu]_,t_,R_]:=Table[RandomVariate[PoissonDistribution[RandomVariate[W]*\[Nu][t]]],{r,1,R}]
(*A trajectory is a sequence of jump times T0,T1,...,TNt, so an instance of Nt and the associated Jump times*)
SampleMPTTTrajectory[W_,\[Nu]_,Inverse\[Nu]_,t_,R_]:=Module[{Nt,JumpTimes,Traj},
Traj={};
Do[
Nt=SampleMPTT[W,\[Nu],t,1];
JumpTimes=Sort[Map[Inverse\[Nu],RandomVariate[UniformDistribution[{0,1}],Nt]]];
(*We add T0=0*)
JumpTimes=PrependTo[JumpTimes,0];
Traj=AppendTo[Traj,JumpTimes]
,{k,1,R}];
Traj
]
(*Function to plot 5 trajectories drawn from a mixed Poisson process with a time transformation and upper boundary*)
PlotTrajectoriesMPTTUpperBoundary[Trajectories_,h_,\[Beta]_,LabelSize_,PicSize_]:=Module[{ArgPiecewise,TrajFunc,TrajPiecewiseFunc,Nt},
TrajFunc={};
Do[
Nt=Length[Trajectories[[r]]]-1;
ArgPiecewise=Table[
{k,Trajectories[[r,k+1]]<=x<=Trajectories[[r,k+2]]},{k,0,Nt-1}];
ArgPiecewise=AppendTo[ArgPiecewise,
{Nt,Trajectories[[r,Nt+1]]<=x<=t}];
TrajFunc=AppendTo[TrajFunc,ArgPiecewise];
,{r,1,Length[Trajectories]}];
TrajPiecewiseFunc=Map[Piecewise,TrajFunc];
TrajPiecewiseFunc=AppendTo[TrajPiecewiseFunc,h[x]+\[Beta]];
Plot[TrajPiecewiseFunc,{x,0,t},PlotRange->All,ImageSize->PicSize,LabelStyle->Directive[Black,LabelSize],Filling->Table[{r->Axis},{r,1,5}],
PlotStyle->{Directive[Dashed,Blue],Directive[Dashed,Blue],Directive[Dashed,Blue],Directive[Dashed,Blue],Directive[Dashed,Blue],Red},
FrameLabel->{"t","","",""},
Frame->True,
FillingStyle->Directive[LightBlue,Opacity[0.1]]]
];
(*Function giving for each trajectory the probability of crossing with CMC and AppellMC*)
CMCvsAppellMCMPTTUpperBoundary[Trajectories_,\[Nu]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{\[Beta]n,Ft,Nts,A0,A1,CMC,Traj,AppellMC},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
Ft=Function[{s},\[Nu][Min[s,t]]/\[Nu][t]];
Nts=Table[Min[Length[Trajectories[[k]]]-1,Floor[\[Beta]+h[t]]+1]
,{k,1,Length[Trajectories]}];
A0[0]=1;
A1[0]=1;
A1[Floor[\[Beta]+h[t]]+1]=0;
Do[A0[k]=-Sum[Binomial[k,j]*(Ft[\[Beta]n[k]])^(j)*A0[k-j],{j,1,k}];
A1[k]=Sum[Binomial[k,j]*A0[k-j],{j,0,k}];
,{k,1,Min[Floor[\[Beta]+h[t]]+1,Max[Nts]]}];
AppellMC=Map[A1,Nts];
CMC={};
Do[
If[Nts[[r]]>h[t]+\[Beta],AppendTo[CMC,0],
Traj=Trajectories[[r]];
If[
Count[Traj-Table[\[Beta]n[i],{i,0,Nts[[r]]}],n_/;n<0]==0
,AppendTo[CMC,1],AppendTo[CMC,0]]]
,{r,1,Length[Trajectories]}];
{AppellMC,CMC}
];
(*Variance Of the CMC estimator*)
VarianceCMCCrossingTimeMPTTUpperBoundary[W_,\[Nu]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{\[Beta]n,Ft,A0,A1,p},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
Ft=Function[{k},\[Nu][\[Beta]n[k]]/\[Nu][t]];
A0[0]=1;A1[0]=1;
Do[
A0[k]=-Sum[
Binomial[k,j]*(Ft[k])^(j)*A0[k-j]
,{j,1,k}];
A1[k]=Sum[
Binomial[k,j]*A0[k-j]
,{j,0,k}];
,{k,1,Floor[\[Beta]+h[t]]}];
p=Sum[
(D[MomentGeneratingFunction[W,\[Nu][t]*(z-1)],{z,n}]/.z->0)/n!*A1[n]
,{n,0,Floor[h[t]+\[Beta]]}];
p(1-p)
];
(*Variance Of the APMC estimator for mixed Poisson process with a time transformation*)
VarianceAPMCCrossingTimeMPTTUpperBoundary[W_,\[Nu]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{\[Beta]n,Ft,A0,A1,p},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
Ft=Function[{k},\[Nu][\[Beta]n[k]]/\[Nu][t]];
A0[0]=1;A1[0]=1;
Do[
A0[k]=-Sum[
Binomial[k,j]*(Ft[k])^(j)*A0[k-j]
,{j,1,k}];
A1[k]=Sum[
Binomial[k,j]*A0[k-j]
,{j,0,k}];
,{k,1,Floor[\[Beta]+h[t]]}];
p=Sum[
(D[MomentGeneratingFunction[W,\[Nu][t]*(z-1)],{z,n}]/.z->0)/n!*A1[n]
,{n,0,Floor[h[t]+\[Beta]]}];
Sum[
(D[MomentGeneratingFunction[W,\[Nu][t]*(z-1)],{z,n}]/.z->0)/n!*(A1[n])^(2)
,{n,0,Floor[h[t]+\[Beta]]}]-p^(2)
];
(*Difference Between the Variance of the two Monte Carlo estimators*)
DeltaVarianceCMCAPMCCrossingTimeMPTTUpperBoundary[W_,\[Nu]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{\[Beta]n,Ft,A0,A1,p,VarCMC,VarAppellMC},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
Ft=Function[{k},\[Nu][\[Beta]n[k]]/\[Nu][t]];
A0[0]=1;A1[0]=1;
Do[
A0[k]=-Sum[
Binomial[k,j]*(Ft[k])^(j)*A0[k-j]
,{j,1,k}];
A1[k]=Sum[
Binomial[k,j]*A0[k-j]
,{j,0,k}];
,{k,1,Floor[\[Beta]+h[t]]}];
p=Sum[
(D[MomentGeneratingFunction[W,\[Nu][t]*(z-1)],{z,n}]/.z->0)/n!*A1[n]
,{n,0,Floor[h[t]+\[Beta]]}];
VarCMC=p(1-p);
VarAppellMC=Sum[
(D[MomentGeneratingFunction[W,\[Nu][t]*(z-1)],{z,n}]/.z->0)/n!*(A1[n])^(2)
,{n,0,Floor[h[t]+\[Beta]]}]-p^(2);
(VarAppellMC-VarCMC)/VarCMC*100

];
(*Function to plot the variances of the Monte Carlo estimators for mixed Poisson process with a time transformation*)
PlotVarianceCMCAPMCrossingTimeMPTTUpperBoundary[W_,\[Nu]_,tMax_,h_,Inverseh_,\[Beta]_,SizeFont_,PicSize_]:=GraphicsGrid[
{{
Plot[{
VarianceCMCCrossingTimeMPTTUpperBoundary[W,\[Nu],t,h,Inverseh,\[Beta]],
VarianceAPMCCrossingTimeMPTTUpperBoundary[W,\[Nu],t,h,Inverseh,\[Beta]]
},{t,0,tMax},PlotRange->All,Filling->None,LabelStyle->{Black,SizeFont},AxesLabel->{"t",""},PlotLegends->Placed[{"\!\(\*SubsuperscriptBox[\(\[Sigma]\), \(CMC\), \(2\)]\)","\!\(\*SubsuperscriptBox[\(\[Sigma]\), \(APMC\), \(2\)]\)"},Above],PlotStyle->{{Red,Dashed},{Blue}},ImageSize->PicSize],
Plot[
DeltaVarianceCMCAPMCCrossingTimeMPTTUpperBoundary[W,\[Nu],t,h,Inverseh,\[Beta]]
,{t,0,tMax},PlotRange->All,Filling->None,LabelStyle->{Black,SizeFont},AxesLabel->{"t",""},PlotLegends->Placed[{"\!\(\*SuperscriptBox[\(\[CapitalDelta]\[Sigma]\), \(2\)]\)"},Above],PlotStyle->Blue,ImageSize->PicSize]
}},ImageSize->Full];
(*Functioon to analyse the variance reduction in the case of a mixed salple process*)
(*Variance Of the CMC estimator*)
VarianceCMCCrossingTimeMSPUpperBoundary[Z_,\[Mu]_,\[Gamma]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{p},
p=SurvivalFunctionCrossingTimeMSPUpperBoundary[Z,\[Mu],\[Gamma],t,h,Inverseh,\[Beta]];
p(1-p)
];
(*Variance of the APMC estimator*)
VarianceAPMCCrossingTimeMSPUpperBoundary[Z_,\[Mu]_,\[Gamma]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{\[Beta]n,Ft,A0,A1,p},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
Ft=Function[{s},\[Mu][Min[s,t]]/\[Mu][t]];
A0[0]=1;A1[0]=1;
Do[
A0[k]=-Sum[
Binomial[k,j]*(Ft[\[Beta]n[k]])^(j)*A0[k-j]
,{j,1,k}];
A1[k]=Sum[
Binomial[k,j]*A0[k-j]
,{j,0,k}];
,{k,1,Floor[\[Beta]+h[t]]}];
p=Sum[(D[MomentGeneratingFunction[Z,Log[1-\[Mu][t]/\[Gamma]+\[Mu][t]/\[Gamma]*z]],{z,n}]/.z->0)/n!*A1[n],{n,0,Floor[\[Beta]+h[t]]}];
Sum[
(D[MomentGeneratingFunction[Z,Log[1-\[Mu][t]/\[Gamma]+\[Mu][t]/\[Gamma]*z]],{z,n}]/.z->0)/n!*(A1[n])^(2)
,{n,0,Floor[h[t]+\[Beta]]}]-p^(2)
];
(*Relative Difference between the variances of the two estimators*)
DeltaVarianceCMCAPMCCrossingTimeMSPUpperBoundary[Z_,\[Mu]_,\[Gamma]_,t_,h_,Inverseh_,\[Beta]_]:=Module[{\[Beta]n,Ft,A0,A1,p,VarCMC,VarAppellMC},
\[Beta]n=Function[{k},Inverseh[Max[k-\[Beta],0]]];
Ft=Function[{k},\[Mu][\[Beta]n[k]]/\[Mu][t]];
A0[0]=1;A1[0]=1;
Do[
A0[k]=-Sum[
Binomial[k,j]*(Ft[k])^(j)*A0[k-j]
,{j,1,k}];
A1[k]=Sum[
Binomial[k,j]*A0[k-j]
,{j,0,k}];
,{k,1,Floor[\[Beta]+h[t]]}];
p=Sum[(D[MomentGeneratingFunction[Z,Log[1-\[Mu][t]/\[Gamma]+\[Mu][t]/\[Gamma]*z]],{z,n}]/.z->0)/n!*A1[n],{n,0,Floor[\[Beta]+h[t]]}];VarCMC=p(1-p);
VarAppellMC=Sum[
(D[MomentGeneratingFunction[Z,Log[1-\[Mu][t]/\[Gamma]+\[Mu][t]/\[Gamma]*z]],{z,n}]/.z->0)/n!*(A1[n])^(2)
,{n,0,Floor[h[t]+\[Beta]]}]-p^(2);
(VarAppellMC-VarCMC)/VarCMC*100

];
(*Function to plot the variances of the Monte Carlo estimators for mixed sample process *)
PlotVarianceCMCAPMCrossingTimeMSPUpperBoundary[Z_,\[Mu]_,\[Gamma]_,tMax_,h_,Inverseh_,\[Beta]_,SizeFont_,PicSize_]:=GraphicsGrid[
{{
Plot[{
VarianceCMCCrossingTimeMSPUpperBoundary[Z,\[Mu],\[Gamma],t,h,Inverseh,\[Beta]],
VarianceAPMCCrossingTimeMSPUpperBoundary[Z,\[Mu],\[Gamma],t,h,Inverseh,\[Beta]]
},{t,0,tMax},PlotRange->All,Filling->None,LabelStyle->{Black,SizeFont},AxesLabel->{"t",""},PlotLegends->Placed[{"\!\(\*SubsuperscriptBox[\(\[Sigma]\), \(CMC\), \(2\)]\)","\!\(\*SubsuperscriptBox[\(\[Sigma]\), \(APMC\), \(2\)]\)"},Above],PlotStyle->{{Red,Dashed},{Blue}},ImageSize->PicSize],
Plot[
DeltaVarianceCMCAPMCCrossingTimeMSPUpperBoundary[Z,\[Mu],\[Gamma],t,h,Inverseh,\[Beta]]
,{t,0,tMax},PlotRange->All,Filling->None,LabelStyle->{Black,SizeFont},AxesLabel->{"t",""},PlotLegends->Placed[{"\!\(\*SuperscriptBox[\(\[CapitalDelta]\[Sigma]\), \(2\)]\)"},Above],PlotStyle->Blue,ImageSize->PicSize]
}},ImageSize->Full];




