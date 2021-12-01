#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Input:: *)
(*(*Solve[{r==b-d,AI/PI==d/b},{b,d}]//FullSimplify*)*)
(*(*beta[x_, \[Eta]_] := 1-1/(\[Eta] + E^(Abs[x-2]^0.5));*)*)
(*beta[x_, \[Eta]_] := \[Eta];*)
(*model="1";*)
(*n=8;*)
(*\[Eta]=0.5;*)
(*\[Lambda]=1;*)
(*\[Delta]=0.375;*)
(*beta[2, \[Eta]]*)
(*test=Table[If[i==j, \[Lambda](2(1-beta[j, \[Eta]])^i-1)-\[Delta],*)
(*		If[2j>=i>j,\[Lambda] Binomial[j,i-j]beta[j, \[Eta]]^(i-j) (1-beta[j, \[Eta]])^(2j-i),If[j>i,\[Lambda] Binomial[j,j-i]beta[j, \[Eta]]^(j-i) (1-beta[j, \[Eta]])^i,0]]],{i,n},{j,n}];*)
(*Export[NotebookDirectory[]<>"test_model" <>model<>".txt",test];*)
(*test*)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ##;
## ## ## #VariableMSRates ## ## ## ##;
## ## ## ## ## ## ## ## ## ## ## ## ## ## ##;


(* ::Input:: *)
(*(*constant beta*)*)
(*beta[x_, \[Eta]_] := \[Eta]*)
(*var="\[Beta]";*)
(*model="constant";*)
(*vars=Table[\[Eta],{\[Eta],0.001,1,0.05}];*)


(* ::Input:: *)
(*(*vary eta -- linear form*)*)
(*var="\[Eta]";*)
(*model="linear";*)
(*c=0.0000001;*)
(*beta[x_, \[Eta]_] :=Min[0.6, \[Eta]* x+c];*)
(*vars=Table[\[Eta],{\[Eta],0.0001,0.12,0.005}];*)


(* ::Input:: *)
(*(*vary eta*)*)
(*var="\[Eta]";*)
(*model="sigmoid";*)
(*\[Theta]=0.8;*)
(*Clear[\[Eta]]*)
(*(*beta[x_, \[Eta]_] := ((Sin[x*3.2-1.7]+\[Theta]*x)/6)^(\[Eta]/x)*)*)
(*beta[x_, \[Eta]_] := 1-1/(\[Eta] + E^(Abs[x-2]^0.5))*)
(*vars=Table[\[Eta],{\[Eta],0.05,5,0.5}];*)


(* ::Input:: *)
(*(*vary theta*)*)
(*var="\[Theta]";*)
(*Clear[\[Theta]];*)
(*model="sin";*)
(*\[Eta] =2;*)
(*beta[x_, \[Theta]_] := 0.25*((Sin[x*3.2-1.7]+\[Theta]*x)/6)^(\[Eta]/x)*)
(*vars=Table[\[Theta],{\[Theta],0.6,2,0.1}];*)
(**)


(* ::Input:: *)
(*xy=Table[{x,beta[x,#]},{x,1,8,0.05}]&/@vars;*)
(*betaplot=ListPlot[#&/@xy,*)
(*Joined->True,*)
(*Frame->True,*)
(*Axes->False,*)
(*PlotRange->{All,All},*)
(*FrameLabel->{"ploidy","Mis-segregation rate (\[Beta])"},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*PlotStyle->ColorData[10,"ColorList"],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*PlotLegends->vars,*)
(*ImageSize->500,*)
(*AspectRatio->1,Epilog->Text[Style[beta[x,\[Eta]],FontSize->16],{5.5,0.1}]*)
(*]*)
(*Export[NotebookDirectory[]<>"betaplot_"<>model<>"_"<>var<>".jpg",betaplot]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*myfun[\[Lambda]_,\[Delta]_?NumberQ,n_?NumberQ, \[Eta]_ ?NumberQ]:=Eigenvalues[Table[If[i==j, \[Lambda](2(1-beta[j, \[Eta]])^i-1)-\[Delta],If[2j>=i>j,\[Lambda] Binomial[j,i-j]beta[j, \[Eta]]^(i-j) (1-beta[j, \[Eta]])^(2j-i),If[j>i,\[Lambda] Binomial[j,j-i]beta[j, \[Eta]]^(j-i) (1-beta[j, \[Eta]])^i,0]]],{i,n},{j,n}],*)
(*1,Method->{"Arnoldi","Criteria"->"RealPart"}][[1]]*)


(* ::Input:: *)
(*nvals={8};*)
(*ploidies={2};*)
(*leg="K = "<>ToString[#]&/@nvals;*)
(*leg2="p = "<>ToString[#]&/@ploidies;*)
(*colorList=ColorData[37,"ColorList"];*)
(*(*vary beta:*)*)
(*(*critTable=Table[{\[Delta],beta[2,critDelta[\[Delta],#]]},{\[Delta],0.0001,1,0.1}]&/@nvals;*)*)
(*(*critTable4=Table[{\[Delta],beta[#,critDelta[\[Delta],8]]},{\[Delta],0.0,1,0.1}]&/@ploidies;*)*)


(* ::Input:: *)
(*(*vary eta*)*)
(*critDelta[\[Delta]_?NumberQ, n_?NumberQ]:=\[Eta]/.FindRoot[myfun[1,\[Delta],n,\[Eta]]==0,{\[Eta],0.9999,0.000005,100}];*)


(* ::Input:: *)
(*(*vary theta*)*)
(*critDelta[\[Delta]_?NumberQ, n_?NumberQ]:=\[Theta]/.FindRoot[myfun[1,\[Delta],n,\[Theta]]==0,{\[Theta],0.6,0.2,100}];*)


(* ::Input:: *)
(*(*vary eta LINEAR*)*)
(*critDelta[\[Delta]_?NumberQ, n_?NumberQ]:=\[Eta]/.FindRoot[myfun[1,\[Delta],n,\[Eta]]==0,{\[Eta],0.01,0.000000000001,1}];*)


(* ::Input:: *)
(*(*Plot*)*)
(*critTablevar=Table[{\[Delta],critDelta[\[Delta],#]},{\[Delta],0.1,0.999,0.1}]&/@nvals;*)
(*critCurve=ListPlot[#&/@critTablevar,*)
(*Joined->True,*)
(*Frame->True,*)
(*Axes->False,*)
(*PlotRange->{All,All},*)
(*(*Epilog\[Rule]points,*)*)
(*PlotStyle->{{Black},{Red}, {Purple}},*)
(*FrameLabel->{"\[Mu]/\[Lambda]",var},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*PlotLegends->leg,*)
(*ImageSize->300,AspectRatio->1*)
(*]*)
(*Export[NotebookDirectory[]<>"varBeta_critcurve_" <>model<>"_"<>var<>".jpg",critCurve]*)


(* ::Input:: *)
(*(*Export[NotebookDirectory[]<>"criticalcurve_newBetafunc_betaconstant_4Ks.csv",critTable[[#]],"CSV"]&/@Range[1]*)*)
(*Export[NotebookDirectory[]<>"beta_criticalcurve_" <>model<>"_"<>var<>ToString[nvals[[#]]]<>".csv",critTablevar[[#]],"CSV","TableHeadings"->{\[Delta]/\[Lambda],\[Eta]}]&/@Range[Length[nvals]];*)
(*Export[NotebookDirectory[]<>"betafunction" <>model<>"_"<>var<>".txt",{beta[x,\[Eta]]}];*)


(* ::Input:: *)
(**)


## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##;
## ## ## #VariableBirthRates ## ## ## ##;
## ## ## ## ## ## ## ## ## ## ## ## ## ## ## ##;


(* ::Input:: *)
(*(*vary eta*)*)
(*var="\[Eta]";*)
(*model="poly";*)
(*(*delta[x_, \[Eta]_] = \[Eta];*)*)
(*Clear[\[Eta]]*)
(*delta[x_,\[Eta]_, n_]:=1-\[Eta]*(Abs[Log[x]-Log[2]]/Abs[Log[n]-Log[2]]) ;*)
(*vars=Table[\[Eta],{\[Eta],0.0001,1,0.1}];*)


(* ::Input:: *)
(*(*vary eta*)*)
(*var="\[Eta]";*)
(*model="sin";*)
(*\[Theta] =0.8;*)
(*Clear[\[Eta]]*)
(*delta[x_, \[Eta]_, n_] := (1/n)*((Sin[x*3.2-1.7]+\[Eta]*x)/6)^(\[Theta]/x)*)
(*vars=Table[\[Eta],{\[Eta],0.6,120,10}];*)
(**)


(* ::Input:: *)
(*xy=Table[{x,delta[x,#,10]},{x,1,10,0.05}]&/@vars;*)
(*deltaplot=ListPlot[#&/@xy,*)
(*Joined->True,*)
(*Frame->True,*)
(*Axes->False,*)
(*PlotRange->{All,All},*)
(*FrameLabel->{"ploidy","Death rate (\[Delta])"},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*PlotStyle->ColorData[10,"ColorList"],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*PlotLegends->vars,*)
(*ImageSize->500,*)
(*AspectRatio->1,Epilog->Text[Style[delta[x, \[Eta], n],FontSize->22],{4.5,0.4}]*)
(*]*)
(**)
(*Export[NotebookDirectory[]<>"deltaplot_"<>model<>"_"<>var<>".jpg",deltaplot]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*myfun[\[Lambda]_,\[Beta]_?NumberQ,n_?NumberQ, \[Eta]_ ?NumberQ]:=Eigenvalues[Table[If[i==j, \[Lambda](2(1-\[Beta])^i-1)-delta[j, \[Eta],n],If[2j>=i>j,\[Lambda] Binomial[j,i-j]\[Beta]^(i-j) (1-\[Beta])^(2j-i),If[j>i,\[Lambda] Binomial[j,j-i]\[Beta]^(j-i) (1-\[Beta])^i,0]]],{i,n},{j,n}],*)
(*1,Method->{"Arnoldi","Criteria"->"RealPart"}][[1]]*)


(* ::Input:: *)
(*nvals={8};*)
(*ploidies={2};*)
(*leg="K = "<>ToString[#]&/@nvals;*)
(*leg2="p = "<>ToString[#]&/@ploidies;*)
(*colorList=ColorData[37,"ColorList"];*)


(* ::Input:: *)
(*(*vary eta poly*)*)
(*critBeta[\[Beta]_?NumberQ, n_?NumberQ]:=\[Eta]/.FindRoot[myfun[1,\[Beta],n,\[Eta]]==0,{\[Eta],0.9,0,2.5}];*)


(* ::Input:: *)
(*(*vary eta sin*)*)
(*critBeta[\[Beta]_?NumberQ, n_?NumberQ]:=\[Eta]/.FindRoot[myfun[1,\[Beta],n,\[Eta]]==0,{\[Eta],0.6,0.2,1200}];*)


(* ::Input:: *)
(*(*Plot*)*)
(*critTablevar=Table[{\[Beta],critBeta[\[Beta],#]},{\[Beta],0.41,0.939,0.01}]&/@nvals;*)
(*critCurve=ListPlot[#&/@critTablevar,*)
(*Joined->True,*)
(*Frame->True,*)
(*Axes->False,*)
(*PlotRange->{All,All},*)
(*(*Epilog\[Rule]points,*)*)
(*PlotStyle->{{Black},{Black, Dashed}, {Black, DotDashed}, {Black, Dotted}},*)
(*FrameLabel->{"\[Beta]",var},*)
(*FrameStyle->Directive[Thickness[0.0035],Black],*)
(*LabelStyle->{FontSize->18,Black,FontFamily->"Arial"},*)
(*PlotLegends->leg,*)
(*ImageSize->300,AspectRatio->1*)
(*]*)
(*Export[NotebookDirectory[]<>"varDelta_critcurve_" <>model<>"_"<>var<>".jpg",critCurve]*)


(* ::InheritFromParent:: *)
(*(*Export[NotebookDirectory[]<>"criticalcurve_newBetafunc_betaconstant_4Ks.csv",critTable[[#]],"CSV"]&/@Range[1]*)*)
(*Export[NotebookDirectory[]<>"delta_criticalcurve_" <>model<>"_"<>var<>ToString[nvals[[#]]]<>".csv",critTablevar[[#]],"CSV","TableHeadings"->{\[Delta]/\[Lambda],\[Eta]}]&/@Range[Length[nvals]]*)
(*Export[NotebookDirectory[]<>"deltafunction" <>model<>"_"<>var<>".txt",{beta[x,\[Eta]]}]*)



